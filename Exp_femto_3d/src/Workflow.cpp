#include "exp_femto_3d/Workflow.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TF3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVirtualPad.h"
#ifdef EXP_FEMTO_3D_HAS_CATS
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#endif
#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  namespace {

    constexpr double kProjection1DWindow = 0.04;
    constexpr double kHbarC = 0.1973269804;
    constexpr double kPiPiLikeSignBohrRadiusFm = 387.5;
    [[maybe_unused]] constexpr double kChargedPionMassMeV = 139.57018;
    constexpr double kFiniteSourceKStarMinMeV = 0.0;
    constexpr double kFiniteSourceKStarMaxMeV = 250.0;
    constexpr int kFiniteSourceKStarBins = 250;
    constexpr double kFullR2MatrixTolerance = 1e-10;
    constexpr double kLevyArgumentTolerance = 1e-12;
    constexpr double kInvalidFullModelCFValue = 1e6;
    constexpr double kFitPenaltyValue = 1e30;
    constexpr const char *kSliceCatalogBuildPhiMappingBranch = "build_uses_symmetric_phi_range";
    constexpr const char *kSliceCatalogSplitMixedEventBranch = "split_mixed_event_by_phi";
    constexpr const char *kSliceCatalogSplitMixedEventQnBranch = "split_mixed_event_by_qn";
    constexpr const char *kQnIntegratedLabel = "qn_all";

    double ComputeKStarMeV(const double q_out, const double q_side, const double q_long) {
      const double q_inv_gev = std::sqrt(q_out * q_out + q_side * q_side + q_long * q_long);
      return 0.5 * q_inv_gev * 1000.0;
    }

    struct CoulombKernelTable {
      CoulombKernelCatalogEntry catalog_entry;
      std::vector<double> kstar_mev;
      std::vector<double> coulomb_factor;

      [[nodiscard]] double Evaluate(const double q_out, const double q_side, const double q_long) const {
        if (kstar_mev.empty() || coulomb_factor.empty() || kstar_mev.size() != coulomb_factor.size()) {
          return 1.0;
        }
        const double kstar = ComputeKStarMeV(q_out, q_side, q_long);
        if (kstar <= kstar_mev.front()) {
          return coulomb_factor.front();
        }
        if (kstar >= kstar_mev.back()) {
          return 1.0;
        }
        const auto upper = std::upper_bound(kstar_mev.begin(), kstar_mev.end(), kstar);
        const std::size_t hi = static_cast<std::size_t>(std::distance(kstar_mev.begin(), upper));
        const std::size_t lo = hi > 0 ? hi - 1 : 0;
        const double x0 = kstar_mev[lo];
        const double x1 = kstar_mev[hi];
        const double y0 = coulomb_factor[lo];
        const double y1 = coulomb_factor[hi];
        if (x1 <= x0) {
          return y0;
        }
        const double fraction = (kstar - x0) / (x1 - x0);
        return y0 + fraction * (y1 - y0);
      }
    };

    const CoulombKernelTable *g_active_coulomb_kernel = nullptr;

    class ActiveCoulombKernelGuard {
     public:
      explicit ActiveCoulombKernelGuard(const CoulombKernelTable *kernel) : previous_(g_active_coulomb_kernel) {
        g_active_coulomb_kernel = kernel;
      }
      ~ActiveCoulombKernelGuard() { g_active_coulomb_kernel = previous_; }

     private:
      const CoulombKernelTable *previous_ = nullptr;
    };

    class BatchModeGuard {
     public:
      BatchModeGuard() : previous_(gROOT->IsBatch()) { gROOT->SetBatch(kTRUE); }
      ~BatchModeGuard() { gROOT->SetBatch(previous_); }

     private:
      bool previous_ = false;
    };

    struct Levy3DPMLContext {
      TH3D *h_se_raw = nullptr;
      TH3D *h_me_raw = nullptr;
      bool use_full_model = false;
      LevyFitOptions fit_options;
      const CoulombKernelTable *coulomb_kernel = nullptr;
      double raw_same_to_mixed_integral_ratio = 1.0;
    };

    Levy3DPMLContext g_levy_3d_pml_context;

    struct PhiSliceCoordinates {
      double low = 0.0;
      double high = 0.0;
      double center = std::numeric_limits<double>::quiet_NaN();
    };

    struct MixedEventProjection {
      std::unique_ptr<TH3D> raw;
      std::unique_ptr<TH3D> norm;
    };

    struct SummaryGraphStats {
      int point_count = 0;
      double mean_y = 0.0;
      double rms_y = 0.0;
    };

    struct ReportSummaryPoint {
      double phi_center = 0.0;
      double phi_error = 0.0;
      double value = 0.0;
      double error = 0.0;
      bool valid = false;
    };

    struct EpsSummaryPoint {
      double mt_center = 0.0;
      double mt_error = 0.0;
      double value = 0.0;
      double error = 0.0;
      bool valid = false;
    };

    struct HarmonicFitResult {
      double intercept = 0.0;
      double harmonic_coefficient = 0.0;
      double intercept_variance = 0.0;
      double harmonic_variance = 0.0;
      double covariance = 0.0;
      bool success = false;
    };

    using FitResultGroupKey = std::tuple<int, int, int>;
    using GroupedFitResults = std::map<FitResultGroupKey, std::vector<LevyFitResult>>;
    using EpsSummaryMap = std::map<std::pair<int, int>, std::vector<EpsSummaryPoint>>;

    struct QnSliceSelection {
      int qn_index = -1;
      double qn_low = std::numeric_limits<double>::quiet_NaN();
      double qn_high = std::numeric_limits<double>::quiet_NaN();
      std::string qn_label = kQnIntegratedLabel;
      bool is_qn_integrated = true;
    };

    std::string FormatDouble(const double value, const int precision = 2) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(precision) << value;
      return stream.str();
    }

    // Keep near-zero values stable so directory names do not depend on negative zero formatting.
    std::string FormatDirectoryValue(const double value, const int precision = 2) {
      const double stable_value = std::abs(value) < 5.0e-7 ? 0.0 : value;
      return FormatDouble(stable_value, precision);
    }

    // Report paths use stable cent/mT labels independent of ROOT's current directory.
    std::string FormatRangeText(const double low, const double high, const int precision = 2) {
      return FormatDouble(low, precision) + "-" + FormatDouble(high, precision);
    }

    std::string BuildReportRangeDirectory(const std::string &prefix, const double low, const double high) {
      return prefix + "_" + FormatDirectoryValue(low, 2) + "-" + FormatDirectoryValue(high, 2);
    }

    std::string BuildReportCentralityDirectory(const LevyFitResult &result) {
      return BuildReportRangeDirectory("cent", result.cent_low, result.cent_high);
    }

    std::string BuildReportMtDirectory(const LevyFitResult &result) {
      return BuildReportRangeDirectory("mt", result.mt_low, result.mt_high);
    }

    std::string BuildReportQnDirectory(const LevyFitResult &result) {
      return result.is_qn_integrated ? "" : result.qn_label;
    }

    std::string BuildBaseGroupId(const RangeBin &centrality_bin, const RangeBin &mt_bin) {
      return "cent_" + FormatDirectoryValue(centrality_bin.min, 2) + "-"
             + FormatDirectoryValue(centrality_bin.max, 2) + "__mt_" + FormatDirectoryValue(mt_bin.min, 2) + "-"
             + FormatDirectoryValue(mt_bin.max, 2);
    }

    std::string BuildGroupId(const RangeBin &centrality_bin,
                             const RangeBin &mt_bin,
                             const QnSliceSelection &qn_selection) {
      std::string group_id = BuildBaseGroupId(centrality_bin, mt_bin);
      if (!qn_selection.is_qn_integrated) {
        group_id += "__" + qn_selection.qn_label;
      }
      return group_id;
    }

    std::string BuildSliceId(const std::string &group_id,
                             const bool is_phi_integrated,
                             const double display_phi_center) {
      if (is_phi_integrated) {
        return group_id + "__phi_all";
      }
      return group_id + "__phi_" + FormatDirectoryValue(display_phi_center, 2);
    }

    std::string BuildSliceDirectory(const std::string &slice_id) {
      return "slices/" + slice_id;
    }

    std::string BuildFitDirectory(const std::string &slice_id) {
      return "fits/" + slice_id;
    }

    // Keep build and fit on the same phi remapping rule so the catalog is the
    // single source of truth for later reinterpretation.
    PhiSliceCoordinates BuildPhiSliceCoordinatesFromRaw(const double raw_phi_low,
                                                        const double raw_phi_high,
                                                        const double raw_phi_center,
                                                        const bool use_symmetric_phi_range) {
      PhiSliceCoordinates coordinates{raw_phi_low, raw_phi_high, raw_phi_center};
      if (use_symmetric_phi_range && raw_phi_center > TMath::Pi() / 2.0) {
        coordinates.low -= TMath::Pi();
        coordinates.high -= TMath::Pi();
        coordinates.center -= TMath::Pi();
      }
      return coordinates;
    }

    // Integrated slices do not carry a phi center, but fit-side reinterpretation
    // still needs a stable helper for non-integrated slices.
    PhiSliceCoordinates ResolveSlicePhiCoordinates(const SliceCatalogEntry &entry, const bool use_symmetric_phi_range) {
      if (entry.is_phi_integrated) {
        return {entry.raw_phi_low, entry.raw_phi_high, std::numeric_limits<double>::quiet_NaN()};
      }
      return BuildPhiSliceCoordinatesFromRaw(
          entry.raw_phi_low, entry.raw_phi_high, entry.raw_phi_center, use_symmetric_phi_range);
    }

    bool PhiCoordinateMatches(const double lhs, const double rhs) {
      if (std::isnan(lhs) && std::isnan(rhs)) {
        return true;
      }
      return NearlyEqual(lhs, rhs);
    }

    bool MatchesStoredDisplayPhiCoordinates(const SliceCatalogEntry &entry, const PhiSliceCoordinates &coordinates) {
      return PhiCoordinateMatches(entry.display_phi_low, coordinates.low)
             && PhiCoordinateMatches(entry.display_phi_high, coordinates.high)
             && PhiCoordinateMatches(entry.display_phi_center, coordinates.center);
    }

    // Legacy SliceCatalog trees do not store the build-side phi mapping flag, so
    // infer one file-level state from the raw/display coordinate pairs.
    bool InferLegacyBuildPhiMappingState(std::vector<SliceCatalogEntry> &entries, const Logger *logger) {
      bool saw_explicit_mapped_slice = false;
      bool saw_explicit_raw_slice = false;
      for (const SliceCatalogEntry &entry : entries) {
        if (entry.is_phi_integrated) {
          continue;
        }
        const PhiSliceCoordinates raw_coordinates = ResolveSlicePhiCoordinates(entry, false);
        const PhiSliceCoordinates mapped_coordinates = ResolveSlicePhiCoordinates(entry, true);
        const bool matches_raw = MatchesStoredDisplayPhiCoordinates(entry, raw_coordinates);
        const bool matches_mapped = MatchesStoredDisplayPhiCoordinates(entry, mapped_coordinates);
        if (matches_mapped && !matches_raw) {
          saw_explicit_mapped_slice = true;
        }
        if (matches_raw && !matches_mapped) {
          saw_explicit_raw_slice = true;
        }
      }

      bool inferred_uses_symmetric_phi_range = false;
      std::string inference_reason = "defaulted to raw phi because no discriminating slice was found";
      if (saw_explicit_mapped_slice && !saw_explicit_raw_slice) {
        inferred_uses_symmetric_phi_range = true;
        inference_reason = "found remapped phi slices";
      } else if (saw_explicit_raw_slice && !saw_explicit_mapped_slice) {
        inference_reason = "found raw phi slices above pi/2";
      } else if (saw_explicit_mapped_slice && saw_explicit_raw_slice) {
        inferred_uses_symmetric_phi_range = true;
        inference_reason = "found conflicting raw/display phi pairs, keeping the mapped interpretation";
      }

      for (SliceCatalogEntry &entry : entries) {
        entry.build_uses_symmetric_phi_range = inferred_uses_symmetric_phi_range;
      }

      if (logger != nullptr) {
        logger->Warn("SliceCatalog is missing '" + std::string(kSliceCatalogBuildPhiMappingBranch)
                     + "'; inferred build phi mapping state from raw/display coordinates (" + inference_reason + ").");
      }
      return inferred_uses_symmetric_phi_range;
    }

    // Build is expected to write one global phi mapping mode per CF file, so
    // normalize inconsistent metadata before fit consumes it.
    bool NormalizeCatalogBuildPhiMappingState(std::vector<SliceCatalogEntry> &entries, const Logger *logger) {
      if (entries.empty()) {
        return false;
      }

      const bool normalized_state = entries.front().build_uses_symmetric_phi_range;
      const bool has_inconsistent_state = std::any_of(entries.begin(), entries.end(), [&](const SliceCatalogEntry &entry) {
        return entry.build_uses_symmetric_phi_range != normalized_state;
      });
      if (has_inconsistent_state && logger != nullptr) {
        logger->Warn("SliceCatalog contains inconsistent build phi mapping metadata; using the first entry as the"
                     " file-level truth.");
      }
      for (SliceCatalogEntry &entry : entries) {
        entry.build_uses_symmetric_phi_range = normalized_state;
      }
      return normalized_state;
    }

    std::string ResolvePath(const std::string &directory, const std::string &file_name) {
      const std::filesystem::path base(directory);
      const std::filesystem::path leaf(file_name);
      if (leaf.is_absolute()) {
        return leaf.string();
      }
      return (base / leaf).string();
    }

    void EnsureDirectoryExists(const std::string &directory) {
      std::filesystem::create_directories(directory);
    }

    std::unique_ptr<TFile> OpenRootFile(const std::string &path, const char *mode) {
      auto file = std::make_unique<TFile>(path.c_str(), mode);
      if (!file || file->IsZombie()) {
        throw std::runtime_error("Cannot open ROOT file: " + path + " with mode " + mode);
      }
      return file;
    }

    void CreateOrResetRootFile(const std::string &path) {
      auto file = OpenRootFile(path, "RECREATE");
      file->Close();
    }

    // Resolve configured file paths before comparing them so report output cannot overwrite source artifacts.
    bool PathsReferToSameLocation(const std::string &left, const std::string &right) {
      try {
        return std::filesystem::weakly_canonical(left) == std::filesystem::weakly_canonical(right);
      } catch (const std::filesystem::filesystem_error &) {
        return std::filesystem::absolute(left).lexically_normal()
               == std::filesystem::absolute(right).lexically_normal();
      }
    }

    TDirectory *GetOrCreateDirectory(TDirectory &parent, const std::string &name) {
      if (auto *existing = parent.GetDirectory(name.c_str())) {
        return existing;
      }

      TDirectory *created = parent.mkdir(name.c_str());
      if (created == nullptr) {
        throw std::runtime_error("Failed to create ROOT directory: " + name);
      }
      return created;
    }

    TDirectory *GetOrCreateDirectoryPath(TDirectory &parent, const std::string &path) {
      TDirectory *current = &parent;
      std::size_t begin = 0;
      while (begin < path.size()) {
        const std::size_t separator = path.find('/', begin);
        const std::string component =
            path.substr(begin, separator == std::string::npos ? std::string::npos : separator - begin);
        if (!component.empty()) {
          current = GetOrCreateDirectory(*current, component);
        }
        if (separator == std::string::npos) {
          break;
        }
        begin = separator + 1U;
      }
      return current;
    }

    // Report canvases keep the existing second-order event-plane convention visible.
    std::string BuildRelativePhiAxisTitle() {
      return "#phi_{pair} - #Psi_{2} [rad]";
    }

    std::string FormatSummaryStat(const double value) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(3) << value;
      return stream.str();
    }

    // Radius-index mapping keeps overview and epsilon code aligned with the fit summary schema.
    double FitResultRadiusValue(const LevyFitResult &result, const std::size_t radius_index) {
      switch (radius_index) {
        case 0:
          return result.rout2;
        case 1:
          return result.rside2;
        case 2:
          return result.rlong2;
        case 3:
          return result.routside2;
        case 4:
          return result.routlong2;
        case 5:
          return result.rsidelong2;
        default:
          return std::numeric_limits<double>::quiet_NaN();
      }
    }

    double FitResultRadiusError(const LevyFitResult &result, const std::size_t radius_index) {
      switch (radius_index) {
        case 0:
          return result.rout2_err;
        case 1:
          return result.rside2_err;
        case 2:
          return result.rlong2_err;
        case 3:
          return result.routside2_err;
        case 4:
          return result.routlong2_err;
        case 5:
          return result.rsidelong2_err;
        default:
          return std::numeric_limits<double>::quiet_NaN();
      }
    }

    ReportSummaryPoint MakeFitResultSummaryPoint(const LevyFitResult &result,
                                                 const double value,
                                                 const double error) {
      ReportSummaryPoint point;
      point.phi_center = result.phi;
      point.phi_error = 0.0;
      point.value = value;
      point.error = std::isfinite(error) && error > 0.0 ? error : 0.0;
      point.valid = std::isfinite(point.phi_center) && std::isfinite(point.value);
      return point;
    }

    // Overview and epsilon summaries are built only from phi-differential slices.
    GroupedFitResults GroupPhiDifferentialResultsByCentMt(const std::vector<LevyFitResult> &results) {
      GroupedFitResults grouped_results;
      for (const LevyFitResult &result : results) {
        if (result.is_phi_integrated) {
          continue;
        }
        grouped_results[{result.centrality_index, result.mt_index, result.qn_index}].push_back(result);
      }
      for (auto &[key, group_results] : grouped_results) {
        (void)key;
        std::sort(group_results.begin(), group_results.end(), [](const LevyFitResult &left,
                                                                 const LevyFitResult &right) {
          return left.phi < right.phi;
        });
      }
      return grouped_results;
    }

    std::vector<ReportSummaryPoint> BuildRadiusSummaryPoints(const std::vector<LevyFitResult> &results,
                                                             const std::size_t radius_index) {
      std::vector<ReportSummaryPoint> points;
      points.reserve(results.size());
      for (const LevyFitResult &result : results) {
        points.push_back(MakeFitResultSummaryPoint(
            result, FitResultRadiusValue(result, radius_index), FitResultRadiusError(result, radius_index)));
      }
      return points;
    }

    std::vector<ReportSummaryPoint> BuildAlphaSummaryPoints(const std::vector<LevyFitResult> &results) {
      std::vector<ReportSummaryPoint> points;
      points.reserve(results.size());
      for (const LevyFitResult &result : results) {
        points.push_back(MakeFitResultSummaryPoint(result, result.alpha, result.alpha_err));
      }
      return points;
    }

    // The eps canvas keeps a compact statistics box for quick visual QA.
    SummaryGraphStats ComputeEpsGraphStats(const std::vector<EpsSummaryPoint> &points) {
      SummaryGraphStats stats;
      double sum_y = 0.0;
      double sum_y2 = 0.0;
      for (const EpsSummaryPoint &point : points) {
        if (!point.valid || !std::isfinite(point.value)) {
          continue;
        }
        ++stats.point_count;
        sum_y += point.value;
        sum_y2 += point.value * point.value;
      }
      if (stats.point_count == 0) {
        return stats;
      }
      const double normalization = static_cast<double>(stats.point_count);
      stats.mean_y = sum_y / normalization;
      stats.rms_y = std::sqrt(std::max(0.0, sum_y2 / normalization - stats.mean_y * stats.mean_y));
      return stats;
    }

    void ApplyReportGraphStyle(TGraphErrors &graph) {
      constexpr int kReportGraphColor = 602;
      constexpr int kReportMarkerStyle = 20;
      graph.SetMarkerStyle(kReportMarkerStyle);
      graph.SetMarkerSize(1.1);
      graph.SetMarkerColor(kReportGraphColor);
      graph.SetLineColor(kReportGraphColor);
      graph.SetLineWidth(2);
      graph.SetFillStyle(0);
      graph.SetDrawOption("ALPE1");
    }

    void ApplyOverviewGraphStyle(TGraphErrors &graph) {
      ApplyReportGraphStyle(graph);
      graph.GetXaxis()->SetTitleSize(0.055);
      graph.GetXaxis()->SetLabelSize(0.045);
      graph.GetYaxis()->SetTitleSize(0.060);
      graph.GetYaxis()->SetLabelSize(0.045);
      graph.GetYaxis()->SetTitleOffset(1.15);
    }

    TPaveText *CreateSummaryInfoBox(const std::string &graph_name, const SummaryGraphStats &stats) {
      auto *info_box = new TPaveText(0.63, 0.72, 0.89, 0.89, "NDC");
      info_box->SetName((graph_name + "_summary_box").c_str());
      info_box->SetBorderSize(1);
      info_box->SetFillColor(0);
      info_box->SetFillStyle(1001);
      info_box->SetShadowColor(0);
      info_box->SetTextAlign(12);
      info_box->SetTextFont(42);
      info_box->AddText(("Points: " + std::to_string(stats.point_count)).c_str());
      info_box->AddText(("Mean(Y): " + FormatSummaryStat(stats.mean_y)).c_str());
      info_box->AddText(("RMS(Y): " + FormatSummaryStat(stats.rms_y)).c_str());
      return info_box;
    }

    TPaveText *CreateSourceParameterOverviewInfoBox(const LevyFitResult &representative_result) {
      auto *info_box = new TPaveText(0.08, 0.08, 0.92, 0.92, "NDC");
      info_box->SetName("source_parameters_overview_info_box");
      info_box->SetBorderSize(0);
      info_box->SetFillColor(0);
      info_box->SetFillStyle(0);
      info_box->SetShadowColor(0);
      info_box->SetTextAlign(12);
      info_box->SetTextFont(42);
      info_box->SetTextSize(0.055);
      info_box->AddText("Exp femto 3D fit source parameters");
      info_box->AddText(("centrality: "
                         + FormatRangeText(representative_result.cent_low, representative_result.cent_high))
                            .c_str());
      info_box->AddText(("m_{T}: " + FormatRangeText(representative_result.mt_low, representative_result.mt_high)
                         + " GeV/c^{2}")
                            .c_str());
      if (representative_result.is_qn_integrated) {
        info_box->AddText("qn: all");
      } else {
        info_box->AddText(("qn: " + representative_result.qn_label + " "
                           + FormatRangeText(representative_result.qn_low, representative_result.qn_high))
                              .c_str());
      }
      info_box->AddText(BuildRelativePhiAxisTitle().c_str());
      return info_box;
    }

    bool HasValidSummaryPoints(const std::vector<ReportSummaryPoint> &points) {
      return std::any_of(points.begin(), points.end(), [](const ReportSummaryPoint &point) {
        return point.valid && std::isfinite(point.phi_center) && std::isfinite(point.value);
      });
    }

    std::unique_ptr<TGraphErrors> MakePhiSummaryGraph(const std::vector<ReportSummaryPoint> &points,
                                                      const std::string &graph_name,
                                                      const std::string &graph_title,
                                                      const std::string &x_axis_title,
                                                      const std::string &y_axis_title) {
      const auto valid_count =
          static_cast<int>(std::count_if(points.begin(), points.end(), [](const ReportSummaryPoint &point) {
            return point.valid && std::isfinite(point.phi_center) && std::isfinite(point.value);
          }));
      if (valid_count == 0) {
        return nullptr;
      }

      auto graph = std::make_unique<TGraphErrors>(valid_count);
      graph->SetName(graph_name.c_str());
      graph->SetTitle(graph_title.c_str());
      graph->GetXaxis()->SetTitle(x_axis_title.c_str());
      graph->GetYaxis()->SetTitle(y_axis_title.c_str());
      ApplyReportGraphStyle(*graph);

      int point_index = 0;
      for (const ReportSummaryPoint &point : points) {
        if (!point.valid || !std::isfinite(point.phi_center) || !std::isfinite(point.value)) {
          continue;
        }
        graph->SetPoint(point_index, point.phi_center, point.value);
        graph->SetPointError(point_index,
                             std::isfinite(point.phi_error) ? point.phi_error : 0.0,
                             std::isfinite(point.error) ? point.error : 0.0);
        ++point_index;
      }
      return graph;
    }

    // Match the Eventgen epsilon extraction: fit Rside2(phi) to a second-harmonic form.
    HarmonicFitResult FitSideRadiusSecondHarmonic(const std::vector<ReportSummaryPoint> &rside_points) {
      constexpr double kMinimumDeterminant = 1.0e-20;
      int used_points = 0;
      bool use_point_errors = true;
      for (const ReportSummaryPoint &point : rside_points) {
        if (!point.valid || !std::isfinite(point.phi_center) || !std::isfinite(point.value)) {
          continue;
        }
        ++used_points;
        if (!std::isfinite(point.error) || point.error <= 0.0) {
          use_point_errors = false;
        }
      }
      if (used_points < 2) {
        return {};
      }

      double sum_w = 0.0;
      double sum_wc = 0.0;
      double sum_wcc = 0.0;
      double sum_wy = 0.0;
      double sum_wcy = 0.0;
      for (const ReportSummaryPoint &point : rside_points) {
        if (!point.valid || !std::isfinite(point.phi_center) || !std::isfinite(point.value)) {
          continue;
        }
        const double cos2phi = std::cos(2.0 * point.phi_center);
        const double weight = use_point_errors ? 1.0 / (point.error * point.error) : 1.0;
        sum_w += weight;
        sum_wc += weight * cos2phi;
        sum_wcc += weight * cos2phi * cos2phi;
        sum_wy += weight * point.value;
        sum_wcy += weight * cos2phi * point.value;
      }

      const double determinant = sum_w * sum_wcc - sum_wc * sum_wc;
      if (!std::isfinite(determinant) || std::abs(determinant) <= kMinimumDeterminant) {
        return {};
      }

      HarmonicFitResult fit;
      fit.intercept = (sum_wy * sum_wcc - sum_wcy * sum_wc) / determinant;
      fit.harmonic_coefficient = (sum_w * sum_wcy - sum_wc * sum_wy) / determinant;
      fit.intercept_variance = sum_wcc / determinant;
      fit.harmonic_variance = sum_w / determinant;
      fit.covariance = -sum_wc / determinant;

      if (!use_point_errors && used_points > 2) {
        double residual_sum2 = 0.0;
        for (const ReportSummaryPoint &point : rside_points) {
          if (!point.valid || !std::isfinite(point.phi_center) || !std::isfinite(point.value)) {
            continue;
          }
          const double residual =
              point.value - fit.intercept - fit.harmonic_coefficient * std::cos(2.0 * point.phi_center);
          residual_sum2 += residual * residual;
        }
        const double residual_variance = residual_sum2 / static_cast<double>(used_points - 2);
        fit.intercept_variance *= residual_variance;
        fit.harmonic_variance *= residual_variance;
        fit.covariance *= residual_variance;
      }

      fit.success = std::isfinite(fit.intercept) && std::isfinite(fit.harmonic_coefficient)
                    && std::isfinite(fit.intercept_variance) && std::isfinite(fit.harmonic_variance)
                    && std::isfinite(fit.covariance);
      return fit;
    }

    // Convert the Rside2 harmonic fit into the plotted epsilon_f(mT) point.
    std::optional<EpsSummaryPoint> ComputeEpsFromRsideSummaryPoints(
        const std::vector<ReportSummaryPoint> &rside_points,
        const LevyFitResult &representative_result) {
      constexpr double kMinimumInterceptMagnitude = 1.0e-12;
      const HarmonicFitResult fit = FitSideRadiusSecondHarmonic(rside_points);
      if (!fit.success || std::abs(fit.intercept) <= kMinimumInterceptMagnitude) {
        return std::nullopt;
      }

      EpsSummaryPoint point;
      point.mt_center = 0.5 * (representative_result.mt_low + representative_result.mt_high);
      point.mt_error = 0.5 * (representative_result.mt_high - representative_result.mt_low);
      point.value = fit.harmonic_coefficient / fit.intercept;
      const double ratio_variance =
          (fit.harmonic_coefficient * fit.harmonic_coefficient * fit.intercept_variance)
              / (fit.intercept * fit.intercept * fit.intercept * fit.intercept)
          + fit.harmonic_variance / (fit.intercept * fit.intercept)
          - (2.0 * fit.harmonic_coefficient * fit.covariance) / (fit.intercept * fit.intercept * fit.intercept);
      point.error = std::isfinite(ratio_variance) ? std::sqrt(std::max(0.0, ratio_variance)) : 0.0;
      point.valid = std::isfinite(point.mt_center) && std::isfinite(point.value) && std::isfinite(point.error);
      return point.valid ? std::optional<EpsSummaryPoint>(point) : std::nullopt;
    }

    // Overview trend lines reuse the same constant/cosine/sine forms as R2 summary fits.
    std::unique_ptr<TF1> BuildOverviewFitFunction(const std::string &name,
                                                  const int parameter_index,
                                                  const double phi_fit_min,
                                                  const double phi_fit_max,
                                                  const double initial_value) {
      const bool alpha_panel = parameter_index < 0;
      const bool sine_panel = parameter_index == 3 || parameter_index == 5;
      const char *formula = alpha_panel ? "[0]" : (sine_panel ? "[0]+2.0*[1]*sin(2.0*x)"
                                                              : "[0]+2.0*[1]*cos(2.0*x)");
      auto fit_function = std::make_unique<TF1>(name.c_str(), formula, phi_fit_min, phi_fit_max);
      fit_function->SetParameter(0, std::isfinite(initial_value) ? initial_value : 0.0);
      if (!alpha_panel) {
        fit_function->SetParameter(1, 0.0);
      }
      fit_function->SetLineColor(kRed + 1);
      fit_function->SetLineWidth(3);
      return fit_function;
    }

    void DrawOverviewGraph(TCanvas &canvas,
                           const int pad_number,
                           std::unique_ptr<TGraphErrors> graph,
                           const int parameter_index,
                           const bool uses_mapped_phi_range,
                           std::vector<std::unique_ptr<TGraphErrors>> &owned_graphs,
                           std::vector<std::unique_ptr<TF1>> &owned_fits) {
      if (graph == nullptr) {
        return;
      }

      TVirtualPad *pad = canvas.cd(pad_number);
      if (pad == nullptr) {
        return;
      }

      pad->SetTicks(1, 1);
      pad->SetLeftMargin(0.16);
      pad->SetRightMargin(0.06);
      pad->SetBottomMargin(0.14);
      pad->SetTopMargin(0.08);
      ApplyOverviewGraphStyle(*graph);
      graph->SetTitle("");
      graph->Draw("ALPE1");

      // Overview canvases carry the same harmonic trend line convention as the stored R2 summary fits.
      double x = 0.0;
      double y = 0.0;
      graph->GetPoint(0, x, y);
      const double phi_fit_min = uses_mapped_phi_range ? -TMath::Pi() / 2.0 : 0.0;
      const double phi_fit_max = uses_mapped_phi_range ? TMath::Pi() / 2.0 : TMath::Pi();
      const bool constant_panel = parameter_index < 0;
      const int minimum_points = constant_panel ? 1 : 2;
      if (graph->GetN() >= minimum_points) {
        auto fit_function = BuildOverviewFitFunction(
            std::string(graph->GetName()) + "_overview_fit", parameter_index, phi_fit_min, phi_fit_max, y);
        graph->Fit(fit_function.get(), "QN");
        fit_function->Draw("L SAME");
        owned_fits.push_back(std::move(fit_function));
      }
      owned_graphs.push_back(std::move(graph));
    }

    double IntegralVisibleRange(TH3D *hist, const bool use_width = false) {
      if (hist == nullptr) {
        return 0.0;
      }
      return hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1, hist->GetNbinsZ(), use_width ? "width" : "");
    }

    std::pair<int, int> GetAxisRangeForWindow(const TAxis *axis, const double q_max) {
      const int first_bin = axis->FindBin(-q_max + 1.0e-9);
      const int last_bin = axis->FindBin(q_max - 1.0e-9);
      return {std::max(first_bin, 1), std::min(last_bin, axis->GetNbins())};
    }

    TH1D *BuildProjectionXWithinWindow(TH3D *hist, const std::string &name, const double q_max) {
      auto [y_min, y_max] = GetAxisRangeForWindow(hist->GetYaxis(), q_max);
      auto [z_min, z_max] = GetAxisRangeForWindow(hist->GetZaxis(), q_max);
      auto *projection = hist->ProjectionX(name.c_str(), y_min, y_max, z_min, z_max);
      const int n_window_bins = std::max(0, y_max - y_min + 1) * std::max(0, z_max - z_min + 1);
      if (n_window_bins > 0) {
        projection->Scale(1.0 / static_cast<double>(n_window_bins));
      }
      return projection;
    }

    TH1D *BuildProjectionYWithinWindow(TH3D *hist, const std::string &name, const double q_max) {
      auto [x_min, x_max] = GetAxisRangeForWindow(hist->GetXaxis(), q_max);
      auto [z_min, z_max] = GetAxisRangeForWindow(hist->GetZaxis(), q_max);
      auto *projection = hist->ProjectionY(name.c_str(), x_min, x_max, z_min, z_max);
      const int n_window_bins = std::max(0, x_max - x_min + 1) * std::max(0, z_max - z_min + 1);
      if (n_window_bins > 0) {
        projection->Scale(1.0 / static_cast<double>(n_window_bins));
      }
      return projection;
    }

    TH1D *BuildProjectionZWithinWindow(TH3D *hist, const std::string &name, const double q_max) {
      auto [x_min, x_max] = GetAxisRangeForWindow(hist->GetXaxis(), q_max);
      auto [y_min, y_max] = GetAxisRangeForWindow(hist->GetYaxis(), q_max);
      auto *projection = hist->ProjectionZ(name.c_str(), x_min, x_max, y_min, y_max);
      const int n_window_bins = std::max(0, x_max - x_min + 1) * std::max(0, y_max - y_min + 1);
      if (n_window_bins > 0) {
        projection->Scale(1.0 / static_cast<double>(n_window_bins));
      }
      return projection;
    }

    void Write1DProjections(TH3D *histogram,
                            TDirectory &directory,
                            const std::string &base_name,
                            const std::string &y_title,
                            const bool use_window = true) {
      TH1D *projection_x = nullptr;
      TH1D *projection_y = nullptr;
      TH1D *projection_z = nullptr;

      if (use_window) {
        projection_x = BuildProjectionXWithinWindow(histogram, base_name + "_ProjX", kProjection1DWindow);
        projection_y = BuildProjectionYWithinWindow(histogram, base_name + "_ProjY", kProjection1DWindow);
        projection_z = BuildProjectionZWithinWindow(histogram, base_name + "_ProjZ", kProjection1DWindow);
      } else {
        projection_x = histogram->ProjectionX(
            (base_name + "_ProjX").c_str(), 1, histogram->GetNbinsY(), 1, histogram->GetNbinsZ());
        projection_y = histogram->ProjectionY(
            (base_name + "_ProjY").c_str(), 1, histogram->GetNbinsX(), 1, histogram->GetNbinsZ());
        projection_z = histogram->ProjectionZ(
            (base_name + "_ProjZ").c_str(), 1, histogram->GetNbinsX(), 1, histogram->GetNbinsY());
      }

      projection_x->SetName((base_name + "_ProjX").c_str());
      projection_y->SetName((base_name + "_ProjY").c_str());
      projection_z->SetName((base_name + "_ProjZ").c_str());

      projection_x->GetXaxis()->SetTitle("q_{out} (GeV/c)");
      projection_y->GetXaxis()->SetTitle("q_{side} (GeV/c)");
      projection_z->GetXaxis()->SetTitle("q_{long} (GeV/c)");
      projection_x->GetYaxis()->SetTitle(y_title.c_str());
      projection_y->GetYaxis()->SetTitle(y_title.c_str());
      projection_z->GetYaxis()->SetTitle(y_title.c_str());

      directory.WriteObject(projection_x, projection_x->GetName());
      directory.WriteObject(projection_y, projection_y->GetName());
      directory.WriteObject(projection_z, projection_z->GetName());

      delete projection_x;
      delete projection_y;
      delete projection_z;
    }

    std::vector<QnSliceSelection> BuildQnSliceSelections(const ApplicationConfig &config, const TAxis &qn_axis) {
      std::vector<QnSliceSelection> selections;
      QnSliceSelection integrated_selection;
      integrated_selection.qn_index = -1;
      integrated_selection.qn_low = qn_axis.GetBinLowEdge(1);
      integrated_selection.qn_high = qn_axis.GetBinUpEdge(qn_axis.GetNbins());
      integrated_selection.qn_label = kQnIntegratedLabel;
      integrated_selection.is_qn_integrated = true;
      selections.push_back(integrated_selection);

      if (!config.build.split_same_event_by_qn) {
        return selections;
      }

      for (std::size_t index = 0; index < config.qn_bins.size(); ++index) {
        const RangeBin &qn_bin = config.qn_bins[index];
        QnSliceSelection selection;
        selection.qn_index = static_cast<int>(index);
        selection.qn_low = qn_bin.min;
        selection.qn_high = qn_bin.max;
        selection.qn_label = qn_bin.label;
        selection.is_qn_integrated = false;
        selections.push_back(selection);
      }
      return selections;
    }

    void ResetAxisVisibleRange(TAxis &axis) {
      axis.SetRange(1, axis.GetNbins());
    }

    // qn bins are configured as half-open physical ranges and mapped to exact THnSparse axis bins.
    void ApplyQnSelection(THnSparseF &sparse, const QnSliceSelection &selection) {
      TAxis *axis = sparse.GetAxis(5);
      if (axis == nullptr) {
        throw std::runtime_error("Input THnSparse is missing qn axis 5.");
      }
      if (selection.is_qn_integrated) {
        ResetAxisVisibleRange(*axis);
        return;
      }
      if (selection.qn_high <= axis->GetBinLowEdge(1) || selection.qn_low >= axis->GetBinUpEdge(axis->GetNbins())) {
        throw std::runtime_error("Configured qn bin " + selection.qn_label + " does not overlap THnSparse qn axis.");
      }

      const double low_inside = std::nextafter(selection.qn_low, std::numeric_limits<double>::infinity());
      const double high_inside = std::nextafter(selection.qn_high, -std::numeric_limits<double>::infinity());
      int first_bin = axis->FindFixBin(low_inside);
      int last_bin = axis->FindFixBin(high_inside);
      first_bin = std::max(1, std::min(first_bin, axis->GetNbins()));
      last_bin = std::max(1, std::min(last_bin, axis->GetNbins()));
      if (last_bin < first_bin) {
        throw std::runtime_error("Configured qn bin " + selection.qn_label + " does not overlap THnSparse qn axis.");
      }
      axis->SetRange(first_bin, last_bin);
    }

    SliceCatalogEntry MakeSliceCatalogEntry(const RangeBin &centrality_bin,
                                            const RangeBin &mt_bin,
                                            const QnSliceSelection &qn_selection,
                                            const int centrality_index,
                                            const int mt_index,
                                            const int phi_index,
                                            const double raw_phi_low,
                                            const double raw_phi_high,
                                            const double raw_phi_center,
                                            const double display_phi_low,
                                            const double display_phi_high,
                                            const double display_phi_center,
                                            const bool build_uses_symmetric_phi_range,
                                            const bool split_mixed_event_by_phi,
                                            const bool split_mixed_event_by_qn,
                                            const bool is_phi_integrated) {
      SliceCatalogEntry entry;
      entry.group_id = BuildGroupId(centrality_bin, mt_bin, qn_selection);
      entry.slice_id = BuildSliceId(entry.group_id, is_phi_integrated, display_phi_center);
      entry.slice_directory = BuildSliceDirectory(entry.slice_id);
      entry.se_object_path = entry.slice_directory + "/SE_raw3d";
      entry.me_object_path = entry.slice_directory + "/ME_raw3d";
      entry.cf_object_path = entry.slice_directory + "/CF3D";
      entry.projection_x_path = entry.slice_directory + "/CF3D_ProjX";
      entry.projection_y_path = entry.slice_directory + "/CF3D_ProjY";
      entry.projection_z_path = entry.slice_directory + "/CF3D_ProjZ";
      entry.centrality_index = centrality_index;
      entry.mt_index = mt_index;
      entry.qn_index = qn_selection.qn_index;
      entry.phi_index = phi_index;
      entry.cent_low = centrality_bin.min;
      entry.cent_high = centrality_bin.max;
      entry.mt_low = mt_bin.min;
      entry.mt_high = mt_bin.max;
      entry.qn_low = qn_selection.qn_low;
      entry.qn_high = qn_selection.qn_high;
      entry.qn_label = qn_selection.qn_label;
      entry.raw_phi_low = raw_phi_low;
      entry.raw_phi_high = raw_phi_high;
      entry.raw_phi_center = raw_phi_center;
      entry.display_phi_low = display_phi_low;
      entry.display_phi_high = display_phi_high;
      entry.display_phi_center = display_phi_center;
      entry.build_uses_symmetric_phi_range = build_uses_symmetric_phi_range;
      entry.split_mixed_event_by_phi = split_mixed_event_by_phi;
      entry.split_mixed_event_by_qn = split_mixed_event_by_qn;
      entry.is_qn_integrated = qn_selection.is_qn_integrated;
      entry.is_phi_integrated = is_phi_integrated;
      return entry;
    }

    void WriteSliceCatalogTree(TFile &output_file, const std::vector<SliceCatalogEntry> &entries) {
      auto *meta_directory = GetOrCreateDirectoryPath(output_file, "meta");
      meta_directory->cd();

      auto tree = std::make_unique<TTree>("SliceCatalog", "SliceCatalog");
      std::string slice_id;
      std::string group_id;
      std::string slice_directory;
      std::string se_object_path;
      std::string me_object_path;
      std::string cf_object_path;
      std::string projection_x_path;
      std::string projection_y_path;
      std::string projection_z_path;
      int centrality_index = -1;
      int mt_index = -1;
      int qn_index = -1;
      int phi_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      double qn_low = std::numeric_limits<double>::quiet_NaN();
      double qn_high = std::numeric_limits<double>::quiet_NaN();
      std::string qn_label;
      double raw_phi_low = 0.0;
      double raw_phi_high = 0.0;
      double raw_phi_center = 0.0;
      double display_phi_low = 0.0;
      double display_phi_high = 0.0;
      double display_phi_center = 0.0;
      int build_uses_symmetric_phi_range = 0;
      int split_mixed_event_by_phi = 0;
      int split_mixed_event_by_qn = 0;
      int is_qn_integrated = 1;
      int is_phi_integrated = 0;

      tree->Branch("slice_id", &slice_id);
      tree->Branch("group_id", &group_id);
      tree->Branch("slice_directory", &slice_directory);
      tree->Branch("se_object_path", &se_object_path);
      tree->Branch("me_object_path", &me_object_path);
      tree->Branch("cf_object_path", &cf_object_path);
      tree->Branch("projection_x_path", &projection_x_path);
      tree->Branch("projection_y_path", &projection_y_path);
      tree->Branch("projection_z_path", &projection_z_path);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("qn_index", &qn_index);
      tree->Branch("phi_index", &phi_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("qn_low", &qn_low);
      tree->Branch("qn_high", &qn_high);
      tree->Branch("qn_label", &qn_label);
      tree->Branch("raw_phi_low", &raw_phi_low);
      tree->Branch("raw_phi_high", &raw_phi_high);
      tree->Branch("raw_phi_center", &raw_phi_center);
      tree->Branch("display_phi_low", &display_phi_low);
      tree->Branch("display_phi_high", &display_phi_high);
      tree->Branch("display_phi_center", &display_phi_center);
      tree->Branch(kSliceCatalogBuildPhiMappingBranch, &build_uses_symmetric_phi_range);
      tree->Branch(kSliceCatalogSplitMixedEventBranch, &split_mixed_event_by_phi);
      tree->Branch(kSliceCatalogSplitMixedEventQnBranch, &split_mixed_event_by_qn);
      tree->Branch("is_qn_integrated", &is_qn_integrated);
      tree->Branch("is_phi_integrated", &is_phi_integrated);

      for (const SliceCatalogEntry &entry : entries) {
        slice_id = entry.slice_id;
        group_id = entry.group_id;
        slice_directory = entry.slice_directory;
        se_object_path = entry.se_object_path;
        me_object_path = entry.me_object_path;
        cf_object_path = entry.cf_object_path;
        projection_x_path = entry.projection_x_path;
        projection_y_path = entry.projection_y_path;
        projection_z_path = entry.projection_z_path;
        centrality_index = entry.centrality_index;
        mt_index = entry.mt_index;
        qn_index = entry.qn_index;
        phi_index = entry.phi_index;
        cent_low = entry.cent_low;
        cent_high = entry.cent_high;
        mt_low = entry.mt_low;
        mt_high = entry.mt_high;
        qn_low = entry.qn_low;
        qn_high = entry.qn_high;
        qn_label = entry.qn_label;
        raw_phi_low = entry.raw_phi_low;
        raw_phi_high = entry.raw_phi_high;
        raw_phi_center = entry.raw_phi_center;
        display_phi_low = entry.display_phi_low;
        display_phi_high = entry.display_phi_high;
        display_phi_center = entry.display_phi_center;
        build_uses_symmetric_phi_range = entry.build_uses_symmetric_phi_range ? 1 : 0;
        split_mixed_event_by_phi = entry.split_mixed_event_by_phi ? 1 : 0;
        split_mixed_event_by_qn = entry.split_mixed_event_by_qn ? 1 : 0;
        is_qn_integrated = entry.is_qn_integrated ? 1 : 0;
        is_phi_integrated = entry.is_phi_integrated ? 1 : 0;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
      output_file.cd();
    }

    std::vector<SliceCatalogEntry> ReadSliceCatalogTree(TFile &input_file, const Logger *logger = nullptr) {
      auto *tree = dynamic_cast<TTree *>(input_file.Get("meta/SliceCatalog"));
      if (tree == nullptr) {
        throw std::runtime_error("Missing meta/SliceCatalog in ROOT file.");
      }

      TTreeReader reader(tree);
      TTreeReaderValue<std::string> slice_id(reader, "slice_id");
      TTreeReaderValue<std::string> group_id(reader, "group_id");
      TTreeReaderValue<std::string> slice_directory(reader, "slice_directory");
      TTreeReaderValue<std::string> se_object_path(reader, "se_object_path");
      TTreeReaderValue<std::string> me_object_path(reader, "me_object_path");
      TTreeReaderValue<std::string> cf_object_path(reader, "cf_object_path");
      TTreeReaderValue<std::string> projection_x_path(reader, "projection_x_path");
      TTreeReaderValue<std::string> projection_y_path(reader, "projection_y_path");
      TTreeReaderValue<std::string> projection_z_path(reader, "projection_z_path");
      TTreeReaderValue<int> centrality_index(reader, "centrality_index");
      TTreeReaderValue<int> mt_index(reader, "mt_index");
      std::unique_ptr<TTreeReaderValue<int>> qn_index;
      if (tree->GetBranch("qn_index") != nullptr) {
        qn_index = std::make_unique<TTreeReaderValue<int>>(reader, "qn_index");
      }
      TTreeReaderValue<int> phi_index(reader, "phi_index");
      TTreeReaderValue<double> cent_low(reader, "cent_low");
      TTreeReaderValue<double> cent_high(reader, "cent_high");
      TTreeReaderValue<double> mt_low(reader, "mt_low");
      TTreeReaderValue<double> mt_high(reader, "mt_high");
      std::unique_ptr<TTreeReaderValue<double>> qn_low;
      if (tree->GetBranch("qn_low") != nullptr) {
        qn_low = std::make_unique<TTreeReaderValue<double>>(reader, "qn_low");
      }
      std::unique_ptr<TTreeReaderValue<double>> qn_high;
      if (tree->GetBranch("qn_high") != nullptr) {
        qn_high = std::make_unique<TTreeReaderValue<double>>(reader, "qn_high");
      }
      std::unique_ptr<TTreeReaderValue<std::string>> qn_label;
      if (tree->GetBranch("qn_label") != nullptr) {
        qn_label = std::make_unique<TTreeReaderValue<std::string>>(reader, "qn_label");
      }
      TTreeReaderValue<double> raw_phi_low(reader, "raw_phi_low");
      TTreeReaderValue<double> raw_phi_high(reader, "raw_phi_high");
      TTreeReaderValue<double> raw_phi_center(reader, "raw_phi_center");
      TTreeReaderValue<double> display_phi_low(reader, "display_phi_low");
      TTreeReaderValue<double> display_phi_high(reader, "display_phi_high");
      TTreeReaderValue<double> display_phi_center(reader, "display_phi_center");
      std::unique_ptr<TTreeReaderValue<int>> build_uses_symmetric_phi_range;
      if (tree->GetBranch(kSliceCatalogBuildPhiMappingBranch) != nullptr) {
        build_uses_symmetric_phi_range =
            std::make_unique<TTreeReaderValue<int>>(reader, kSliceCatalogBuildPhiMappingBranch);
      }
      std::unique_ptr<TTreeReaderValue<int>> split_mixed_event_by_phi;
      if (tree->GetBranch(kSliceCatalogSplitMixedEventBranch) != nullptr) {
        split_mixed_event_by_phi =
            std::make_unique<TTreeReaderValue<int>>(reader, kSliceCatalogSplitMixedEventBranch);
      }
      std::unique_ptr<TTreeReaderValue<int>> split_mixed_event_by_qn;
      if (tree->GetBranch(kSliceCatalogSplitMixedEventQnBranch) != nullptr) {
        split_mixed_event_by_qn =
            std::make_unique<TTreeReaderValue<int>>(reader, kSliceCatalogSplitMixedEventQnBranch);
      }
      std::unique_ptr<TTreeReaderValue<int>> is_qn_integrated;
      if (tree->GetBranch("is_qn_integrated") != nullptr) {
        is_qn_integrated = std::make_unique<TTreeReaderValue<int>>(reader, "is_qn_integrated");
      }
      TTreeReaderValue<int> is_phi_integrated(reader, "is_phi_integrated");

      std::vector<SliceCatalogEntry> entries;
      while (reader.Next()) {
        SliceCatalogEntry entry;
        entry.slice_id = *slice_id;
        entry.group_id = *group_id;
        entry.slice_directory = *slice_directory;
        entry.se_object_path = *se_object_path;
        entry.me_object_path = *me_object_path;
        entry.cf_object_path = *cf_object_path;
        entry.projection_x_path = *projection_x_path;
        entry.projection_y_path = *projection_y_path;
        entry.projection_z_path = *projection_z_path;
        entry.centrality_index = *centrality_index;
        entry.mt_index = *mt_index;
        entry.qn_index = qn_index ? **qn_index : -1;
        entry.phi_index = *phi_index;
        entry.cent_low = *cent_low;
        entry.cent_high = *cent_high;
        entry.mt_low = *mt_low;
        entry.mt_high = *mt_high;
        entry.qn_low = qn_low ? **qn_low : std::numeric_limits<double>::quiet_NaN();
        entry.qn_high = qn_high ? **qn_high : std::numeric_limits<double>::quiet_NaN();
        entry.qn_label = qn_label ? **qn_label : kQnIntegratedLabel;
        entry.raw_phi_low = *raw_phi_low;
        entry.raw_phi_high = *raw_phi_high;
        entry.raw_phi_center = *raw_phi_center;
        entry.display_phi_low = *display_phi_low;
        entry.display_phi_high = *display_phi_high;
        entry.display_phi_center = *display_phi_center;
        entry.build_uses_symmetric_phi_range =
            build_uses_symmetric_phi_range ? (**build_uses_symmetric_phi_range != 0) : false;
        entry.split_mixed_event_by_phi =
            split_mixed_event_by_phi ? (**split_mixed_event_by_phi != 0) : false;
        entry.split_mixed_event_by_qn = split_mixed_event_by_qn ? (**split_mixed_event_by_qn != 0) : false;
        entry.is_qn_integrated = is_qn_integrated ? (**is_qn_integrated != 0) : true;
        entry.is_phi_integrated = (*is_phi_integrated != 0);
        entries.push_back(entry);
      }
      if (!build_uses_symmetric_phi_range) {
        InferLegacyBuildPhiMappingState(entries, logger);
      }
      NormalizeCatalogBuildPhiMappingState(entries, logger);
      return entries;
    }

    TH3D *LoadStoredHistogram3D(TFile &input_file, const std::string &object_path, const std::string &clone_name) {
      auto *histogram = dynamic_cast<TH3D *>(input_file.Get(object_path.c_str()));
      if (histogram == nullptr) {
        return nullptr;
      }
      auto *clone = static_cast<TH3D *>(histogram->Clone(clone_name.c_str()));
      clone->SetDirectory(nullptr);
      return clone;
    }

    bool HasCatsFiniteSourceSupport() {
#ifdef EXP_FEMTO_3D_HAS_CATS
      return true;
#else
      return false;
#endif
    }

#ifdef EXP_FEMTO_3D_HAS_CATS
    double CatsSphericalGaussianSource(double *parameters) {
      const double radius = parameters[1];
      const double source_size = parameters[3];
      return 4.0 * TMath::Pi() * radius * radius
             * std::pow(4.0 * TMath::Pi() * source_size * source_size, -1.5)
             * std::exp(-(radius * radius) / (4.0 * source_size * source_size));
    }
#endif

    CoulombKernelTable BuildFiniteSourceCoulombKernel(const SliceCatalogEntry &entry,
                                                      const FiniteSourceMode finite_source_mode,
                                                      const double seed_radius_fm,
                                                      const double final_radius_fm) {
      if (!std::isfinite(final_radius_fm) || final_radius_fm <= 0.0) {
        throw std::runtime_error("Finite-source Coulomb kernel requires a finite positive R_eff.");
      }
      CoulombKernelTable table;
      table.catalog_entry.group_id = entry.group_id;
      table.catalog_entry.centrality_index = entry.centrality_index;
      table.catalog_entry.mt_index = entry.mt_index;
      table.catalog_entry.qn_index = entry.qn_index;
      table.catalog_entry.cent_low = entry.cent_low;
      table.catalog_entry.cent_high = entry.cent_high;
      table.catalog_entry.mt_low = entry.mt_low;
      table.catalog_entry.mt_high = entry.mt_high;
      table.catalog_entry.qn_low = entry.qn_low;
      table.catalog_entry.qn_high = entry.qn_high;
      table.catalog_entry.qn_label = entry.qn_label;
      table.catalog_entry.is_qn_integrated = entry.is_qn_integrated;
      table.catalog_entry.finite_source_mode = ToString(finite_source_mode);
      table.catalog_entry.seed_radius_fm = seed_radius_fm;
      table.catalog_entry.final_radius_fm = final_radius_fm;
      table.catalog_entry.kstar_min_mev = kFiniteSourceKStarMinMeV;
      table.catalog_entry.kstar_max_mev = kFiniteSourceKStarMaxMeV;
      table.catalog_entry.kstar_bin_count = kFiniteSourceKStarBins;

#ifdef EXP_FEMTO_3D_HAS_CATS
      table.catalog_entry.cats_enabled = true;
      CATS kitty;
      kitty.SetNotifications(CATS::nError);
      kitty.SetMomBins(kFiniteSourceKStarBins, kFiniteSourceKStarMinMeV, kFiniteSourceKStarMaxMeV);
      CATSparameters source_parameters(CATSparameters::tSource, 1, true);
      source_parameters.SetParameter(0, final_radius_fm);
      kitty.SetAnaSource(CatsSphericalGaussianSource, source_parameters);
      kitty.SetUseAnalyticSource(true);
      kitty.SetNumChannels(1);
      kitty.SetNumPW(0, 0);
      kitty.SetSpin(0, 0);
      kitty.SetChannelWeight(0, 1.0);
      kitty.SetQ1Q2(1);
      kitty.SetPdgId(211, 211);
      kitty.SetQuantumStatistics(0);
      kitty.SetRedMass(0.5 * kChargedPionMassMeV);
      kitty.KillTheCat();

      table.kstar_mev.reserve(kFiniteSourceKStarBins);
      table.coulomb_factor.reserve(kFiniteSourceKStarBins);
      for (int bin = 0; bin < kFiniteSourceKStarBins; ++bin) {
        const double momentum = kitty.GetMomentum(static_cast<unsigned>(bin));
        double value = kitty.GetCorrFun(static_cast<unsigned>(bin));
        if (!std::isfinite(value) || value < 0.0) {
          throw std::runtime_error("CATS returned a non-finite or negative finite-source Coulomb value inside the kernel table.");
        }
        table.kstar_mev.push_back(momentum);
        table.coulomb_factor.push_back(value);
      }
      return table;
#else
      (void)entry;
      (void)finite_source_mode;
      (void)seed_radius_fm;
      throw std::runtime_error("fit.coulomb_mode=\"finite_source\" requires CATS support, but this build was configured without EXP_FEMTO_3D_HAS_CATS.");
#endif
    }

    std::string BuildSparseObjectPath(const ApplicationConfig &config, const std::string &subtask_name) {
      return config.input.task_name + "/" + subtask_name + "/" + config.input.sparse_object_name;
    }

    bool MatchSelectedBin(const SliceCatalogEntry &entry,
                          const std::vector<RangeBin> &centrality_bins,
                          const std::vector<RangeBin> &mt_bins) {
      const bool centrality_matched =
          std::any_of(centrality_bins.begin(), centrality_bins.end(), [&](const RangeBin &bin) {
            return NearlyEqual(entry.cent_low, bin.min) && NearlyEqual(entry.cent_high, bin.max);
          });
      if (!centrality_matched) {
        return false;
      }

      return std::any_of(mt_bins.begin(), mt_bins.end(), [&](const RangeBin &bin) {
        return NearlyEqual(entry.mt_low, bin.min) && NearlyEqual(entry.mt_high, bin.max);
      });
    }

    double ComputeLikeSignPiPiGamowFactor(const double q_out, const double q_side, const double q_long) {
      const double q_magnitude = std::sqrt(q_out * q_out + q_side * q_side + q_long * q_long);
      const double k_star_fm = 0.5 * q_magnitude / kHbarC;
      if (k_star_fm <= 1.0e-12) {
        return 0.0;
      }

      const double eta = 1.0 / (k_star_fm * kPiPiLikeSignBohrRadiusFm);
      const double two_pi_eta = 2.0 * TMath::Pi() * eta;
      if (two_pi_eta > 700.0) {
        return 0.0;
      }

      const double denominator = std::exp(two_pi_eta) - 1.0;
      if (denominator <= 0.0) {
        return 0.0;
      }
      return std::max(0.0, std::min(two_pi_eta / denominator, 1.0));
    }

    CoulombMode CoulombModeFromCode(const double code) {
      if (code < 0.5) {
        return CoulombMode::kNone;
      }
      if (code < 1.5) {
        return CoulombMode::kGamow;
      }
      return CoulombMode::kFiniteSource;
    }

    double EvaluateCoulombFactor(const double q_out,
                                 const double q_side,
                                 const double q_long,
                                 const LevyFitOptions &fit_options) {
      switch (fit_options.coulomb_mode) {
        case CoulombMode::kNone:
          return 1.0;
        case CoulombMode::kGamow:
          return ComputeLikeSignPiPiGamowFactor(q_out, q_side, q_long);
        case CoulombMode::kFiniteSource:
          if (g_active_coulomb_kernel == nullptr) {
            return kInvalidFullModelCFValue;
          }
          return g_active_coulomb_kernel->Evaluate(q_out, q_side, q_long);
      }
      return 1.0;
    }

    double ComputeBowlerSinyukovLikeSignPiPiValue(const double norm,
                                                  const double lambda,
                                                  const double levy_exponent,
                                                  const LevyFitOptions &fit_options,
                                                  const double q_out,
                                                  const double q_side,
                                                  const double q_long) {
      const double lambda_eff = fit_options.use_core_halo_lambda ? lambda : 1.0;
      const double coulomb_factor = EvaluateCoulombFactor(q_out, q_side, q_long, fit_options);
      const double quantum_stat_term = std::exp(-levy_exponent);
      return norm * ((1.0 - lambda_eff) + lambda_eff * coulomb_factor * (1.0 + quantum_stat_term));
    }

    double ComputeQ2BaselineFactor(const double q_out,
                                   const double q_side,
                                   const double q_long,
                                   const double baseline_q2,
                                   const LevyFitOptions &fit_options) {
      if (!fit_options.use_q2_baseline) {
        return 1.0;
      }
      const double q2 = q_out * q_out + q_side * q_side + q_long * q_long;
      return 1.0 + baseline_q2 * q2;
    }

    bool IsFullR2MatrixPositiveSemiDefinite(const double rout2,
                                            const double rside2,
                                            const double rlong2,
                                            const double routside2,
                                            const double routlong2,
                                            const double rsidelong2,
                                            const double tolerance = kFullR2MatrixTolerance) {
      if (rout2 < -tolerance || rside2 < -tolerance || rlong2 < -tolerance) {
        return false;
      }

      const double det_out_side = rout2 * rside2 - routside2 * routside2;
      const double det_out_long = rout2 * rlong2 - routlong2 * routlong2;
      const double det_side_long = rside2 * rlong2 - rsidelong2 * rsidelong2;
      if (det_out_side < -tolerance || det_out_long < -tolerance || det_side_long < -tolerance) {
        return false;
      }

      const double determinant = rout2 * (rside2 * rlong2 - rsidelong2 * rsidelong2)
                                 - routside2 * (routside2 * rlong2 - routlong2 * rsidelong2)
                                 + routlong2 * (routside2 * rsidelong2 - routlong2 * rside2);
      return determinant >= -tolerance;
    }

    bool HasValidFullR2MatrixFromParameterArray(const double *parameters) {
      if (parameters == nullptr) {
        return false;
      }
      return IsFullR2MatrixPositiveSemiDefinite(
          parameters[2], parameters[3], parameters[4], parameters[5], parameters[6], parameters[7]);
    }

    double EvaluateDiagonalLevyCF(const double q_out,
                                  const double q_side,
                                  const double q_long,
                                  const double norm,
                                  const double lambda,
                                  const double rout2,
                                  const double rside2,
                                  const double rlong2,
                                  const double alpha,
                                  const double baseline_q2,
                                  const LevyFitOptions &fit_options) {
      const double q_out2 = q_out * q_out;
      const double q_side2 = q_side * q_side;
      const double q_long2 = q_long * q_long;
      const double argument = (rout2 * q_out2 + rside2 * q_side2 + rlong2 * q_long2) / (kHbarC * kHbarC);
      const double levy_exponent = std::pow(std::max(argument, 0.0), alpha / 2.0);
      const double femto_value =
          ComputeBowlerSinyukovLikeSignPiPiValue(norm, lambda, levy_exponent, fit_options, q_out, q_side, q_long);
      return femto_value * ComputeQ2BaselineFactor(q_out, q_side, q_long, baseline_q2, fit_options);
    }

    double EvaluateFullLevyCF(const double q_out,
                              const double q_side,
                              const double q_long,
                              const double norm,
                              const double lambda,
                              const double rout2,
                              const double rside2,
                              const double rlong2,
                              const double routside2,
                              const double routlong2,
                              const double rsidelong2,
                              const double alpha,
                              const double baseline_q2,
                              const LevyFitOptions &fit_options) {
      if (!IsFullR2MatrixPositiveSemiDefinite(rout2, rside2, rlong2, routside2, routlong2, rsidelong2)) {
        return kInvalidFullModelCFValue;
      }

      const double argument =
          (rout2 * q_out * q_out + rside2 * q_side * q_side + rlong2 * q_long * q_long
           + 2.0 * routside2 * q_out * q_side + 2.0 * routlong2 * q_out * q_long + 2.0 * rsidelong2 * q_side * q_long)
          / (kHbarC * kHbarC);
      if (argument < -kLevyArgumentTolerance) {
        return kInvalidFullModelCFValue;
      }
      const double protected_argument = argument < 0.0 ? 0.0 : argument;
      const double levy_exponent = std::pow(protected_argument, alpha / 2.0);
      const double femto_value =
          ComputeBowlerSinyukovLikeSignPiPiValue(norm, lambda, levy_exponent, fit_options, q_out, q_side, q_long);
      return femto_value * ComputeQ2BaselineFactor(q_out, q_side, q_long, baseline_q2, fit_options);
    }

    double Levy3DModel(double *x, double *parameters) {
      LevyFitOptions fit_options;
      fit_options.use_q2_baseline = parameters[7] > 0.5;
      fit_options.coulomb_mode = CoulombModeFromCode(parameters[8]);
      fit_options.use_core_halo_lambda = parameters[9] > 0.5;
      return EvaluateDiagonalLevyCF(x[0],
                                    x[1],
                                    x[2],
                                    parameters[0],
                                    parameters[1],
                                    parameters[2],
                                    parameters[3],
                                    parameters[4],
                                    parameters[5],
                                    parameters[6],
                                    fit_options);
    }

    double Levy3DFullModel(double *x, double *parameters) {
      LevyFitOptions fit_options;
      fit_options.use_q2_baseline = parameters[10] > 0.5;
      fit_options.coulomb_mode = CoulombModeFromCode(parameters[11]);
      fit_options.use_core_halo_lambda = parameters[12] > 0.5;
      return EvaluateFullLevyCF(x[0],
                                x[1],
                                x[2],
                                parameters[0],
                                parameters[1],
                                parameters[2],
                                parameters[3],
                                parameters[4],
                                parameters[5],
                                parameters[6],
                                parameters[7],
                                parameters[8],
                                parameters[9],
                                fit_options);
    }

    struct EffectiveLevyFitParameter {
      double initial = 0.0;
      double lower = 0.0;
      double upper = 0.0;
      bool has_limits = false;
      std::optional<double> fixed_value;
    };

    // Merge sparse TOML overrides with the legacy defaults at the last moment so
    // omitted fields keep the historical fit seed, bounds, and fixed-parameter policy.
    EffectiveLevyFitParameter MakeLevyFitParameter(const double default_initial,
                                                   const std::optional<std::pair<double, double>> default_limits,
                                                   const LevyFitParameterOverride &parameter) {
      EffectiveLevyFitParameter effective;
      effective.initial = parameter.fixed_value.value_or(parameter.initial.value_or(default_initial));
      effective.fixed_value = parameter.fixed_value;
      if (parameter.min.has_value() && parameter.max.has_value()) {
        effective.lower = *parameter.min;
        effective.upper = *parameter.max;
        effective.has_limits = true;
      } else if (default_limits.has_value()) {
        effective.lower = default_limits->first;
        effective.upper = default_limits->second;
        effective.has_limits = true;
      }
      return effective;
    }

    void ApplyLevyFitParameter(TF3 *fit_function, const int index, const EffectiveLevyFitParameter &parameter) {
      fit_function->SetParameter(index, parameter.initial);
      if (parameter.has_limits) {
        fit_function->SetParLimits(index, parameter.lower, parameter.upper);
      }
      if (parameter.fixed_value.has_value()) {
        fit_function->FixParameter(index, *parameter.fixed_value);
      }
    }

    bool IsUserFixedLevyFitParameter(const int parameter_index,
                                     const bool use_full_model,
                                     const LevyFitOptions &fit_options) {
      if (parameter_index == 1) {
        return fit_options.parameters.lambda.fixed_value.has_value();
      }
      if (!use_full_model && parameter_index == 5) {
        return fit_options.parameters.alpha.fixed_value.has_value();
      }
      if (use_full_model && parameter_index == 8) {
        return fit_options.parameters.alpha.fixed_value.has_value();
      }
      return false;
    }

    TF3 *BuildLevyFitFunction(const std::string &function_name, const LevyFitOptions &fit_options) {
      const double q2_max = 3.0 * fit_options.fit_q_max * fit_options.fit_q_max;
      const double baseline_min = q2_max > 0.0 ? -0.9 / q2_max : -10.0;
      const double baseline_max = q2_max > 0.0 ? 2.0 / q2_max : 10.0;
      auto *fit_function = new TF3(function_name.c_str(),
                                   Levy3DModel,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   10);
      fit_function->SetParName(0, "Norm");
      fit_function->SetParName(1, "lambda");
      fit_function->SetParName(2, "Rout2");
      fit_function->SetParName(3, "Rside2");
      fit_function->SetParName(4, "Rlong2");
      fit_function->SetParName(5, "alpha");
      fit_function->SetParName(6, "BaselineQ2");
      fit_function->SetParName(7, "UseQ2Baseline");
      fit_function->SetParName(8, "CoulombModeCode");
      fit_function->SetParName(9, "UseCoreHaloLambda");
      ApplyLevyFitParameter(
          fit_function, 0, MakeLevyFitParameter(1.0, std::make_pair(0.5, 1.5), fit_options.parameters.norm));
      ApplyLevyFitParameter(
          fit_function, 1, MakeLevyFitParameter(0.5, std::make_pair(0.0, 1.0), fit_options.parameters.lambda));
      ApplyLevyFitParameter(
          fit_function, 2, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rout2));
      ApplyLevyFitParameter(
          fit_function, 3, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rside2));
      ApplyLevyFitParameter(
          fit_function, 4, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rlong2));
      ApplyLevyFitParameter(
          fit_function, 5, MakeLevyFitParameter(1.5, std::make_pair(0.5, 2.0), fit_options.parameters.alpha));
      ApplyLevyFitParameter(fit_function,
                            6,
                            MakeLevyFitParameter(
                                0.0, std::make_pair(baseline_min, baseline_max), fit_options.parameters.baseline_q2));
      fit_function->SetParameter(7, 0.0);
      fit_function->SetParameter(8, 0.0);
      fit_function->SetParameter(9, 1.0);
      fit_function->FixParameter(7, fit_options.use_q2_baseline ? 1.0 : 0.0);
      fit_function->FixParameter(8, static_cast<double>(CoulombModeCode(fit_options.coulomb_mode)));
      fit_function->FixParameter(9, fit_options.use_core_halo_lambda ? 1.0 : 0.0);
      if (!fit_options.use_q2_baseline) {
        fit_function->FixParameter(6, 0.0);
      }
      if (!fit_options.use_core_halo_lambda) {
        fit_function->FixParameter(1, 1.0);
      }
      fit_function->SetNpx(60);
      fit_function->SetNpy(60);
      fit_function->SetNpz(60);
      return fit_function;
    }

    TF3 *BuildFullLevyFitFunction(const std::string &function_name, const LevyFitOptions &fit_options) {
      const double q2_max = 3.0 * fit_options.fit_q_max * fit_options.fit_q_max;
      const double baseline_min = q2_max > 0.0 ? -0.9 / q2_max : -10.0;
      const double baseline_max = q2_max > 0.0 ? 2.0 / q2_max : 10.0;
      auto *fit_function = new TF3(function_name.c_str(),
                                   Levy3DFullModel,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   -fit_options.fit_q_max,
                                   fit_options.fit_q_max,
                                   13);
      fit_function->SetParName(0, "Norm");
      fit_function->SetParName(1, "lambda");
      fit_function->SetParName(2, "Rout2");
      fit_function->SetParName(3, "Rside2");
      fit_function->SetParName(4, "Rlong2");
      fit_function->SetParName(5, "Routside2");
      fit_function->SetParName(6, "Routlong2");
      fit_function->SetParName(7, "Rsidelong2");
      fit_function->SetParName(8, "alpha");
      fit_function->SetParName(9, "BaselineQ2");
      fit_function->SetParName(10, "UseQ2Baseline");
      fit_function->SetParName(11, "CoulombModeCode");
      fit_function->SetParName(12, "UseCoreHaloLambda");
      ApplyLevyFitParameter(
          fit_function, 0, MakeLevyFitParameter(1.0, std::make_pair(0.5, 1.5), fit_options.parameters.norm));
      ApplyLevyFitParameter(
          fit_function, 1, MakeLevyFitParameter(0.5, std::make_pair(0.0, 1.0), fit_options.parameters.lambda));
      ApplyLevyFitParameter(
          fit_function, 2, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rout2));
      ApplyLevyFitParameter(
          fit_function, 3, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rside2));
      ApplyLevyFitParameter(
          fit_function, 4, MakeLevyFitParameter(25.0, std::make_pair(0.01, 400.0), fit_options.parameters.rlong2));
      ApplyLevyFitParameter(fit_function, 5, MakeLevyFitParameter(0.0, std::nullopt, fit_options.parameters.routside2));
      ApplyLevyFitParameter(fit_function, 6, MakeLevyFitParameter(0.0, std::nullopt, fit_options.parameters.routlong2));
      ApplyLevyFitParameter(fit_function, 7, MakeLevyFitParameter(0.0, std::nullopt, fit_options.parameters.rsidelong2));
      ApplyLevyFitParameter(
          fit_function, 8, MakeLevyFitParameter(1.5, std::make_pair(0.5, 2.0), fit_options.parameters.alpha));
      ApplyLevyFitParameter(fit_function,
                            9,
                            MakeLevyFitParameter(
                                0.0, std::make_pair(baseline_min, baseline_max), fit_options.parameters.baseline_q2));
      fit_function->SetParameter(10, 0.0);
      fit_function->SetParameter(11, 0.0);
      fit_function->SetParameter(12, 1.0);
      fit_function->FixParameter(10, fit_options.use_q2_baseline ? 1.0 : 0.0);
      fit_function->FixParameter(11, static_cast<double>(CoulombModeCode(fit_options.coulomb_mode)));
      fit_function->FixParameter(12, fit_options.use_core_halo_lambda ? 1.0 : 0.0);
      if (!fit_options.use_q2_baseline) {
        fit_function->FixParameter(9, 0.0);
      }
      if (!fit_options.use_core_halo_lambda) {
        fit_function->FixParameter(1, 1.0);
      }
      fit_function->SetNpx(60);
      fit_function->SetNpy(60);
      fit_function->SetNpz(60);
      return fit_function;
    }

    double EvaluateLevyModelFromParameterArray(const double q_out,
                                               const double q_side,
                                               const double q_long,
                                               const double *parameters,
                                               const bool use_full_model) {
      double x[3] = {q_out, q_side, q_long};
      return use_full_model ? Levy3DFullModel(x, const_cast<double *>(parameters))
                            : Levy3DModel(x, const_cast<double *>(parameters));
    }

    double ComputeRawToNormalizedCFScale(TH3D *h_se_raw, TH3D *h_me_raw) {
      if (h_se_raw == nullptr || h_me_raw == nullptr) {
        return 0.0;
      }
      const double int_se = IntegralVisibleRange(h_se_raw, true);
      const double int_me = IntegralVisibleRange(h_me_raw, true);
      if (int_se <= 0.0 || int_me <= 0.0) {
        return 0.0;
      }
      return int_se / int_me;
    }

    double ComputePMLNeg2LogLContribution(const double same_counts,
                                          const double mixed_counts,
                                          const double model_ratio) {
      if (same_counts < 0.0 || mixed_counts < 0.0 || model_ratio <= 0.0 || !std::isfinite(model_ratio)) {
        return kFitPenaltyValue;
      }
      if (same_counts == 0.0 && mixed_counts == 0.0) {
        return 0.0;
      }
      if (same_counts == 0.0) {
        return 2.0 * mixed_counts * std::log1p(model_ratio);
      }
      if (mixed_counts == 0.0) {
        const double log_term =
            model_ratio >= 1.0 ? std::log1p(1.0 / model_ratio) : std::log1p(model_ratio) - std::log(model_ratio);
        return 2.0 * same_counts * log_term;
      }

      const double total_counts = same_counts + mixed_counts;
      const double arg1 = model_ratio * total_counts / (same_counts * (model_ratio + 1.0));
      const double arg2 = total_counts / (mixed_counts * (model_ratio + 1.0));
      if (arg1 <= 0.0 || arg2 <= 0.0 || !std::isfinite(arg1) || !std::isfinite(arg2)) {
        return kFitPenaltyValue;
      }
      return -2.0 * (same_counts * std::log(arg1) + mixed_counts * std::log(arg2));
    }

    void Levy3DPMLFCN(Int_t &npar, Double_t *grad, Double_t &f, Double_t *parameters, Int_t flag) {
      (void)npar;
      (void)grad;
      (void)flag;
      if (g_levy_3d_pml_context.h_se_raw == nullptr || g_levy_3d_pml_context.h_me_raw == nullptr) {
        f = kFitPenaltyValue;
        return;
      }
      if (g_levy_3d_pml_context.use_full_model && !HasValidFullR2MatrixFromParameterArray(parameters)) {
        f = kFitPenaltyValue;
        return;
      }

      double neg2_log_l = 0.0;
      for (int ix = 1; ix <= g_levy_3d_pml_context.h_se_raw->GetNbinsX(); ++ix) {
        const double q_out = g_levy_3d_pml_context.h_se_raw->GetXaxis()->GetBinCenter(ix);
        if (std::abs(q_out) > g_levy_3d_pml_context.fit_options.fit_q_max) {
          continue;
        }
        for (int iy = 1; iy <= g_levy_3d_pml_context.h_se_raw->GetNbinsY(); ++iy) {
          const double q_side = g_levy_3d_pml_context.h_se_raw->GetYaxis()->GetBinCenter(iy);
          if (std::abs(q_side) > g_levy_3d_pml_context.fit_options.fit_q_max) {
            continue;
          }
          for (int iz = 1; iz <= g_levy_3d_pml_context.h_se_raw->GetNbinsZ(); ++iz) {
            const double q_long = g_levy_3d_pml_context.h_se_raw->GetZaxis()->GetBinCenter(iz);
            if (std::abs(q_long) > g_levy_3d_pml_context.fit_options.fit_q_max) {
              continue;
            }

            const double same_counts = g_levy_3d_pml_context.h_se_raw->GetBinContent(ix, iy, iz);
            const double mixed_counts = g_levy_3d_pml_context.h_me_raw->GetBinContent(ix, iy, iz);
            if (same_counts == 0.0 && mixed_counts == 0.0) {
              continue;
            }

            double model_ratio = EvaluateLevyModelFromParameterArray(
                q_out, q_side, q_long, parameters, g_levy_3d_pml_context.use_full_model);
            model_ratio *= g_levy_3d_pml_context.raw_same_to_mixed_integral_ratio;
            const double contribution = ComputePMLNeg2LogLContribution(same_counts, mixed_counts, model_ratio);
            if (!std::isfinite(contribution) || contribution >= kFitPenaltyValue) {
              f = kFitPenaltyValue;
              return;
            }
            neg2_log_l += contribution;
          }
        }
      }

      f = std::isfinite(neg2_log_l) ? neg2_log_l : kFitPenaltyValue;
    }

    int CountPMLUsableBins(TH3D *h_se_raw, TH3D *h_me_raw, const double q_max) {
      if (h_se_raw == nullptr || h_me_raw == nullptr) {
        return 0;
      }
      int n_points = 0;
      for (int ix = 1; ix <= h_se_raw->GetNbinsX(); ++ix) {
        const double q_out = h_se_raw->GetXaxis()->GetBinCenter(ix);
        if (std::abs(q_out) > q_max) {
          continue;
        }
        for (int iy = 1; iy <= h_se_raw->GetNbinsY(); ++iy) {
          const double q_side = h_se_raw->GetYaxis()->GetBinCenter(iy);
          if (std::abs(q_side) > q_max) {
            continue;
          }
          for (int iz = 1; iz <= h_se_raw->GetNbinsZ(); ++iz) {
            const double q_long = h_se_raw->GetZaxis()->GetBinCenter(iz);
            if (std::abs(q_long) > q_max) {
              continue;
            }
            if (h_se_raw->GetBinContent(ix, iy, iz) + h_me_raw->GetBinContent(ix, iy, iz) > 0.0) {
              ++n_points;
            }
          }
        }
      }
      return n_points;
    }

    double EstimatePMLStepSize(const int parameter_index, const bool use_full_model) {
      if ((!use_full_model && (parameter_index >= 2 && parameter_index <= 4))
          || (use_full_model && (parameter_index >= 2 && parameter_index <= 7))) {
        return 0.1;
      }
      return 0.01;
    }

    bool IsPMLParameterFixed(const int parameter_index, const bool use_full_model, const LevyFitOptions &fit_options) {
      if (!use_full_model) {
        if (parameter_index >= 7) {
          return true;
        }
        if (IsUserFixedLevyFitParameter(parameter_index, use_full_model, fit_options)) {
          return true;
        }
        if (parameter_index == 1 && !fit_options.use_core_halo_lambda) {
          return true;
        }
        if (parameter_index == 6 && !fit_options.use_q2_baseline) {
          return true;
        }
        return false;
      }

      if (parameter_index >= 10) {
        return true;
      }
      if (IsUserFixedLevyFitParameter(parameter_index, use_full_model, fit_options)) {
        return true;
      }
      if (parameter_index == 1 && !fit_options.use_core_halo_lambda) {
        return true;
      }
      if (parameter_index == 9 && !fit_options.use_q2_baseline) {
        return true;
      }
      return false;
    }

    void ConfigurePMLMinuit(TMinuit &minuit,
                            TF3 *fit_function,
                            const bool use_full_model,
                            const LevyFitOptions &fit_options) {
      const int n_parameters = fit_function->GetNpar();
      Int_t error_code = 0;
      for (int index = 0; index < n_parameters; ++index) {
        double lower = 0.0;
        double upper = 0.0;
        fit_function->GetParLimits(index, lower, upper);
        const double value = fit_function->GetParameter(index);
        const double step = EstimatePMLStepSize(index, use_full_model);
        const bool fixed = IsPMLParameterFixed(index, use_full_model, fit_options);
        if (fixed) {
          minuit.mnparm(index, fit_function->GetParName(index), value, step, 0.0, 0.0, error_code);
          minuit.FixParameter(index);
        } else {
          minuit.mnparm(index, fit_function->GetParName(index), value, step, lower, upper, error_code);
        }
      }
    }

    bool RunPMLFit(TF3 *fit_function,
                   TH3D *h_se_raw,
                   TH3D *h_me_raw,
                   const bool use_full_model,
                   const LevyFitOptions &fit_options,
                   double &fit_statistic,
                   int &ndf,
                   int &fit_status,
                   double &edm,
                   int &minuit_istat) {
      if (fit_function == nullptr || h_se_raw == nullptr || h_me_raw == nullptr) {
        return false;
      }
      const double raw_to_normalized_scale = ComputeRawToNormalizedCFScale(h_se_raw, h_me_raw);
      if (raw_to_normalized_scale <= 0.0 || !std::isfinite(raw_to_normalized_scale)) {
        return false;
      }

      g_levy_3d_pml_context.h_se_raw = h_se_raw;
      g_levy_3d_pml_context.h_me_raw = h_me_raw;
      g_levy_3d_pml_context.use_full_model = use_full_model;
      g_levy_3d_pml_context.fit_options = fit_options;
      g_levy_3d_pml_context.coulomb_kernel = g_active_coulomb_kernel;
      g_levy_3d_pml_context.raw_same_to_mixed_integral_ratio = raw_to_normalized_scale;

      ActiveCoulombKernelGuard kernel_guard(g_levy_3d_pml_context.coulomb_kernel);
      TMinuit minuit(fit_function->GetNpar());
      minuit.SetFCN(Levy3DPMLFCN);
      minuit.SetPrintLevel(-1);
      minuit.SetErrorDef(1.0);
      ConfigurePMLMinuit(minuit, fit_function, use_full_model, fit_options);

      Int_t error_code = 0;
      Double_t arglist[2];
      arglist[0] = 100000;
      arglist[1] = 0.1;
      minuit.mnexcm("MIGRAD", arglist, 2, error_code);
      const Int_t migrad_error = error_code;
      arglist[0] = 0;
      minuit.mnexcm("HESSE", arglist, 1, error_code);

      Double_t fmin = 0.0;
      Double_t fedm = 0.0;
      Double_t errdef = 0.0;
      Int_t npari = 0;
      Int_t nparx = 0;
      Int_t istat = 0;
      minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
      (void)errdef;
      (void)nparx;

      for (int index = 0; index < fit_function->GetNpar(); ++index) {
        double value = 0.0;
        double error = 0.0;
        minuit.GetParameter(index, value, error);
        fit_function->SetParameter(index, value);
        fit_function->SetParError(index, error);
      }

      fit_statistic = fmin;
      edm = fedm;
      ndf = std::max(0, CountPMLUsableBins(h_se_raw, h_me_raw, fit_options.fit_q_max) - npari);
      minuit_istat = istat;
      fit_status = migrad_error != 0 ? -migrad_error : istat;
      return true;
    }

    void StyleProjectionHistogram(TH1D *histogram,
                                  const int color,
                                  const int marker_style,
                                  const std::string &x_title) {
      histogram->SetLineColor(color);
      histogram->SetMarkerColor(color);
      histogram->SetMarkerStyle(marker_style);
      histogram->SetMarkerSize(0.9);
      histogram->SetLineWidth(2);
      histogram->GetXaxis()->SetTitle(x_title.c_str());
      histogram->GetYaxis()->SetTitle("C(q)");
    }

    void StyleProjectionCurve(TGraph *graph, const int color) {
      graph->SetLineColor(color);
      graph->SetLineWidth(3);
      graph->SetMarkerSize(0.0);
    }

    double GetGraphMaximum(const TGraph *graph) {
      if (graph == nullptr || graph->GetN() <= 0) {
        return 0.0;
      }
      double x = 0.0;
      double y = 0.0;
      graph->GetPoint(0, x, y);
      double maximum = y;
      for (int index = 1; index < graph->GetN(); ++index) {
        graph->GetPoint(index, x, y);
        maximum = std::max(maximum, y);
      }
      return maximum;
    }

    double EvaluateProjectionXWindowAverage(TF3 *fit_function,
                                            const TH3D *reference_histogram,
                                            const double q_out,
                                            const double q_max) {
      auto [y_min, y_max] = GetAxisRangeForWindow(reference_histogram->GetYaxis(), q_max);
      auto [z_min, z_max] = GetAxisRangeForWindow(reference_histogram->GetZaxis(), q_max);
      double value_sum = 0.0;
      int n_points = 0;
      for (int iy = y_min; iy <= y_max; ++iy) {
        const double q_side = reference_histogram->GetYaxis()->GetBinCenter(iy);
        for (int iz = z_min; iz <= z_max; ++iz) {
          const double q_long = reference_histogram->GetZaxis()->GetBinCenter(iz);
          value_sum += fit_function->Eval(q_out, q_side, q_long);
          ++n_points;
        }
      }
      return n_points > 0 ? value_sum / static_cast<double>(n_points) : 0.0;
    }

    double EvaluateProjectionYWindowAverage(TF3 *fit_function,
                                            const TH3D *reference_histogram,
                                            const double q_side,
                                            const double q_max) {
      auto [x_min, x_max] = GetAxisRangeForWindow(reference_histogram->GetXaxis(), q_max);
      auto [z_min, z_max] = GetAxisRangeForWindow(reference_histogram->GetZaxis(), q_max);
      double value_sum = 0.0;
      int n_points = 0;
      for (int ix = x_min; ix <= x_max; ++ix) {
        const double q_out = reference_histogram->GetXaxis()->GetBinCenter(ix);
        for (int iz = z_min; iz <= z_max; ++iz) {
          const double q_long = reference_histogram->GetZaxis()->GetBinCenter(iz);
          value_sum += fit_function->Eval(q_out, q_side, q_long);
          ++n_points;
        }
      }
      return n_points > 0 ? value_sum / static_cast<double>(n_points) : 0.0;
    }

    double EvaluateProjectionZWindowAverage(TF3 *fit_function,
                                            const TH3D *reference_histogram,
                                            const double q_long,
                                            const double q_max) {
      auto [x_min, x_max] = GetAxisRangeForWindow(reference_histogram->GetXaxis(), q_max);
      auto [y_min, y_max] = GetAxisRangeForWindow(reference_histogram->GetYaxis(), q_max);
      double value_sum = 0.0;
      int n_points = 0;
      for (int ix = x_min; ix <= x_max; ++ix) {
        const double q_out = reference_histogram->GetXaxis()->GetBinCenter(ix);
        for (int iy = y_min; iy <= y_max; ++iy) {
          const double q_side = reference_histogram->GetYaxis()->GetBinCenter(iy);
          value_sum += fit_function->Eval(q_out, q_side, q_long);
          ++n_points;
        }
      }
      return n_points > 0 ? value_sum / static_cast<double>(n_points) : 0.0;
    }

    TGraph *BuildProjectionCurveXWithinWindow(TF3 *fit_function,
                                              const TH3D *reference_histogram,
                                              const std::string &name,
                                              const double q_max,
                                              const int n_samples = 400) {
      auto *graph = new TGraph(n_samples);
      graph->SetName(name.c_str());
      const double x_min = reference_histogram->GetXaxis()->GetBinLowEdge(1);
      const double x_max = reference_histogram->GetXaxis()->GetBinUpEdge(reference_histogram->GetNbinsX());
      for (int index = 0; index < n_samples; ++index) {
        const double x = x_min + (x_max - x_min) * static_cast<double>(index) / static_cast<double>(n_samples - 1);
        graph->SetPoint(index, x, EvaluateProjectionXWindowAverage(fit_function, reference_histogram, x, q_max));
      }
      return graph;
    }

    TGraph *BuildProjectionCurveYWithinWindow(TF3 *fit_function,
                                              const TH3D *reference_histogram,
                                              const std::string &name,
                                              const double q_max,
                                              const int n_samples = 400) {
      auto *graph = new TGraph(n_samples);
      graph->SetName(name.c_str());
      const double x_min = reference_histogram->GetYaxis()->GetBinLowEdge(1);
      const double x_max = reference_histogram->GetYaxis()->GetBinUpEdge(reference_histogram->GetNbinsY());
      for (int index = 0; index < n_samples; ++index) {
        const double x = x_min + (x_max - x_min) * static_cast<double>(index) / static_cast<double>(n_samples - 1);
        graph->SetPoint(index, x, EvaluateProjectionYWindowAverage(fit_function, reference_histogram, x, q_max));
      }
      return graph;
    }

    TGraph *BuildProjectionCurveZWithinWindow(TF3 *fit_function,
                                              const TH3D *reference_histogram,
                                              const std::string &name,
                                              const double q_max,
                                              const int n_samples = 400) {
      auto *graph = new TGraph(n_samples);
      graph->SetName(name.c_str());
      const double x_min = reference_histogram->GetZaxis()->GetBinLowEdge(1);
      const double x_max = reference_histogram->GetZaxis()->GetBinUpEdge(reference_histogram->GetNbinsZ());
      for (int index = 0; index < n_samples; ++index) {
        const double x = x_min + (x_max - x_min) * static_cast<double>(index) / static_cast<double>(n_samples - 1);
        graph->SetPoint(index, x, EvaluateProjectionZWindowAverage(fit_function, reference_histogram, x, q_max));
      }
      return graph;
    }

    std::string FormatParameterLine(const std::string &label,
                                    const double value,
                                    const double error,
                                    const std::string &unit = "") {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(3) << label << " = " << value << " #pm " << error;
      if (!unit.empty()) {
        stream << " " << unit;
      }
      return stream.str();
    }

    std::string BuildFitModeTitle(const LevyFitResult &fit_result) {
      return fit_result.has_off_diagonal ? "Full Levy fit" : "Diagonal Levy fit";
    }

    std::string BuildFitSwitchLine(const LevyFitResult &fit_result) {
      std::string coulomb_text = "Coulomb: " + fit_result.coulomb_mode;
      if (!fit_result.finite_source_mode.empty()) {
        coulomb_text += "/" + fit_result.finite_source_mode;
      }
      return coulomb_text + ", core-halo: " + (fit_result.uses_core_halo_lambda ? "on" : "off")
             + ", q^{2} baseline: " + (fit_result.uses_q2_baseline ? "on" : "off");
    }

    std::string DescribeCovarianceQuality(const int istat) {
      switch (istat) {
        case 0:
          return "not available";
        case 1:
          return "approximate";
        case 2:
          return "forced pos-def";
        case 3:
          return "full, accurate";
        default:
          return "not applicable";
      }
    }

    std::string CovarianceQualityToken(const int istat) {
      switch (istat) {
        case 0:
          return "not_available";
        case 1:
          return "approximate";
        case 2:
          return "forced_pos_def";
        case 3:
          return "full_accurate";
        default:
          return "not_applicable";
      }
    }

    std::string BuildFitStatisticLine(const LevyFitResult &fit_result) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(2) << (fit_result.uses_pml ? "-2 ln L/NDF = " : "#chi^{2}/NDF = ")
             << fit_result.fit_statistic << "/" << fit_result.ndf;
      return stream.str();
    }

    TPaveText *BuildFitParameterBox(
        const LevyFitResult &fit_result, const double x1, const double y1, const double x2, const double y2) {
      auto *box = new TPaveText(x1, y1, x2, y2, "NDC");
      box->SetFillColor(0);
      box->SetFillStyle(1001);
      box->SetBorderSize(1);
      box->SetTextAlign(12);
      box->SetTextFont(42);
      double text_size = fit_result.has_off_diagonal ? 0.024 : 0.028;
      if (fit_result.uses_pml) {
        text_size = fit_result.has_off_diagonal ? 0.022 : 0.024;
      }
      box->SetTextSize(text_size);
      box->AddText(BuildFitModeTitle(fit_result).c_str());
      box->AddText(BuildFitSwitchLine(fit_result).c_str());
      box->AddText(FormatParameterLine("N", fit_result.norm, fit_result.norm_err).c_str());
      if (fit_result.uses_core_halo_lambda) {
        box->AddText(FormatParameterLine("#lambda", fit_result.lambda, fit_result.lambda_err).c_str());
      } else {
        box->AddText("#lambda fixed = 1.000");
      }
      box->AddText(FormatParameterLine("#alpha", fit_result.alpha, fit_result.alpha_err).c_str());
      if (fit_result.uses_q2_baseline) {
        box->AddText(
            FormatParameterLine("b_{q^{2}}", fit_result.baseline_q2, fit_result.baseline_q2_err, "(GeV/c)^{-2}")
                .c_str());
      } else {
        box->AddText("b_{q^{2}} fixed = 0.000");
      }
      box->AddText(FormatParameterLine("R_{out}^{2}", fit_result.rout2, fit_result.rout2_err, "fm^{2}").c_str());
      box->AddText(FormatParameterLine("R_{side}^{2}", fit_result.rside2, fit_result.rside2_err, "fm^{2}").c_str());
      box->AddText(FormatParameterLine("R_{long}^{2}", fit_result.rlong2, fit_result.rlong2_err, "fm^{2}").c_str());
      if (fit_result.has_off_diagonal) {
        box->AddText(
            FormatParameterLine("R_{outside}^{2}", fit_result.routside2, fit_result.routside2_err, "fm^{2}").c_str());
        box->AddText(
            FormatParameterLine("R_{outlong}^{2}", fit_result.routlong2, fit_result.routlong2_err, "fm^{2}").c_str());
        box->AddText(FormatParameterLine("R_{sidelong}^{2}", fit_result.rsidelong2, fit_result.rsidelong2_err, "fm^{2}")
                         .c_str());
      }
      box->AddText((std::string("Fit method: ") + (fit_result.uses_pml ? "PML" : "chi2")).c_str());
      box->AddText(BuildFitStatisticLine(fit_result).c_str());
      if (fit_result.uses_pml) {
        std::ostringstream edm_line;
        edm_line << std::scientific << std::setprecision(3) << "EDM = " << fit_result.edm;
        box->AddText(edm_line.str().c_str());
        box->AddText((std::string("istat = ") + std::to_string(fit_result.minuit_istat)).c_str());
        box->AddText((std::string("Cov quality: ") + DescribeCovarianceQuality(fit_result.minuit_istat)).c_str());
      }
      return box;
    }

    TCanvas *BuildProjectionCanvas(const std::string &canvas_name,
                                   TH1D *data_histogram,
                                   TGraph *fit_graph,
                                   const std::string &x_title,
                                   const LevyFitResult &fit_result) {
      StyleProjectionHistogram(data_histogram, kBlack, 20, x_title);
      StyleProjectionCurve(fit_graph, kRed + 1);
      const double max_value = std::max(data_histogram->GetMaximum(), GetGraphMaximum(fit_graph));
      data_histogram->SetMaximum(max_value * 1.15);

      auto *canvas = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
      canvas->SetMargin(0.13, 0.05, 0.12, 0.07);
      canvas->SetGrid();
      data_histogram->Draw("E1");
      fit_graph->Draw("L SAME");

      auto *legend = new TLegend(0.62, 0.72, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->AddEntry(data_histogram, "Data projection", "lep");
      legend->AddEntry(fit_graph, "Levy fit projection", "l");
      legend->Draw();

      const double y1 =
          fit_result.has_off_diagonal ? (fit_result.uses_pml ? 0.32 : 0.40) : (fit_result.uses_pml ? 0.42 : 0.50);
      auto *parameter_box = BuildFitParameterBox(fit_result, 0.16, y1, 0.58, 0.88);
      parameter_box->Draw();
      canvas->Update();
      return canvas;
    }

    TCanvas *Build3DComparisonCanvas(const std::string &canvas_name,
                                     TH3D *data_histogram,
                                     TF3 *fit_function,
                                     const LevyFitResult &fit_result) {
      auto *canvas = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 1400, 600);
      canvas->Divide(2, 1);
      canvas->cd(1);
      gPad->SetTheta(24);
      gPad->SetPhi(32);
      data_histogram->Draw("BOX2Z");

      canvas->cd(2);
      gPad->SetTheta(24);
      gPad->SetPhi(32);
      fit_function->GetXaxis()->SetTitle("q_{out} (GeV/c)");
      fit_function->GetYaxis()->SetTitle("q_{side} (GeV/c)");
      fit_function->GetZaxis()->SetTitle("q_{long} (GeV/c)");
      fit_function->SetLineColor(kRed + 1);
      fit_function->SetLineWidth(2);
      fit_function->Draw("ISO");

      const double y1 =
          fit_result.has_off_diagonal ? (fit_result.uses_pml ? 0.18 : 0.28) : (fit_result.uses_pml ? 0.30 : 0.40);
      auto *parameter_box = BuildFitParameterBox(fit_result, 0.12, y1, 0.58, 0.88);
      parameter_box->Draw();
      canvas->Update();
      return canvas;
    }

    struct SingleSliceFitOutput {
      std::unique_ptr<TF3> fit_function;
      LevyFitResult result;
    };

    double ComputeEffectiveRadiusFromResult(const LevyFitResult &result) {
      if (!std::isfinite(result.rout2) || !std::isfinite(result.rside2) || !std::isfinite(result.rlong2)
          || result.rout2 <= 0.0 || result.rside2 <= 0.0 || result.rlong2 <= 0.0) {
        throw std::runtime_error("Finite-source Coulomb seed fit produced non-finite or non-positive diagonal R^2 values for "
                                 + result.slice_id + ".");
      }
      return std::sqrt((result.rout2 + result.rside2 + result.rlong2) / 3.0);
    }

    void FillFitResultFromFunction(const SliceCatalogEntry &entry,
                                   TF3 *fit_function,
                                   const FitModel model,
                                   const LevyFitOptions &fit_options,
                                   const bool fit_uses_symmetric_phi_range,
                                   const double fit_statistic,
                                   const int fit_ndf,
                                   const int fit_status,
                                   const double fit_edm,
                                   const int fit_minuit_istat,
                                   const double finite_source_radius_fm,
                                   LevyFitResult &fit_result) {
      fit_result.fit_model = ToString(model);
      fit_result.slice_id = entry.slice_id;
      fit_result.group_id = entry.group_id;
      fit_result.slice_directory = BuildFitDirectory(entry.slice_id);
      fit_result.centrality_index = entry.centrality_index;
      fit_result.mt_index = entry.mt_index;
      fit_result.qn_index = entry.qn_index;
      fit_result.phi_index = entry.phi_index;
      fit_result.cent_low = entry.cent_low;
      fit_result.cent_high = entry.cent_high;
      fit_result.mt_low = entry.mt_low;
      fit_result.mt_high = entry.mt_high;
      fit_result.qn_low = entry.qn_low;
      fit_result.qn_high = entry.qn_high;
      fit_result.qn_label = entry.qn_label;
      const PhiSliceCoordinates fit_phi_coordinates = ResolveSlicePhiCoordinates(entry, fit_uses_symmetric_phi_range);
      fit_result.is_qn_integrated = entry.is_qn_integrated;
      fit_result.is_phi_integrated = entry.is_phi_integrated;
      fit_result.phi = fit_phi_coordinates.center;
      fit_result.fit_uses_symmetric_phi_range = fit_uses_symmetric_phi_range;
      fit_result.has_off_diagonal = model == FitModel::kFull;
      fit_result.uses_coulomb = UsesCoulomb(fit_options.coulomb_mode);
      fit_result.coulomb_mode = ToString(fit_options.coulomb_mode);
      fit_result.finite_source_mode = fit_options.coulomb_mode == CoulombMode::kFiniteSource
                                          ? ToString(fit_options.finite_source_mode)
                                          : "";
      fit_result.finite_source_radius_fm = finite_source_radius_fm;
      fit_result.uses_core_halo_lambda = fit_options.use_core_halo_lambda;
      fit_result.uses_q2_baseline = fit_options.use_q2_baseline;
      fit_result.uses_pml = fit_options.use_pml;
      fit_result.norm = fit_function->GetParameter(0);
      fit_result.norm_err = fit_function->GetParError(0);
      fit_result.lambda = fit_options.use_core_halo_lambda ? fit_function->GetParameter(1) : 1.0;
      fit_result.lambda_err = fit_options.use_core_halo_lambda ? fit_function->GetParError(1) : 0.0;
      fit_result.rout2 = fit_function->GetParameter(2);
      fit_result.rout2_err = fit_function->GetParError(2);
      fit_result.rside2 = fit_function->GetParameter(3);
      fit_result.rside2_err = fit_function->GetParError(3);
      fit_result.rlong2 = fit_function->GetParameter(4);
      fit_result.rlong2_err = fit_function->GetParError(4);
      if (model == FitModel::kFull) {
        fit_result.routside2 = fit_function->GetParameter(5);
        fit_result.routside2_err = fit_function->GetParError(5);
        fit_result.routlong2 = fit_function->GetParameter(6);
        fit_result.routlong2_err = fit_function->GetParError(6);
        fit_result.rsidelong2 = fit_function->GetParameter(7);
        fit_result.rsidelong2_err = fit_function->GetParError(7);
        fit_result.alpha = fit_function->GetParameter(8);
        fit_result.alpha_err = fit_function->GetParError(8);
        fit_result.baseline_q2 = fit_options.use_q2_baseline ? fit_function->GetParameter(9) : 0.0;
        fit_result.baseline_q2_err = fit_options.use_q2_baseline ? fit_function->GetParError(9) : 0.0;
      } else {
        fit_result.alpha = fit_function->GetParameter(5);
        fit_result.alpha_err = fit_function->GetParError(5);
        fit_result.baseline_q2 = fit_options.use_q2_baseline ? fit_function->GetParameter(6) : 0.0;
        fit_result.baseline_q2_err = fit_options.use_q2_baseline ? fit_function->GetParError(6) : 0.0;
      }
      fit_result.fit_statistic = fit_statistic;
      fit_result.edm = fit_options.use_pml ? fit_edm : std::numeric_limits<double>::quiet_NaN();
      fit_result.ndf = fit_ndf;
      fit_result.status = fit_status;
      fit_result.minuit_istat = fit_options.use_pml ? fit_minuit_istat : -1;
    }

    std::optional<SingleSliceFitOutput> FitSingleSlice(TH3D *h_cf,
                                                       TH3D *h_se_raw,
                                                       TH3D *h_me_raw,
                                                       const SliceCatalogEntry &entry,
                                                       const FitModel model,
                                                       const LevyFitOptions &fit_options,
                                                       const bool fit_uses_symmetric_phi_range,
                                                       const CoulombKernelTable *coulomb_kernel) {
      if (fit_options.coulomb_mode == CoulombMode::kFiniteSource && coulomb_kernel == nullptr) {
        throw std::runtime_error("Finite-source Coulomb fit has no prepared kernel for group " + entry.group_id + ".");
      }
      ActiveCoulombKernelGuard kernel_guard(coulomb_kernel);
      SingleSliceFitOutput output;
      output.fit_function.reset(model == FitModel::kFull
                                    ? BuildFullLevyFitFunction(entry.slice_id + "_levy3d_full_fit", fit_options)
                                    : BuildLevyFitFunction(entry.slice_id + "_levy3d_fit", fit_options));

      double fit_statistic = 0.0;
      double fit_edm = std::numeric_limits<double>::quiet_NaN();
      int fit_ndf = 0;
      int fit_status = -1;
      int fit_minuit_istat = -1;
      bool fit_succeeded = false;
      if (fit_options.use_pml) {
        fit_succeeded = RunPMLFit(output.fit_function.get(),
                                  h_se_raw,
                                  h_me_raw,
                                  model == FitModel::kFull,
                                  fit_options,
                                  fit_statistic,
                                  fit_ndf,
                                  fit_status,
                                  fit_edm,
                                  fit_minuit_istat);
      } else {
        const auto fit_status_object = h_cf->Fit(output.fit_function.get(), "RSMNQ0");
        fit_statistic = output.fit_function->GetChisquare();
        fit_ndf = output.fit_function->GetNDF();
        fit_status = static_cast<int>(fit_status_object);
        fit_succeeded = true;
      }

      const double finite_source_radius_fm =
          fit_options.coulomb_mode == CoulombMode::kFiniteSource && coulomb_kernel != nullptr
              ? coulomb_kernel->catalog_entry.final_radius_fm
              : std::numeric_limits<double>::quiet_NaN();
      FillFitResultFromFunction(entry,
                                output.fit_function.get(),
                                model,
                                fit_options,
                                fit_uses_symmetric_phi_range,
                                fit_statistic,
                                fit_ndf,
                                fit_status,
                                fit_edm,
                                fit_minuit_istat,
                                finite_source_radius_fm,
                                output.result);
      if (!fit_succeeded) {
        return std::nullopt;
      }
      return std::optional<SingleSliceFitOutput>(std::move(output));
    }

    void WriteSingleSliceFitArtifacts(TH3D *h_cf,
                                      const SliceCatalogEntry &entry,
                                      TF3 *fit_function,
                                      const LevyFitResult &fit_result,
                                      const std::string &fit_root_path,
                                      TFile *shared_output_file,
                                      const CoulombKernelTable *coulomb_kernel) {
      ActiveCoulombKernelGuard kernel_guard(coulomb_kernel);
      BatchModeGuard batch_guard;

      auto *projection_x_data = BuildProjectionXWithinWindow(h_cf, entry.slice_id + "_data_ProjX", kProjection1DWindow);
      auto *projection_y_data = BuildProjectionYWithinWindow(h_cf, entry.slice_id + "_data_ProjY", kProjection1DWindow);
      auto *projection_z_data = BuildProjectionZWithinWindow(h_cf, entry.slice_id + "_data_ProjZ", kProjection1DWindow);

      auto *projection_x_fit =
          BuildProjectionCurveXWithinWindow(fit_function, h_cf, entry.slice_id + "_fit_ProjX", kProjection1DWindow);
      auto *projection_y_fit =
          BuildProjectionCurveYWithinWindow(fit_function, h_cf, entry.slice_id + "_fit_ProjY", kProjection1DWindow);
      auto *projection_z_fit =
          BuildProjectionCurveZWithinWindow(fit_function, h_cf, entry.slice_id + "_fit_ProjZ", kProjection1DWindow);

      auto *canvas_x = BuildProjectionCanvas(
          entry.slice_id + "_canvas_ProjX", projection_x_data, projection_x_fit, "q_{out} (GeV/c)", fit_result);
      auto *canvas_y = BuildProjectionCanvas(
          entry.slice_id + "_canvas_ProjY", projection_y_data, projection_y_fit, "q_{side} (GeV/c)", fit_result);
      auto *canvas_z = BuildProjectionCanvas(
          entry.slice_id + "_canvas_ProjZ", projection_z_data, projection_z_fit, "q_{long} (GeV/c)", fit_result);
      auto *comparison_canvas = Build3DComparisonCanvas(entry.slice_id + "_canvas_3D", h_cf, fit_function, fit_result);

      TFile *output_file = shared_output_file;
      std::unique_ptr<TFile> owned_output_file;
      if (output_file == nullptr) {
        owned_output_file = OpenRootFile(fit_root_path, "UPDATE");
        output_file = owned_output_file.get();
      }

      auto *directory = GetOrCreateDirectoryPath(*output_file, BuildFitDirectory(entry.slice_id));
      directory->cd();

      h_cf->SetName("CF3D");
      fit_function->SetName("LevyFit3D");
      projection_x_data->SetName("Data_ProjX");
      projection_y_data->SetName("Data_ProjY");
      projection_z_data->SetName("Data_ProjZ");
      projection_x_fit->SetName("Fit_ProjX");
      projection_y_fit->SetName("Fit_ProjY");
      projection_z_fit->SetName("Fit_ProjZ");
      canvas_x->SetName("Canvas_ProjX");
      canvas_y->SetName("Canvas_ProjY");
      canvas_z->SetName("Canvas_ProjZ");
      comparison_canvas->SetName("Canvas_3D");

      h_cf->Write("", TObject::kOverwrite);
      fit_function->Write("", TObject::kOverwrite);
      projection_x_data->Write("", TObject::kOverwrite);
      projection_y_data->Write("", TObject::kOverwrite);
      projection_z_data->Write("", TObject::kOverwrite);
      projection_x_fit->Write("", TObject::kOverwrite);
      projection_y_fit->Write("", TObject::kOverwrite);
      projection_z_fit->Write("", TObject::kOverwrite);
      canvas_x->Write("", TObject::kOverwrite);
      canvas_y->Write("", TObject::kOverwrite);
      canvas_z->Write("", TObject::kOverwrite);
      comparison_canvas->Write("", TObject::kOverwrite);
      output_file->cd();

      delete projection_x_data;
      delete projection_y_data;
      delete projection_z_data;
      delete projection_x_fit;
      delete projection_y_fit;
      delete projection_z_fit;
      delete canvas_x;
      delete canvas_y;
      delete canvas_z;
      delete comparison_canvas;
    }

    void WriteFitResultsSummaryTsv(const std::string &path, const std::vector<LevyFitResult> &results) {
      std::ofstream output(path);
      output << "sliceId\tgroupId\tfitModel\tusesCoulomb\tcoulombMode\tfiniteSourceMode\tfiniteSourceRadiusFm"
                "\tusesCoreHaloLambda\tusesQ2Baseline\tusesPML\tcentLow\tcentHigh\tmTLow\tmTHigh"
                "\tqnIndex\tqnLow\tqnHigh\tqnLabel\tisQnIntegrated\tphi\tisPhiIntegrated"
                "\tnorm\tnormErr\tlambda\tlambdaErr\tRout2\tRout2Err\tRside2\tRside2Err"
                "\tRlong2\tRlong2Err\tRoutside2\tRoutside2Err\tRoutlong2\tRoutlong2Err"
                "\tRsidelong2\tRsidelong2Err\talpha\talphaErr\tbaselineQ2\tbaselineQ2Err"
                "\tfitStatistic\tedm\tndf\tstatus\tminuitIstat\tcovarianceQuality\n";
      output << std::fixed << std::setprecision(6);
      for (const LevyFitResult &result : results) {
        output << result.slice_id << "\t" << result.group_id << "\t" << result.fit_model << "\t"
               << (result.uses_coulomb ? 1 : 0) << "\t" << result.coulomb_mode << "\t" << result.finite_source_mode
               << "\t" << result.finite_source_radius_fm << "\t" << (result.uses_core_halo_lambda ? 1 : 0) << "\t"
               << (result.uses_q2_baseline ? 1 : 0) << "\t" << (result.uses_pml ? 1 : 0) << "\t" << result.cent_low
               << "\t" << result.cent_high << "\t" << result.mt_low << "\t" << result.mt_high << "\t"
               << result.qn_index << "\t" << result.qn_low << "\t" << result.qn_high << "\t" << result.qn_label
               << "\t" << (result.is_qn_integrated ? 1 : 0) << "\t" << result.phi
               << "\t" << (result.is_phi_integrated ? 1 : 0) << "\t" << result.norm << "\t" << result.norm_err << "\t"
               << result.lambda << "\t" << result.lambda_err << "\t" << result.rout2 << "\t" << result.rout2_err << "\t"
               << result.rside2 << "\t" << result.rside2_err << "\t" << result.rlong2 << "\t" << result.rlong2_err
               << "\t" << result.routside2 << "\t" << result.routside2_err << "\t" << result.routlong2 << "\t"
               << result.routlong2_err << "\t" << result.rsidelong2 << "\t" << result.rsidelong2_err << "\t"
               << result.alpha << "\t" << result.alpha_err << "\t" << result.baseline_q2 << "\t"
               << result.baseline_q2_err << "\t" << result.fit_statistic << "\t" << result.edm << "\t" << result.ndf
               << "\t" << result.status << "\t" << result.minuit_istat << "\t"
               << CovarianceQualityToken(result.minuit_istat) << "\n";
      }
    }

    void WriteFitCatalogTree(TFile &output_file, const std::vector<LevyFitResult> &results) {
      auto *meta_directory = GetOrCreateDirectoryPath(output_file, "meta");
      meta_directory->cd();

      auto tree = std::make_unique<TTree>("FitCatalog", "FitCatalog");
      std::string slice_id;
      std::string group_id;
      std::string fit_model;
      int centrality_index = -1;
      int mt_index = -1;
      int qn_index = -1;
      int phi_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      double qn_low = std::numeric_limits<double>::quiet_NaN();
      double qn_high = std::numeric_limits<double>::quiet_NaN();
      std::string qn_label;
      double phi = 0.0;
      int fit_uses_symmetric_phi_range = 0;
      int is_qn_integrated = 1;
      int is_phi_integrated = 0;
      double norm = 0.0;
      double norm_err = 0.0;
      double lambda = 0.0;
      double lambda_err = 0.0;
      double rout2 = 0.0;
      double rout2_err = 0.0;
      double rside2 = 0.0;
      double rside2_err = 0.0;
      double rlong2 = 0.0;
      double rlong2_err = 0.0;
      double routside2 = 0.0;
      double routside2_err = 0.0;
      double routlong2 = 0.0;
      double routlong2_err = 0.0;
      double rsidelong2 = 0.0;
      double rsidelong2_err = 0.0;
      double alpha = 0.0;
      double alpha_err = 0.0;
      double baseline_q2 = 0.0;
      double baseline_q2_err = 0.0;
      double fit_statistic = 0.0;
      double edm = 0.0;
      int ndf = 0;
      int status = 0;
      int minuit_istat = 0;
      int uses_coulomb = 0;
      std::string coulomb_mode;
      std::string finite_source_mode;
      double finite_source_radius_fm = 0.0;
      int uses_core_halo_lambda = 0;
      int uses_q2_baseline = 0;
      int uses_pml = 0;
      int has_off_diagonal = 0;

      tree->Branch("slice_id", &slice_id);
      tree->Branch("group_id", &group_id);
      tree->Branch("fit_model", &fit_model);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("qn_index", &qn_index);
      tree->Branch("phi_index", &phi_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("qn_low", &qn_low);
      tree->Branch("qn_high", &qn_high);
      tree->Branch("qn_label", &qn_label);
      tree->Branch("phi", &phi);
      tree->Branch("fit_uses_symmetric_phi_range", &fit_uses_symmetric_phi_range);
      tree->Branch("is_qn_integrated", &is_qn_integrated);
      tree->Branch("is_phi_integrated", &is_phi_integrated);
      tree->Branch("norm", &norm);
      tree->Branch("norm_err", &norm_err);
      tree->Branch("lambda", &lambda);
      tree->Branch("lambda_err", &lambda_err);
      tree->Branch("rout2", &rout2);
      tree->Branch("rout2_err", &rout2_err);
      tree->Branch("rside2", &rside2);
      tree->Branch("rside2_err", &rside2_err);
      tree->Branch("rlong2", &rlong2);
      tree->Branch("rlong2_err", &rlong2_err);
      tree->Branch("routside2", &routside2);
      tree->Branch("routside2_err", &routside2_err);
      tree->Branch("routlong2", &routlong2);
      tree->Branch("routlong2_err", &routlong2_err);
      tree->Branch("rsidelong2", &rsidelong2);
      tree->Branch("rsidelong2_err", &rsidelong2_err);
      tree->Branch("alpha", &alpha);
      tree->Branch("alpha_err", &alpha_err);
      tree->Branch("baseline_q2", &baseline_q2);
      tree->Branch("baseline_q2_err", &baseline_q2_err);
      tree->Branch("fit_statistic", &fit_statistic);
      tree->Branch("edm", &edm);
      tree->Branch("ndf", &ndf);
      tree->Branch("status", &status);
      tree->Branch("minuit_istat", &minuit_istat);
      tree->Branch("uses_coulomb", &uses_coulomb);
      tree->Branch("coulomb_mode", &coulomb_mode);
      tree->Branch("coulombMode", &coulomb_mode);
      tree->Branch("finite_source_mode", &finite_source_mode);
      tree->Branch("finiteSourceMode", &finite_source_mode);
      tree->Branch("finite_source_radius_fm", &finite_source_radius_fm);
      tree->Branch("finiteSourceRadiusFm", &finite_source_radius_fm);
      tree->Branch("uses_core_halo_lambda", &uses_core_halo_lambda);
      tree->Branch("uses_q2_baseline", &uses_q2_baseline);
      tree->Branch("uses_pml", &uses_pml);
      tree->Branch("has_off_diagonal", &has_off_diagonal);

      for (const LevyFitResult &result : results) {
        slice_id = result.slice_id;
        group_id = result.group_id;
        fit_model = result.fit_model;
        centrality_index = result.centrality_index;
        mt_index = result.mt_index;
        qn_index = result.qn_index;
        phi_index = result.phi_index;
        cent_low = result.cent_low;
        cent_high = result.cent_high;
        mt_low = result.mt_low;
        mt_high = result.mt_high;
        qn_low = result.qn_low;
        qn_high = result.qn_high;
        qn_label = result.qn_label;
        phi = result.phi;
        fit_uses_symmetric_phi_range = result.fit_uses_symmetric_phi_range ? 1 : 0;
        is_qn_integrated = result.is_qn_integrated ? 1 : 0;
        is_phi_integrated = result.is_phi_integrated ? 1 : 0;
        norm = result.norm;
        norm_err = result.norm_err;
        lambda = result.lambda;
        lambda_err = result.lambda_err;
        rout2 = result.rout2;
        rout2_err = result.rout2_err;
        rside2 = result.rside2;
        rside2_err = result.rside2_err;
        rlong2 = result.rlong2;
        rlong2_err = result.rlong2_err;
        routside2 = result.routside2;
        routside2_err = result.routside2_err;
        routlong2 = result.routlong2;
        routlong2_err = result.routlong2_err;
        rsidelong2 = result.rsidelong2;
        rsidelong2_err = result.rsidelong2_err;
        alpha = result.alpha;
        alpha_err = result.alpha_err;
        baseline_q2 = result.baseline_q2;
        baseline_q2_err = result.baseline_q2_err;
        fit_statistic = result.fit_statistic;
        edm = result.edm;
        ndf = result.ndf;
        status = result.status;
        minuit_istat = result.minuit_istat;
        uses_coulomb = result.uses_coulomb ? 1 : 0;
        coulomb_mode = result.coulomb_mode;
        finite_source_mode = result.finite_source_mode;
        finite_source_radius_fm = result.finite_source_radius_fm;
        uses_core_halo_lambda = result.uses_core_halo_lambda ? 1 : 0;
        uses_q2_baseline = result.uses_q2_baseline ? 1 : 0;
        uses_pml = result.uses_pml ? 1 : 0;
        has_off_diagonal = result.has_off_diagonal ? 1 : 0;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
      output_file.cd();
    }

    void WriteCoulombKernelCatalogTree(TFile &output_file, const std::vector<CoulombKernelCatalogEntry> &entries) {
      auto *meta_directory = GetOrCreateDirectoryPath(output_file, "meta");
      meta_directory->cd();

      auto tree = std::make_unique<TTree>("CoulombKernelCatalog", "CoulombKernelCatalog");
      std::string group_id;
      int centrality_index = -1;
      int mt_index = -1;
      int qn_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      double qn_low = std::numeric_limits<double>::quiet_NaN();
      double qn_high = std::numeric_limits<double>::quiet_NaN();
      std::string qn_label;
      int is_qn_integrated = 1;
      std::string finite_source_mode;
      double seed_radius_fm = 0.0;
      double final_radius_fm = 0.0;
      int cats_enabled = 0;
      double kstar_min_mev = 0.0;
      double kstar_max_mev = 0.0;
      int kstar_bin_count = 0;

      tree->Branch("group_id", &group_id);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("qn_index", &qn_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("qn_low", &qn_low);
      tree->Branch("qn_high", &qn_high);
      tree->Branch("qn_label", &qn_label);
      tree->Branch("is_qn_integrated", &is_qn_integrated);
      tree->Branch("finite_source_mode", &finite_source_mode);
      tree->Branch("seed_radius_fm", &seed_radius_fm);
      tree->Branch("final_radius_fm", &final_radius_fm);
      tree->Branch("cats_enabled", &cats_enabled);
      tree->Branch("kstar_min_mev", &kstar_min_mev);
      tree->Branch("kstar_max_mev", &kstar_max_mev);
      tree->Branch("kstar_bin_count", &kstar_bin_count);

      for (const CoulombKernelCatalogEntry &entry : entries) {
        group_id = entry.group_id;
        centrality_index = entry.centrality_index;
        mt_index = entry.mt_index;
        qn_index = entry.qn_index;
        cent_low = entry.cent_low;
        cent_high = entry.cent_high;
        mt_low = entry.mt_low;
        mt_high = entry.mt_high;
        qn_low = entry.qn_low;
        qn_high = entry.qn_high;
        qn_label = entry.qn_label;
        is_qn_integrated = entry.is_qn_integrated ? 1 : 0;
        finite_source_mode = entry.finite_source_mode;
        seed_radius_fm = entry.seed_radius_fm;
        final_radius_fm = entry.final_radius_fm;
        cats_enabled = entry.cats_enabled ? 1 : 0;
        kstar_min_mev = entry.kstar_min_mev;
        kstar_max_mev = entry.kstar_max_mev;
        kstar_bin_count = entry.kstar_bin_count;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
      output_file.cd();
    }

    void WriteR2Graphs(TFile &output_file, const std::vector<LevyFitResult> &results) {
      std::map<std::string, std::vector<LevyFitResult>> grouped_results;
      for (const LevyFitResult &result : results) {
        if (result.is_phi_integrated) {
          continue;
        }
        grouped_results[result.group_id].push_back(result);
      }

      for (auto &[group_id, group_results] : grouped_results) {
        if (group_results.empty()) {
          continue;
        }
        std::sort(group_results.begin(), group_results.end(), [](const LevyFitResult &lhs, const LevyFitResult &rhs) {
          return lhs.phi < rhs.phi;
        });

        auto *directory = GetOrCreateDirectoryPath(output_file, "summary/R2_vs_phi/" + group_id);
        directory->cd();

        const bool uses_mapped_phi_range = group_results.front().fit_uses_symmetric_phi_range;
        const double phi_fit_min = uses_mapped_phi_range ? -TMath::Pi() / 2.0 : 0.0;
        const double phi_fit_max = uses_mapped_phi_range ? TMath::Pi() / 2.0 : TMath::Pi();
        const int n_points = static_cast<int>(group_results.size());
        const bool has_off_diagonal =
            std::any_of(group_results.begin(), group_results.end(), [](const LevyFitResult &result) {
              return result.has_off_diagonal;
            });
        const bool uses_core_halo_lambda =
            std::any_of(group_results.begin(), group_results.end(), [](const LevyFitResult &result) {
              return result.uses_core_halo_lambda;
            });
        const bool uses_q2_baseline =
            std::any_of(group_results.begin(), group_results.end(), [](const LevyFitResult &result) {
              return result.uses_q2_baseline;
            });

        auto *g_rout2 = new TGraphErrors(n_points);
        auto *g_rside2 = new TGraphErrors(n_points);
        auto *g_rlong2 = new TGraphErrors(n_points);
        TGraphErrors *g_routside2 = has_off_diagonal ? new TGraphErrors(n_points) : nullptr;
        TGraphErrors *g_routlong2 = has_off_diagonal ? new TGraphErrors(n_points) : nullptr;
        TGraphErrors *g_rsidelong2 = has_off_diagonal ? new TGraphErrors(n_points) : nullptr;
        auto *g_alpha = new TGraphErrors(n_points);
        TGraphErrors *g_lambda = uses_core_halo_lambda ? new TGraphErrors(n_points) : nullptr;
        TGraphErrors *g_baseline_q2 = uses_q2_baseline ? new TGraphErrors(n_points) : nullptr;

        for (int index = 0; index < n_points; ++index) {
          const LevyFitResult &result = group_results[index];
          g_rout2->SetPoint(index, result.phi, result.rout2);
          g_rout2->SetPointError(index, 0.0, result.rout2_err);
          g_rside2->SetPoint(index, result.phi, result.rside2);
          g_rside2->SetPointError(index, 0.0, result.rside2_err);
          g_rlong2->SetPoint(index, result.phi, result.rlong2);
          g_rlong2->SetPointError(index, 0.0, result.rlong2_err);
          g_alpha->SetPoint(index, result.phi, result.alpha);
          g_alpha->SetPointError(index, 0.0, result.alpha_err);

          if (has_off_diagonal) {
            g_routside2->SetPoint(index, result.phi, result.routside2);
            g_routside2->SetPointError(index, 0.0, result.routside2_err);
            g_routlong2->SetPoint(index, result.phi, result.routlong2);
            g_routlong2->SetPointError(index, 0.0, result.routlong2_err);
            g_rsidelong2->SetPoint(index, result.phi, result.rsidelong2);
            g_rsidelong2->SetPointError(index, 0.0, result.rsidelong2_err);
          }
          if (uses_core_halo_lambda) {
            g_lambda->SetPoint(index, result.phi, result.lambda);
            g_lambda->SetPointError(index, 0.0, result.lambda_err);
          }
          if (uses_q2_baseline) {
            g_baseline_q2->SetPoint(index, result.phi, result.baseline_q2);
            g_baseline_q2->SetPointError(index, 0.0, result.baseline_q2_err);
          }
        }

        auto *fit_cos_rout2 = new TF1("Rout2PhiFit", "[0]+2.0*[1]*cos(2.0*x)", phi_fit_min, phi_fit_max);
        auto *fit_cos_rside2 = new TF1("Rside2PhiFit", "[0]+2.0*[1]*cos(2.0*x)", phi_fit_min, phi_fit_max);
        auto *fit_cos_rlong2 = new TF1("Rlong2PhiFit", "[0]+2.0*[1]*cos(2.0*x)", phi_fit_min, phi_fit_max);
        fit_cos_rout2->SetParameters(group_results.front().rout2, 0.0);
        fit_cos_rside2->SetParameters(group_results.front().rside2, 0.0);
        fit_cos_rlong2->SetParameters(group_results.front().rlong2, 0.0);

        g_rout2->SetName("Rout2_vs_phi");
        g_rside2->SetName("Rside2_vs_phi");
        g_rlong2->SetName("Rlong2_vs_phi");
        g_rout2->Fit(fit_cos_rout2, "QN");
        g_rside2->Fit(fit_cos_rside2, "QN");
        g_rlong2->Fit(fit_cos_rlong2, "QN");
        g_rout2->Write("", TObject::kOverwrite);
        g_rside2->Write("", TObject::kOverwrite);
        g_rlong2->Write("", TObject::kOverwrite);
        fit_cos_rout2->Write("Rout2_phi_fit", TObject::kOverwrite);
        fit_cos_rside2->Write("Rside2_phi_fit", TObject::kOverwrite);
        fit_cos_rlong2->Write("Rlong2_phi_fit", TObject::kOverwrite);

        if (has_off_diagonal) {
          auto *fit_sin_routside2 = new TF1("Routside2PhiFit", "[0]+2.0*[1]*sin(2.0*x)", phi_fit_min, phi_fit_max);
          auto *fit_cos_routlong2 = new TF1("Routlong2PhiFit", "[0]+2.0*[1]*cos(2.0*x)", phi_fit_min, phi_fit_max);
          auto *fit_sin_rsidelong2 = new TF1("Rsidelong2PhiFit", "[0]+2.0*[1]*sin(2.0*x)", phi_fit_min, phi_fit_max);
          fit_sin_routside2->SetParameters(group_results.front().routside2, 0.0);
          fit_cos_routlong2->SetParameters(group_results.front().routlong2, 0.0);
          fit_sin_rsidelong2->SetParameters(group_results.front().rsidelong2, 0.0);
          g_routside2->SetName("Routside2_vs_phi");
          g_routlong2->SetName("Routlong2_vs_phi");
          g_rsidelong2->SetName("Rsidelong2_vs_phi");
          g_routside2->Fit(fit_sin_routside2, "QN");
          g_routlong2->Fit(fit_cos_routlong2, "QN");
          g_rsidelong2->Fit(fit_sin_rsidelong2, "QN");
          g_routside2->Write("", TObject::kOverwrite);
          g_routlong2->Write("", TObject::kOverwrite);
          g_rsidelong2->Write("", TObject::kOverwrite);
          fit_sin_routside2->Write("Routside2_phi_fit", TObject::kOverwrite);
          fit_cos_routlong2->Write("Routlong2_phi_fit", TObject::kOverwrite);
          fit_sin_rsidelong2->Write("Rsidelong2_phi_fit", TObject::kOverwrite);
          delete fit_sin_routside2;
          delete fit_cos_routlong2;
          delete fit_sin_rsidelong2;
          delete g_routside2;
          delete g_routlong2;
          delete g_rsidelong2;
        }

        auto *fit_const_alpha = new TF1("AlphaPhiFit", "[0]", phi_fit_min, phi_fit_max);
        fit_const_alpha->SetParameter(0, group_results.front().alpha);
        g_alpha->SetName("alpha_vs_phi");
        g_alpha->Fit(fit_const_alpha, "QN");
        g_alpha->Write("", TObject::kOverwrite);
        fit_const_alpha->Write("alpha_phi_fit", TObject::kOverwrite);
        delete fit_const_alpha;

        if (uses_core_halo_lambda) {
          auto *fit_const_lambda = new TF1("LambdaPhiFit", "[0]", phi_fit_min, phi_fit_max);
          fit_const_lambda->SetParameter(0, group_results.front().lambda);
          g_lambda->SetName("lambda_vs_phi");
          g_lambda->Fit(fit_const_lambda, "QN");
          g_lambda->Write("", TObject::kOverwrite);
          fit_const_lambda->Write("lambda_phi_fit", TObject::kOverwrite);
          delete fit_const_lambda;
          delete g_lambda;
        }
        if (uses_q2_baseline) {
          auto *fit_const_baseline = new TF1("BaselineQ2PhiFit", "[0]", phi_fit_min, phi_fit_max);
          fit_const_baseline->SetParameter(0, group_results.front().baseline_q2);
          g_baseline_q2->SetName("baselineQ2_vs_phi");
          g_baseline_q2->Fit(fit_const_baseline, "QN");
          g_baseline_q2->Write("", TObject::kOverwrite);
          fit_const_baseline->Write("baselineQ2_phi_fit", TObject::kOverwrite);
          delete fit_const_baseline;
          delete g_baseline_q2;
        }

        delete g_rout2;
        delete g_rside2;
        delete g_rlong2;
        delete g_alpha;
        delete fit_cos_rout2;
        delete fit_cos_rside2;
        delete fit_cos_rlong2;
        output_file.cd();
      }
    }

    // Write one publication-style source-parameter canvas for each fitted cent/mT group.
    void WriteSourceParameterOverviewCanvases(TDirectory &directory, const std::vector<LevyFitResult> &results) {
      const GroupedFitResults grouped_results = GroupPhiDifferentialResultsByCentMt(results);
      const std::array<std::string, 6> graph_names = {
          "Rout2_vs_phi", "Rside2_vs_phi", "Rlong2_vs_phi", "Routside2_vs_phi", "Routlong2_vs_phi",
          "Rsidelong2_vs_phi"};
      const std::array<std::string, 6> graph_titles = {"R_{out}^{2}",
                                                       "R_{side}^{2}",
                                                       "R_{long}^{2}",
                                                       "R_{outside}^{2}",
                                                       "R_{outlong}^{2}",
                                                       "R_{sidelong}^{2}"};
      const std::string relative_phi_title = BuildRelativePhiAxisTitle();

      for (const auto &entry : grouped_results) {
        const std::vector<LevyFitResult> &group_results = entry.second;
        if (group_results.empty()) {
          continue;
        }

        bool has_drawable_panel = HasValidSummaryPoints(BuildAlphaSummaryPoints(group_results));
        for (std::size_t radius_index = 0; radius_index < graph_names.size(); ++radius_index) {
          has_drawable_panel =
              has_drawable_panel || HasValidSummaryPoints(BuildRadiusSummaryPoints(group_results, radius_index));
        }
        if (!has_drawable_panel) {
          continue;
        }

        TDirectory *cent_directory =
            GetOrCreateDirectory(directory, BuildReportCentralityDirectory(group_results.front()));
        TDirectory *mt_directory = GetOrCreateDirectory(*cent_directory, BuildReportMtDirectory(group_results.front()));
        TDirectory *target_directory = mt_directory;
        const std::string qn_directory_name = BuildReportQnDirectory(group_results.front());
        if (!qn_directory_name.empty()) {
          target_directory = GetOrCreateDirectory(*mt_directory, qn_directory_name);
        }
        target_directory->cd();

        auto canvas =
            std::make_unique<TCanvas>("source_parameters_overview_canvas", "Source parameter overview", 1200, 1800);
        canvas->Divide(2, 4, 0.0, 0.0);
        std::vector<std::unique_ptr<TGraphErrors>> owned_graphs;
        std::vector<std::unique_ptr<TF1>> owned_fits;

        TVirtualPad *info_pad = canvas->cd(1);
        if (info_pad != nullptr) {
          info_pad->SetLeftMargin(0.04);
          info_pad->SetRightMargin(0.04);
          info_pad->SetBottomMargin(0.04);
          info_pad->SetTopMargin(0.04);
          CreateSourceParameterOverviewInfoBox(group_results.front())->Draw();
        }

        DrawOverviewGraph(*canvas,
                          2,
                          MakePhiSummaryGraph(BuildAlphaSummaryPoints(group_results),
                                              "alpha_vs_phi",
                                              "#alpha vs " + relative_phi_title,
                                              relative_phi_title,
                                              "#alpha"),
                          -1,
                          group_results.front().fit_uses_symmetric_phi_range,
                          owned_graphs,
                          owned_fits);

        // Keep the same visual panel order as the Eventgen overview: alpha first,
        // then out/cross, side/cross, and long/cross rows.
        constexpr std::array<std::size_t, 6> kOverviewRadiusOrder = {0U, 3U, 1U, 4U, 2U, 5U};
        for (std::size_t order_index = 0; order_index < kOverviewRadiusOrder.size(); ++order_index) {
          const std::size_t radius_index = kOverviewRadiusOrder[order_index];
          DrawOverviewGraph(*canvas,
                            static_cast<int>(order_index + 3U),
                            MakePhiSummaryGraph(BuildRadiusSummaryPoints(group_results, radius_index),
                                                graph_names[radius_index],
                                                graph_titles[radius_index] + " vs " + relative_phi_title,
                                                relative_phi_title,
                                                graph_titles[radius_index] + " [fm^{2}]"),
                            static_cast<int>(radius_index),
                            group_results.front().fit_uses_symmetric_phi_range,
                            owned_graphs,
                            owned_fits);
        }

        canvas->Modified();
        canvas->Update();
        canvas->Write("", TObject::kOverwrite);
        directory.cd();
      }
    }

    // Write one epsilon_f(mT) summary graph per centrality from the side-radius harmonic fits.
    void WriteEpsVsMtGraphs(TDirectory &directory, const std::vector<LevyFitResult> &results) {
      EpsSummaryMap eps_summary_points;
      std::map<std::pair<int, int>, std::string> centrality_directory_names;
      std::map<std::pair<int, int>, std::string> qn_directory_names;
      const GroupedFitResults grouped_results = GroupPhiDifferentialResultsByCentMt(results);
      for (const auto &entry : grouped_results) {
        const std::vector<LevyFitResult> &group_results = entry.second;
        if (group_results.empty()) {
          continue;
        }
        constexpr std::size_t kSideRadiusIndex = 1U;
        const std::optional<EpsSummaryPoint> eps_point =
            ComputeEpsFromRsideSummaryPoints(BuildRadiusSummaryPoints(group_results, kSideRadiusIndex),
                                             group_results.front());
        if (eps_point.has_value() && eps_point->valid) {
          const auto summary_key = std::make_pair(group_results.front().centrality_index, group_results.front().qn_index);
          eps_summary_points[summary_key].push_back(*eps_point);
          centrality_directory_names[summary_key] = BuildReportCentralityDirectory(group_results.front());
          qn_directory_names[summary_key] = BuildReportQnDirectory(group_results.front());
        }
      }

      for (auto &entry : eps_summary_points) {
        auto &points = entry.second;
        std::sort(points.begin(), points.end(), [](const EpsSummaryPoint &left, const EpsSummaryPoint &right) {
          return left.mt_center < right.mt_center;
        });
        const auto valid_count =
            static_cast<int>(std::count_if(points.begin(), points.end(), [](const EpsSummaryPoint &point) {
              return point.valid && std::isfinite(point.mt_center) && std::isfinite(point.value);
            }));
        if (valid_count == 0) {
          continue;
        }

        TDirectory *cent_directory = GetOrCreateDirectory(directory, centrality_directory_names[entry.first]);
        TDirectory *target_directory = cent_directory;
        const std::string qn_directory_name = qn_directory_names[entry.first];
        if (!qn_directory_name.empty()) {
          target_directory = GetOrCreateDirectory(*cent_directory, qn_directory_name);
        }
        auto graph = std::make_unique<TGraphErrors>(valid_count);
        graph->SetName("epsf_vs_mt");
        graph->SetTitle("#epsilon_{f} vs m_{T}");
        graph->GetXaxis()->SetTitle("m_{T} [GeV/c^{2}]");
        graph->GetYaxis()->SetTitle("#epsilon_{f} = 2R_{s,2}^{2}/R_{s,0}^{2}");
        ApplyReportGraphStyle(*graph);

        int point_index = 0;
        for (const EpsSummaryPoint &point : points) {
          if (!point.valid || !std::isfinite(point.mt_center) || !std::isfinite(point.value)) {
            continue;
          }
          graph->SetPoint(point_index, point.mt_center, point.value);
          graph->SetPointError(point_index,
                               std::isfinite(point.mt_error) ? point.mt_error : 0.0,
                               std::isfinite(point.error) ? point.error : 0.0);
          ++point_index;
        }

        target_directory->cd();
        auto canvas = std::make_unique<TCanvas>("epsf_vs_mt_canvas", "#epsilon_{f} vs m_{T}", 800, 600);
        canvas->SetTicks(1, 1);
        graph->Draw("ALPE1");
        CreateSummaryInfoBox("epsf_vs_mt", ComputeEpsGraphStats(points))->Draw();
        canvas->Modified();
        canvas->Update();
        graph->Write("", TObject::kOverwrite);
        canvas->Write("", TObject::kOverwrite);
        directory.cd();
      }
    }

    // The standalone report file is independent of the detailed per-slice fit ROOT file.
    void WriteFitReportRootFile(const std::string &fit_report_root_path,
                                const std::vector<LevyFitResult> &fit_results,
                                const std::vector<CoulombKernelCatalogEntry> &kernel_catalog_entries) {
      CreateOrResetRootFile(fit_report_root_path);
      auto output_file = OpenRootFile(fit_report_root_path, "UPDATE");

      // The report file is summary-only: it mirrors the legacy summary products
      // and adds canvases grouped by centrality and mT.
      WriteFitCatalogTree(*output_file, fit_results);
      WriteCoulombKernelCatalogTree(*output_file, kernel_catalog_entries);
      WriteR2Graphs(*output_file, fit_results);
      auto *source_parameters_directory = GetOrCreateDirectoryPath(*output_file, "source_parameters");
      WriteSourceParameterOverviewCanvases(*source_parameters_directory, fit_results);
      auto *eps_directory = GetOrCreateDirectoryPath(*output_file, "eps_vs_mt");
      WriteEpsVsMtGraphs(*eps_directory, fit_results);
      output_file->Write();
      output_file->Close();
    }

  }  // namespace

  BuildCfRunStatistics RunBuildCf(const ApplicationConfig &config, const Logger &logger) {
    EnsureDirectoryExists(config.output.output_directory);
    const std::string cf_root_path = ResolvePath(config.output.output_directory, config.output.cf_root_name);
    CreateOrResetRootFile(cf_root_path);

    const bool old_add_directory = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    auto input_file = OpenRootFile(config.input.input_root, "READ");
    const std::string se_path = BuildSparseObjectPath(config, config.input.same_event_subtask);
    const std::string me_path = BuildSparseObjectPath(config, config.input.mixed_event_subtask);

    auto *h_se_origin = dynamic_cast<THnSparseF *>(input_file->Get(se_path.c_str()));
    auto *h_me_origin = dynamic_cast<THnSparseF *>(input_file->Get(me_path.c_str()));
    if (h_se_origin == nullptr || h_me_origin == nullptr) {
      TH1::AddDirectory(old_add_directory);
      throw std::runtime_error("Cannot resolve required THnSparse inputs from ROOT file.");
    }

    h_se_origin->Sumw2();
    h_me_origin->Sumw2();

    TAxis *phi_axis = h_se_origin->GetAxis(6);
    const int n_phi_bins = phi_axis->GetNbins();
    TAxis *qn_axis = h_se_origin->GetAxis(5);
    if (qn_axis == nullptr || h_me_origin->GetAxis(5) == nullptr) {
      TH1::AddDirectory(old_add_directory);
      throw std::runtime_error("Cannot resolve qn axis 5 from required THnSparse inputs.");
    }
    const std::vector<QnSliceSelection> qn_slice_selections = BuildQnSliceSelections(config, *qn_axis);

    BuildCfRunStatistics statistics;
    statistics.requested_groups = config.centrality_bins.size() * config.mt_bins.size() * qn_slice_selections.size();
    std::vector<SliceCatalogEntry> catalog_entries;
    // Count planned slices up front so skipped branches still advance the CLI progress bar.
    const std::size_t slices_per_group = static_cast<std::size_t>(n_phi_bins + 1);
    const std::size_t planned_slices = statistics.requested_groups * slices_per_group;
    std::size_t processed_slices = 0;
    ProgressReporter progress(logger, "build-cf", planned_slices, config.build.progress);

    std::unique_ptr<TFile> shared_output_file;
    if (!config.build.reopen_output_file_per_slice) {
      shared_output_file = OpenRootFile(cf_root_path, "UPDATE");
    }

    logger.Info("Starting build-cf stage.");
    progress.Update(0);

    for (std::size_t centrality_index = 0; centrality_index < config.centrality_bins.size(); ++centrality_index) {
      const RangeBin &centrality_bin = config.centrality_bins[centrality_index];
      for (std::size_t mt_index = 0; mt_index < config.mt_bins.size(); ++mt_index) {
        const RangeBin &mt_bin = config.mt_bins[mt_index];
        const std::string base_group_id = BuildBaseGroupId(centrality_bin, mt_bin);
        logger.Debug("Building base group " + base_group_id);

        h_se_origin->GetAxis(4)->SetRangeUser(centrality_bin.min, centrality_bin.max);
        h_se_origin->GetAxis(3)->SetRangeUser(mt_bin.min, mt_bin.max);
        h_se_origin->GetAxis(6)->SetRange(1, n_phi_bins);
        ResetAxisVisibleRange(*h_se_origin->GetAxis(5));
        h_me_origin->GetAxis(4)->SetRangeUser(centrality_bin.min, centrality_bin.max);
        h_me_origin->GetAxis(3)->SetRangeUser(mt_bin.min, mt_bin.max);
        ResetAxisVisibleRange(*h_me_origin->GetAxis(5));

        // ME denominator splitting is independently opt-in in phi and qn, so the projection helper applies
        // the requested qn policy before each raw 3D projection.
        auto build_me_projection = [&](const int first_phi_bin,
                                       const int last_phi_bin,
                                       const QnSliceSelection &qn_selection,
                                       const std::string &norm_name) -> std::unique_ptr<MixedEventProjection> {
          h_me_origin->GetAxis(6)->SetRange(first_phi_bin, last_phi_bin);
          if (config.build.split_mixed_event_by_qn) {
            ApplyQnSelection(*h_me_origin, qn_selection);
          } else {
            ResetAxisVisibleRange(*h_me_origin->GetAxis(5));
          }
          auto raw = std::unique_ptr<TH3D>(static_cast<TH3D *>(h_me_origin->Projection(0, 1, 2)));
          raw->SetDirectory(nullptr);
          auto norm = std::unique_ptr<TH3D>(static_cast<TH3D *>(raw->Clone(norm_name.c_str())));
          norm->SetDirectory(nullptr);
          const double integral = IntegralVisibleRange(norm.get(), true);
          if (integral == 0.0) {
            return nullptr;
          }
          norm->Scale(1.0 / integral);
          auto projection = std::make_unique<MixedEventProjection>();
          projection->raw = std::move(raw);
          projection->norm = std::move(norm);
          return projection;
        };

        std::unique_ptr<MixedEventProjection> integrated_me_projection;
        if (!config.build.split_mixed_event_by_phi && !config.build.split_mixed_event_by_qn) {
          integrated_me_projection =
              build_me_projection(1, n_phi_bins, qn_slice_selections.front(), base_group_id + "_ME_norm");
          if (integrated_me_projection == nullptr) {
            logger.Warn("Zero mixed-event integral for " + base_group_id + "; skipping group.");
            statistics.skipped_zero_mixed_event_groups += qn_slice_selections.size();
            processed_slices += slices_per_group * qn_slice_selections.size();
            progress.Update(processed_slices);
            continue;
          }
        }

        auto write_slice = [&](TH3D *h_se_raw,
                               TH3D *h_se_norm,
                               const MixedEventProjection &me_projection,
                               const SliceCatalogEntry &entry) {
          auto *h_cf = static_cast<TH3D *>(h_se_norm->Clone("CF3D"));
          h_cf->SetDirectory(nullptr);
          h_cf->Divide(me_projection.norm.get());
          h_cf->GetXaxis()->SetTitle("q_{out} (GeV/c)");
          h_cf->GetYaxis()->SetTitle("q_{side} (GeV/c)");
          h_cf->GetZaxis()->SetTitle("q_{long} (GeV/c)");

          auto *h_me_write = static_cast<TH3D *>(me_projection.raw->Clone("ME_raw3d"));
          h_me_write->SetDirectory(nullptr);
          h_se_raw->SetName("SE_raw3d");
          h_cf->SetName("CF3D");

          TFile *output_file = shared_output_file.get();
          std::unique_ptr<TFile> owned_output_file;
          if (output_file == nullptr) {
            owned_output_file = OpenRootFile(cf_root_path, "UPDATE");
            output_file = owned_output_file.get();
          }

          auto *directory = GetOrCreateDirectoryPath(*output_file, entry.slice_directory);
          directory->cd();
          directory->WriteObject(h_se_raw, "SE_raw3d");
          directory->WriteObject(h_me_write, "ME_raw3d");
          directory->WriteObject(h_cf, "CF3D");
          Write1DProjections(h_cf, *directory, "CF3D", "C(q)", true);
          if (config.build.write_normalized_se_me_1d_projections) {
            Write1DProjections(h_se_norm, *directory, "SE_norm3d", "Normalized density", true);
            Write1DProjections(me_projection.norm.get(), *directory, "ME_norm3d", "Normalized density", true);
          }
          output_file->cd();
          catalog_entries.push_back(entry);
          ++statistics.stored_slices;
          ++processed_slices;
          progress.Update(processed_slices);

          delete h_me_write;
          delete h_cf;
        };

        auto write_qn_selection = [&](const QnSliceSelection &qn_selection) {
          const std::string group_id = BuildGroupId(centrality_bin, mt_bin, qn_selection);
          logger.Debug("Building group " + group_id);
          ApplyQnSelection(*h_se_origin, qn_selection);

          std::unique_ptr<MixedEventProjection> qn_integrated_me_projection;
          const MixedEventProjection *integrated_me_for_selection = integrated_me_projection.get();
          if (!config.build.split_mixed_event_by_phi && config.build.split_mixed_event_by_qn) {
            qn_integrated_me_projection = build_me_projection(1, n_phi_bins, qn_selection, group_id + "_ME_norm");
            if (qn_integrated_me_projection == nullptr) {
              logger.Warn("Zero mixed-event integral for " + group_id + "; skipping group.");
              ++statistics.skipped_zero_mixed_event_groups;
              processed_slices += slices_per_group;
              progress.Update(processed_slices);
              return;
            }
            integrated_me_for_selection = qn_integrated_me_projection.get();
          }

          h_se_origin->GetAxis(6)->SetRange(1, n_phi_bins);
          auto *h_se_all_raw = static_cast<TH3D *>(h_se_origin->Projection(0, 1, 2));
          h_se_all_raw->SetDirectory(nullptr);
          auto *h_se_all_norm = static_cast<TH3D *>(h_se_all_raw->Clone((group_id + "_SE_all_norm").c_str()));
          h_se_all_norm->SetDirectory(nullptr);
          const double int_se_all = IntegralVisibleRange(h_se_all_norm, true);
          if (int_se_all == 0.0) {
            logger.Warn("Zero same-event integral for " + group_id + " phi=all.");
            ++statistics.skipped_zero_same_event_slices;
            ++processed_slices;
            progress.Update(processed_slices);
          } else {
            h_se_all_norm->Scale(1.0 / int_se_all);
            std::unique_ptr<MixedEventProjection> split_me_projection;
            const MixedEventProjection *me_projection = integrated_me_for_selection;
            if (config.build.split_mixed_event_by_phi) {
              split_me_projection = build_me_projection(1, n_phi_bins, qn_selection, group_id + "_ME_all_norm");
              if (split_me_projection == nullptr) {
                logger.Warn("Zero mixed-event integral for " + group_id + " phi=all.");
                ++statistics.skipped_zero_mixed_event_slices;
                ++processed_slices;
                progress.Update(processed_slices);
                me_projection = nullptr;
              } else {
                me_projection = split_me_projection.get();
              }
            }
            if (me_projection != nullptr) {
              const double raw_phi_low = phi_axis->GetBinLowEdge(1);
              const double raw_phi_high = phi_axis->GetBinUpEdge(phi_axis->GetNbins());
              const auto entry = MakeSliceCatalogEntry(centrality_bin,
                                                       mt_bin,
                                                       qn_selection,
                                                       static_cast<int>(centrality_index),
                                                       static_cast<int>(mt_index),
                                                       -1,
                                                       raw_phi_low,
                                                       raw_phi_high,
                                                       0.5 * (raw_phi_low + raw_phi_high),
                                                       raw_phi_low,
                                                       raw_phi_high,
                                                       std::numeric_limits<double>::quiet_NaN(),
                                                       config.build.map_pair_phi_to_symmetric_range,
                                                       config.build.split_mixed_event_by_phi,
                                                       config.build.split_mixed_event_by_qn,
                                                       true);
              write_slice(h_se_all_raw, h_se_all_norm, *me_projection, entry);
            }
          }
          delete h_se_all_raw;
          delete h_se_all_norm;

          for (int phi_index = 1; phi_index <= n_phi_bins; ++phi_index) {
            h_se_origin->GetAxis(6)->SetRange(phi_index, phi_index);
            auto *h_se_raw = static_cast<TH3D *>(h_se_origin->Projection(0, 1, 2));
            h_se_raw->SetDirectory(nullptr);
            auto *h_se_norm = static_cast<TH3D *>(h_se_raw->Clone((group_id + "_SE_slice_norm").c_str()));
            h_se_norm->SetDirectory(nullptr);
            const double int_se = IntegralVisibleRange(h_se_norm, true);
            if (int_se == 0.0) {
              logger.Warn("Zero same-event integral for " + group_id + " phi bin " + std::to_string(phi_index)
                          + "; skipping slice.");
              ++statistics.skipped_zero_same_event_slices;
              ++processed_slices;
              progress.Update(processed_slices);
              delete h_se_raw;
              delete h_se_norm;
              continue;
            }
            h_se_norm->Scale(1.0 / int_se);

            std::unique_ptr<MixedEventProjection> split_me_projection;
            const MixedEventProjection *me_projection = integrated_me_for_selection;
            if (config.build.split_mixed_event_by_phi) {
              split_me_projection = build_me_projection(phi_index,
                                                        phi_index,
                                                        qn_selection,
                                                        group_id + "_ME_phi" + std::to_string(phi_index) + "_norm");
              if (split_me_projection == nullptr) {
                logger.Warn("Zero mixed-event integral for " + group_id + " phi bin " + std::to_string(phi_index)
                            + "; skipping slice.");
                ++statistics.skipped_zero_mixed_event_slices;
                ++processed_slices;
                progress.Update(processed_slices);
                delete h_se_raw;
                delete h_se_norm;
                continue;
              }
              me_projection = split_me_projection.get();
            }

            const double raw_phi_low = phi_axis->GetBinLowEdge(phi_index);
            const double raw_phi_high = phi_axis->GetBinUpEdge(phi_index);
            const double raw_phi_center = phi_axis->GetBinCenter(phi_index);
            // Build records both raw and display phi so fit can later follow or override
            // the original mapping choice without rebuilding the CF histograms.
            const PhiSliceCoordinates display_phi_coordinates = BuildPhiSliceCoordinatesFromRaw(
                raw_phi_low, raw_phi_high, raw_phi_center, config.build.map_pair_phi_to_symmetric_range);

            const auto entry = MakeSliceCatalogEntry(centrality_bin,
                                                     mt_bin,
                                                     qn_selection,
                                                     static_cast<int>(centrality_index),
                                                     static_cast<int>(mt_index),
                                                     phi_index - 1,
                                                     raw_phi_low,
                                                     raw_phi_high,
                                                     raw_phi_center,
                                                     display_phi_coordinates.low,
                                                     display_phi_coordinates.high,
                                                     display_phi_coordinates.center,
                                                     config.build.map_pair_phi_to_symmetric_range,
                                                     config.build.split_mixed_event_by_phi,
                                                     config.build.split_mixed_event_by_qn,
                                                     false);
            write_slice(h_se_raw, h_se_norm, *me_projection, entry);
            delete h_se_raw;
            delete h_se_norm;
          }
        };

        for (const QnSliceSelection &qn_selection : qn_slice_selections) {
          write_qn_selection(qn_selection);
        }

        h_se_origin->GetAxis(6)->SetRange(1, n_phi_bins);
        h_me_origin->GetAxis(6)->SetRange(1, n_phi_bins);
        ResetAxisVisibleRange(*h_se_origin->GetAxis(5));
        ResetAxisVisibleRange(*h_me_origin->GetAxis(5));
      }
    }

    if (shared_output_file) {
      WriteSliceCatalogTree(*shared_output_file, catalog_entries);
      shared_output_file->Close();
      shared_output_file.reset();
    } else {
      auto catalog_output_file = OpenRootFile(cf_root_path, "UPDATE");
      WriteSliceCatalogTree(*catalog_output_file, catalog_entries);
      catalog_output_file->Close();
    }

    progress.Finish();
    TH1::AddDirectory(old_add_directory);
    logger.Info("Completed build-cf stage: stored " + std::to_string(statistics.stored_slices) + " slices.");
    return statistics;
  }

  FitRunStatistics RunFit(const ApplicationConfig &config,
                          const Logger &logger,
                          const std::optional<FitModel> override_model,
                          const std::optional<std::string> input_cf_root_path) {
    EnsureDirectoryExists(config.output.output_directory);
    EnsureDirectoryExists(config.output.fit_report_directory);
    const std::string cf_root_path = input_cf_root_path.has_value()
                                         ? ResolvePath(config.output.output_directory, *input_cf_root_path)
                                         : ResolvePath(config.output.output_directory, config.output.cf_root_name);
    const std::string fit_root_path = ResolvePath(config.output.output_directory, config.output.fit_root_name);
    const std::string fit_summary_path = ResolvePath(config.output.output_directory, config.output.fit_summary_name);
    const std::string fit_report_root_path =
        ResolvePath(config.output.fit_report_directory, config.output.fit_report_root_name);
    if (PathsReferToSameLocation(fit_report_root_path, fit_root_path)
        || PathsReferToSameLocation(fit_report_root_path, cf_root_path)) {
      throw std::runtime_error(
          "output.fit_report_root_name must resolve to a file distinct from the CF and detailed fit ROOT files.");
    }
    if (config.fit.options.coulomb_mode == CoulombMode::kFiniteSource && !HasCatsFiniteSourceSupport()) {
      throw std::runtime_error(
          "fit.coulomb_mode=\"finite_source\" requires CATS support; reconfigure with EXP_FEMTO_3D_HAS_CATS=1.");
    }

    CreateOrResetRootFile(fit_root_path);
    const bool old_add_directory = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    auto input_file = OpenRootFile(cf_root_path, "READ");
    const std::vector<SliceCatalogEntry> catalog_entries = ReadSliceCatalogTree(*input_file, &logger);
    const bool input_cf_uses_symmetric_phi_range =
        !catalog_entries.empty() && catalog_entries.front().build_uses_symmetric_phi_range;
    if (!catalog_entries.empty() && input_cf_uses_symmetric_phi_range != config.build.map_pair_phi_to_symmetric_range) {
      logger.Warn("Input CF build phi mapping metadata does not match config.build.map_pair_phi_to_symmetric_range;"
                  " using the CF metadata as the build-stage truth.");
    }
    const bool fit_uses_symmetric_phi_range =
        config.fit.map_pair_phi_to_symmetric_range.value_or(input_cf_uses_symmetric_phi_range);

    std::vector<const SliceCatalogEntry *> selected_entries;
    selected_entries.reserve(catalog_entries.size());
    for (const SliceCatalogEntry &entry : catalog_entries) {
      if (MatchSelectedBin(entry, config.fit_centrality_bins, config.fit_mt_bins)) {
        selected_entries.push_back(&entry);
      }
    }

    FitRunStatistics statistics;
    statistics.catalog_slices = catalog_entries.size();
    const std::size_t total_selected_slices = selected_entries.size();
    std::size_t processed_selected_slices = 0;
    ProgressReporter progress(logger, "fit", total_selected_slices, config.fit.progress);

    const FitModel model = override_model.value_or(config.fit.model);
    logger.Info("Starting fit stage with model=" + ToString(model) + ", coulombMode="
                + ToString(config.fit.options.coulomb_mode) + ", finiteSourceMode="
                + (config.fit.options.coulomb_mode == CoulombMode::kFiniteSource
                       ? ToString(config.fit.options.finite_source_mode)
                       : std::string(""))
                + ", inputCFPhiMapping=" + std::string(input_cf_uses_symmetric_phi_range ? "symmetric" : "raw")
                + ", fitPhiMapping=" + std::string(fit_uses_symmetric_phi_range ? "symmetric" : "raw") + ".");
    progress.Update(0);

    std::map<FitResultGroupKey, CoulombKernelTable> finite_source_kernels;
    std::vector<CoulombKernelCatalogEntry> kernel_catalog_entries;

    auto load_raw_histograms_if_needed = [&](const SliceCatalogEntry &entry,
                                             std::unique_ptr<TH3D> &h_se_raw,
                                             std::unique_ptr<TH3D> &h_me_raw,
                                             const bool throw_on_missing) -> bool {
      if (!config.fit.options.use_pml) {
        return true;
      }
      h_se_raw.reset(LoadStoredHistogram3D(*input_file, entry.se_object_path, entry.slice_id + "_se_raw"));
      h_me_raw.reset(LoadStoredHistogram3D(*input_file, entry.me_object_path, entry.slice_id + "_me_raw"));
      if (h_se_raw && h_me_raw) {
        return true;
      }
      if (throw_on_missing) {
        throw std::runtime_error("PML requested but raw SE/ME histograms are missing for finite-source seed slice "
                                 + entry.slice_id + ".");
      }
      return false;
    };

    if (config.fit.options.coulomb_mode == CoulombMode::kFiniteSource) {
      std::vector<FitResultGroupKey> group_order;
      std::map<FitResultGroupKey, const SliceCatalogEntry *> phi_integrated_by_group;
      for (const SliceCatalogEntry *entry : selected_entries) {
        const FitResultGroupKey key{entry->centrality_index, entry->mt_index, entry->qn_index};
        if (phi_integrated_by_group.find(key) == phi_integrated_by_group.end()
            && std::find(group_order.begin(), group_order.end(), key) == group_order.end()) {
          group_order.push_back(key);
        }
        if (entry->is_phi_integrated) {
          phi_integrated_by_group[key] = entry;
        }
      }

      logger.Info("Preparing finite-source Coulomb kernels for " + std::to_string(group_order.size())
                  + " cent/mT groups.");
      for (const FitResultGroupKey &key : group_order) {
        const auto seed_entry_iter = phi_integrated_by_group.find(key);
        if (seed_entry_iter == phi_integrated_by_group.end()) {
          throw std::runtime_error("Finite-source Coulomb requires a phi-integrated seed slice for selected group.");
        }
        const SliceCatalogEntry &seed_entry = *seed_entry_iter->second;
        std::unique_ptr<TH3D> h_cf(
            LoadStoredHistogram3D(*input_file, seed_entry.cf_object_path, seed_entry.slice_id + "_finite_seed_cf"));
        if (!h_cf) {
          throw std::runtime_error("Missing CF histogram for finite-source seed slice " + seed_entry.slice_id + ".");
        }
        std::unique_ptr<TH3D> h_se_raw;
        std::unique_ptr<TH3D> h_me_raw;
        load_raw_histograms_if_needed(seed_entry, h_se_raw, h_me_raw, true);

        LevyFitOptions seed_options = config.fit.options;
        seed_options.coulomb_mode = CoulombMode::kGamow;
        auto seed_fit = FitSingleSlice(h_cf.get(),
                                       h_se_raw.get(),
                                       h_me_raw.get(),
                                       seed_entry,
                                       model,
                                       seed_options,
                                       fit_uses_symmetric_phi_range,
                                       nullptr);
        if (!seed_fit.has_value()) {
          throw std::runtime_error("Gamow seed fit failed for finite-source group " + seed_entry.group_id + ".");
        }
        const double seed_radius_fm = ComputeEffectiveRadiusFromResult(seed_fit->result);
        double final_radius_fm = seed_radius_fm;
        CoulombKernelTable working_kernel = BuildFiniteSourceCoulombKernel(
            seed_entry, config.fit.options.finite_source_mode, seed_radius_fm, seed_radius_fm);

        if (config.fit.options.finite_source_mode == FiniteSourceMode::kIterative1D) {
          auto iterative_fit = FitSingleSlice(h_cf.get(),
                                             h_se_raw.get(),
                                             h_me_raw.get(),
                                             seed_entry,
                                             model,
                                             config.fit.options,
                                             fit_uses_symmetric_phi_range,
                                             &working_kernel);
          if (!iterative_fit.has_value()) {
            throw std::runtime_error("Finite-source iterative seed fit failed for group " + seed_entry.group_id + ".");
          }
          final_radius_fm = ComputeEffectiveRadiusFromResult(iterative_fit->result);
          working_kernel = BuildFiniteSourceCoulombKernel(
              seed_entry, config.fit.options.finite_source_mode, seed_radius_fm, final_radius_fm);
        }

        logger.Info("Finite-source kernel " + seed_entry.group_id + ": seed R_eff="
                    + FormatDouble(seed_radius_fm, 3) + " fm, final R_eff=" + FormatDouble(final_radius_fm, 3)
                    + " fm, k*=" + FormatDouble(working_kernel.catalog_entry.kstar_min_mev, 1) + "-"
                    + FormatDouble(working_kernel.catalog_entry.kstar_max_mev, 1) + " MeV/c.");
        kernel_catalog_entries.push_back(working_kernel.catalog_entry);
        finite_source_kernels.emplace(key, std::move(working_kernel));
      }
    }

    std::unique_ptr<TFile> shared_output_file;
    if (!config.fit.reopen_output_file_per_slice) {
      shared_output_file = OpenRootFile(fit_root_path, "UPDATE");
    }

    std::vector<LevyFitResult> fit_results;
    for (const SliceCatalogEntry *entry_ptr : selected_entries) {
      const SliceCatalogEntry &entry = *entry_ptr;
      ++statistics.selected_slices;
      logger.Debug("Fitting slice " + entry.slice_id);

      std::unique_ptr<TH3D> h_cf(LoadStoredHistogram3D(*input_file, entry.cf_object_path, entry.slice_id + "_cf_data"));
      if (!h_cf) {
        logger.Warn("Missing CF histogram for slice " + entry.slice_id);
        ++statistics.skipped_missing_objects;
        ++processed_selected_slices;
        progress.Update(processed_selected_slices);
        continue;
      }

      std::unique_ptr<TH3D> h_se_raw;
      std::unique_ptr<TH3D> h_me_raw;
      if (!load_raw_histograms_if_needed(entry, h_se_raw, h_me_raw, false)) {
        logger.Warn("PML requested but raw SE/ME histograms are missing for " + entry.slice_id);
        ++statistics.skipped_missing_raw_histograms;
        ++processed_selected_slices;
        progress.Update(processed_selected_slices);
        continue;
      }

      const CoulombKernelTable *coulomb_kernel = nullptr;
      if (config.fit.options.coulomb_mode == CoulombMode::kFiniteSource) {
        const FitResultGroupKey key{entry.centrality_index, entry.mt_index, entry.qn_index};
        const auto kernel_iter = finite_source_kernels.find(key);
        if (kernel_iter == finite_source_kernels.end()) {
          throw std::runtime_error("Missing finite-source Coulomb kernel for group " + entry.group_id + ".");
        }
        coulomb_kernel = &kernel_iter->second;
      }

      auto fit_output = FitSingleSlice(h_cf.get(),
                                       h_se_raw.get(),
                                       h_me_raw.get(),
                                       entry,
                                       model,
                                       config.fit.options,
                                       fit_uses_symmetric_phi_range,
                                       coulomb_kernel);
      if (fit_output.has_value()) {
        WriteSingleSliceFitArtifacts(h_cf.get(),
                                     entry,
                                     fit_output->fit_function.get(),
                                     fit_output->result,
                                     fit_root_path,
                                     shared_output_file.get(),
                                     coulomb_kernel);
        fit_results.push_back(fit_output->result);
        ++statistics.fitted_slices;
      }
      ++processed_selected_slices;
      progress.Update(processed_selected_slices);
    }

    if (shared_output_file) {
      WriteR2Graphs(*shared_output_file, fit_results);
      WriteFitCatalogTree(*shared_output_file, fit_results);
      WriteCoulombKernelCatalogTree(*shared_output_file, kernel_catalog_entries);
      shared_output_file->Close();
      shared_output_file.reset();
    } else {
      auto output_file = OpenRootFile(fit_root_path, "UPDATE");
      WriteR2Graphs(*output_file, fit_results);
      WriteFitCatalogTree(*output_file, fit_results);
      WriteCoulombKernelCatalogTree(*output_file, kernel_catalog_entries);
      output_file->Close();
    }
    WriteFitResultsSummaryTsv(fit_summary_path, fit_results);
    WriteFitReportRootFile(fit_report_root_path, fit_results, kernel_catalog_entries);

    progress.Finish();
    TH1::AddDirectory(old_add_directory);
    logger.Info("Completed fit stage: fitted " + std::to_string(statistics.fitted_slices) + " slices.");
    return statistics;
  }

  std::vector<SliceCatalogEntry> LoadSliceCatalog(const std::string &cf_root_path) {
    auto input_file = OpenRootFile(cf_root_path, "READ");
    return ReadSliceCatalogTree(*input_file);
  }

}  // namespace exp_femto_3d
