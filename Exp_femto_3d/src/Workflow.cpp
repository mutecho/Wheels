#include "exp_femto_3d/Workflow.h"

#include <algorithm>
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
#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  namespace {

    constexpr double kProjection1DWindow = 0.04;
    constexpr double kHbarC = 0.1973269804;
    constexpr double kPiPiLikeSignBohrRadiusFm = 387.5;
    constexpr double kFullR2MatrixTolerance = 1e-10;
    constexpr double kLevyArgumentTolerance = 1e-12;
    constexpr double kInvalidFullModelCFValue = 1e6;
    constexpr double kFitPenaltyValue = 1e30;

    struct Levy3DPMLContext {
      TH3D *h_se_raw = nullptr;
      TH3D *h_me_raw = nullptr;
      bool use_full_model = false;
      LevyFitOptions fit_options;
      double raw_same_to_mixed_integral_ratio = 1.0;
    };

    Levy3DPMLContext g_levy_3d_pml_context;

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

    std::string BuildGroupId(const RangeBin &centrality_bin, const RangeBin &mt_bin) {
      return "cent_" + FormatDirectoryValue(centrality_bin.min, 2) + "-"
             + FormatDirectoryValue(centrality_bin.max, 2) + "__mt_" + FormatDirectoryValue(mt_bin.min, 2) + "-"
             + FormatDirectoryValue(mt_bin.max, 2);
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

    SliceCatalogEntry MakeSliceCatalogEntry(const RangeBin &centrality_bin,
                                            const RangeBin &mt_bin,
                                            const int centrality_index,
                                            const int mt_index,
                                            const int phi_index,
                                            const double raw_phi_low,
                                            const double raw_phi_high,
                                            const double raw_phi_center,
                                            const double display_phi_low,
                                            const double display_phi_high,
                                            const double display_phi_center,
                                            const bool is_phi_integrated) {
      SliceCatalogEntry entry;
      entry.group_id = BuildGroupId(centrality_bin, mt_bin);
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
      entry.phi_index = phi_index;
      entry.cent_low = centrality_bin.min;
      entry.cent_high = centrality_bin.max;
      entry.mt_low = mt_bin.min;
      entry.mt_high = mt_bin.max;
      entry.raw_phi_low = raw_phi_low;
      entry.raw_phi_high = raw_phi_high;
      entry.raw_phi_center = raw_phi_center;
      entry.display_phi_low = display_phi_low;
      entry.display_phi_high = display_phi_high;
      entry.display_phi_center = display_phi_center;
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
      int phi_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      double raw_phi_low = 0.0;
      double raw_phi_high = 0.0;
      double raw_phi_center = 0.0;
      double display_phi_low = 0.0;
      double display_phi_high = 0.0;
      double display_phi_center = 0.0;
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
      tree->Branch("phi_index", &phi_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("raw_phi_low", &raw_phi_low);
      tree->Branch("raw_phi_high", &raw_phi_high);
      tree->Branch("raw_phi_center", &raw_phi_center);
      tree->Branch("display_phi_low", &display_phi_low);
      tree->Branch("display_phi_high", &display_phi_high);
      tree->Branch("display_phi_center", &display_phi_center);
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
        phi_index = entry.phi_index;
        cent_low = entry.cent_low;
        cent_high = entry.cent_high;
        mt_low = entry.mt_low;
        mt_high = entry.mt_high;
        raw_phi_low = entry.raw_phi_low;
        raw_phi_high = entry.raw_phi_high;
        raw_phi_center = entry.raw_phi_center;
        display_phi_low = entry.display_phi_low;
        display_phi_high = entry.display_phi_high;
        display_phi_center = entry.display_phi_center;
        is_phi_integrated = entry.is_phi_integrated ? 1 : 0;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
      output_file.cd();
    }

    std::vector<SliceCatalogEntry> ReadSliceCatalogTree(TFile &input_file) {
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
      TTreeReaderValue<int> phi_index(reader, "phi_index");
      TTreeReaderValue<double> cent_low(reader, "cent_low");
      TTreeReaderValue<double> cent_high(reader, "cent_high");
      TTreeReaderValue<double> mt_low(reader, "mt_low");
      TTreeReaderValue<double> mt_high(reader, "mt_high");
      TTreeReaderValue<double> raw_phi_low(reader, "raw_phi_low");
      TTreeReaderValue<double> raw_phi_high(reader, "raw_phi_high");
      TTreeReaderValue<double> raw_phi_center(reader, "raw_phi_center");
      TTreeReaderValue<double> display_phi_low(reader, "display_phi_low");
      TTreeReaderValue<double> display_phi_high(reader, "display_phi_high");
      TTreeReaderValue<double> display_phi_center(reader, "display_phi_center");
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
        entry.phi_index = *phi_index;
        entry.cent_low = *cent_low;
        entry.cent_high = *cent_high;
        entry.mt_low = *mt_low;
        entry.mt_high = *mt_high;
        entry.raw_phi_low = *raw_phi_low;
        entry.raw_phi_high = *raw_phi_high;
        entry.raw_phi_center = *raw_phi_center;
        entry.display_phi_low = *display_phi_low;
        entry.display_phi_high = *display_phi_high;
        entry.display_phi_center = *display_phi_center;
        entry.is_phi_integrated = (*is_phi_integrated != 0);
        entries.push_back(entry);
      }
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

    double ComputeBowlerSinyukovLikeSignPiPiValue(const double norm,
                                                  const double lambda,
                                                  const double levy_exponent,
                                                  const LevyFitOptions &fit_options,
                                                  const double q_out,
                                                  const double q_side,
                                                  const double q_long) {
      const double lambda_eff = fit_options.use_core_halo_lambda ? lambda : 1.0;
      const double coulomb_factor =
          fit_options.use_coulomb ? ComputeLikeSignPiPiGamowFactor(q_out, q_side, q_long) : 1.0;
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
      fit_options.use_coulomb = parameters[8] > 0.5;
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
      fit_options.use_coulomb = parameters[11] > 0.5;
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
      fit_function->SetParName(8, "UseCoulomb");
      fit_function->SetParName(9, "UseCoreHaloLambda");
      fit_function->SetParameters(1.0, 0.5, 25.0, 25.0, 25.0, 1.5, 0.0, 0.0, 0.0, 1.0);
      fit_function->SetParLimits(0, 0.5, 1.5);
      fit_function->SetParLimits(1, 0.0, 1.0);
      fit_function->SetParLimits(2, 0.01, 400.0);
      fit_function->SetParLimits(3, 0.01, 400.0);
      fit_function->SetParLimits(4, 0.01, 400.0);
      fit_function->SetParLimits(5, 0.5, 2.0);
      fit_function->SetParLimits(6, baseline_min, baseline_max);
      fit_function->FixParameter(7, fit_options.use_q2_baseline ? 1.0 : 0.0);
      fit_function->FixParameter(8, fit_options.use_coulomb ? 1.0 : 0.0);
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
      fit_function->SetParName(11, "UseCoulomb");
      fit_function->SetParName(12, "UseCoreHaloLambda");
      const double initial_parameters[13] = {1.0, 0.5, 25.0, 25.0, 25.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.0};
      fit_function->SetParameters(initial_parameters);
      fit_function->SetParLimits(0, 0.5, 1.5);
      fit_function->SetParLimits(1, 0.0, 1.0);
      fit_function->SetParLimits(2, 0.01, 400.0);
      fit_function->SetParLimits(3, 0.01, 400.0);
      fit_function->SetParLimits(4, 0.01, 400.0);
      fit_function->SetParLimits(8, 0.5, 2.0);
      fit_function->SetParLimits(9, baseline_min, baseline_max);
      fit_function->FixParameter(10, fit_options.use_q2_baseline ? 1.0 : 0.0);
      fit_function->FixParameter(11, fit_options.use_coulomb ? 1.0 : 0.0);
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
      g_levy_3d_pml_context.raw_same_to_mixed_integral_ratio = raw_to_normalized_scale;

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
      return std::string("Coulomb(#pi^{#pm}#pi^{#pm}): ") + (fit_result.uses_coulomb ? "on" : "off")
             + ", core-halo: " + (fit_result.uses_core_halo_lambda ? "on" : "off")
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

    void FillFitResultFromFunction(const SliceCatalogEntry &entry,
                                   TF3 *fit_function,
                                   const FitModel model,
                                   const LevyFitOptions &fit_options,
                                   const double fit_statistic,
                                   const int fit_ndf,
                                   const int fit_status,
                                   const double fit_edm,
                                   const int fit_minuit_istat,
                                   LevyFitResult &fit_result) {
      fit_result.fit_model = ToString(model);
      fit_result.slice_id = entry.slice_id;
      fit_result.group_id = entry.group_id;
      fit_result.slice_directory = BuildFitDirectory(entry.slice_id);
      fit_result.centrality_index = entry.centrality_index;
      fit_result.mt_index = entry.mt_index;
      fit_result.phi_index = entry.phi_index;
      fit_result.cent_low = entry.cent_low;
      fit_result.cent_high = entry.cent_high;
      fit_result.mt_low = entry.mt_low;
      fit_result.mt_high = entry.mt_high;
      fit_result.is_phi_integrated = entry.is_phi_integrated;
      fit_result.phi = entry.is_phi_integrated ? std::numeric_limits<double>::quiet_NaN() : entry.display_phi_center;
      fit_result.has_off_diagonal = model == FitModel::kFull;
      fit_result.uses_coulomb = fit_options.use_coulomb;
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

    bool FitAndWriteSingleSlice(TH3D *h_cf,
                                TH3D *h_se_raw,
                                TH3D *h_me_raw,
                                const SliceCatalogEntry &entry,
                                const FitModel model,
                                const LevyFitOptions &fit_options,
                                const std::string &fit_root_path,
                                LevyFitResult &fit_result,
                                TFile *shared_output_file) {
      TF3 *fit_function = model == FitModel::kFull
                              ? BuildFullLevyFitFunction(entry.slice_id + "_levy3d_full_fit", fit_options)
                              : BuildLevyFitFunction(entry.slice_id + "_levy3d_fit", fit_options);

      double fit_statistic = 0.0;
      double fit_edm = std::numeric_limits<double>::quiet_NaN();
      int fit_ndf = 0;
      int fit_status = -1;
      int fit_minuit_istat = -1;
      bool fit_succeeded = false;
      if (fit_options.use_pml) {
        fit_succeeded = RunPMLFit(fit_function,
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
        const auto fit_status_object = h_cf->Fit(fit_function, "RSMNQ0");
        fit_statistic = fit_function->GetChisquare();
        fit_ndf = fit_function->GetNDF();
        fit_status = static_cast<int>(fit_status_object);
        fit_succeeded = true;
      }

      FillFitResultFromFunction(entry,
                                fit_function,
                                model,
                                fit_options,
                                fit_statistic,
                                fit_ndf,
                                fit_status,
                                fit_edm,
                                fit_minuit_istat,
                                fit_result);
      if (!fit_succeeded) {
        delete fit_function;
        return false;
      }

      const bool old_batch_mode = gROOT->IsBatch();
      gROOT->SetBatch(kTRUE);

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
      delete fit_function;

      gROOT->SetBatch(old_batch_mode);
      return true;
    }

    void WriteFitResultsSummaryTsv(const std::string &path, const std::vector<LevyFitResult> &results) {
      std::ofstream output(path);
      output << "sliceId\tgroupId\tfitModel\tusesCoulomb\tusesCoreHaloLambda\tusesQ2Baseline"
                "\tusesPML\tcentLow\tcentHigh\tmTLow\tmTHigh\tphi\tisPhiIntegrated"
                "\tnorm\tnormErr\tlambda\tlambdaErr\tRout2\tRout2Err\tRside2\tRside2Err"
                "\tRlong2\tRlong2Err\tRoutside2\tRoutside2Err\tRoutlong2\tRoutlong2Err"
                "\tRsidelong2\tRsidelong2Err\talpha\talphaErr\tbaselineQ2\tbaselineQ2Err"
                "\tfitStatistic\tedm\tndf\tstatus\tminuitIstat\tcovarianceQuality\n";
      output << std::fixed << std::setprecision(6);
      for (const LevyFitResult &result : results) {
        output << result.slice_id << "\t" << result.group_id << "\t" << result.fit_model << "\t"
               << (result.uses_coulomb ? 1 : 0) << "\t" << (result.uses_core_halo_lambda ? 1 : 0) << "\t"
               << (result.uses_q2_baseline ? 1 : 0) << "\t" << (result.uses_pml ? 1 : 0) << "\t" << result.cent_low
               << "\t" << result.cent_high << "\t" << result.mt_low << "\t" << result.mt_high << "\t" << result.phi
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
      int phi_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      double phi = 0.0;
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
      int uses_core_halo_lambda = 0;
      int uses_q2_baseline = 0;
      int uses_pml = 0;
      int has_off_diagonal = 0;

      tree->Branch("slice_id", &slice_id);
      tree->Branch("group_id", &group_id);
      tree->Branch("fit_model", &fit_model);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("phi_index", &phi_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("phi", &phi);
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
        phi_index = result.phi_index;
        cent_low = result.cent_low;
        cent_high = result.cent_high;
        mt_low = result.mt_low;
        mt_high = result.mt_high;
        phi = result.phi;
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
        uses_core_halo_lambda = result.uses_core_halo_lambda ? 1 : 0;
        uses_q2_baseline = result.uses_q2_baseline ? 1 : 0;
        uses_pml = result.uses_pml ? 1 : 0;
        has_off_diagonal = result.has_off_diagonal ? 1 : 0;
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

        const bool uses_mapped_phi_range =
            std::any_of(group_results.begin(), group_results.end(), [](const LevyFitResult &result) {
              return result.phi < -1.0e-6;
            });
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

    BuildCfRunStatistics statistics;
    statistics.requested_groups = config.centrality_bins.size() * config.mt_bins.size();
    std::vector<SliceCatalogEntry> catalog_entries;
    // Count planned slices up front so skipped branches still advance the CLI progress bar.
    const std::size_t slices_per_group = static_cast<std::size_t>(n_phi_bins + 1);
    const std::size_t planned_slices = statistics.requested_groups * slices_per_group;
    std::size_t processed_slices = 0;
    ProgressReporter progress(logger, "build-cf", planned_slices);

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
        const std::string group_id = BuildGroupId(centrality_bin, mt_bin);
        logger.Debug("Building group " + group_id);

        h_se_origin->GetAxis(4)->SetRangeUser(centrality_bin.min, centrality_bin.max);
        h_se_origin->GetAxis(3)->SetRangeUser(mt_bin.min, mt_bin.max);
        h_se_origin->GetAxis(6)->SetRange(1, n_phi_bins);
        h_me_origin->GetAxis(4)->SetRangeUser(centrality_bin.min, centrality_bin.max);
        h_me_origin->GetAxis(3)->SetRangeUser(mt_bin.min, mt_bin.max);

        auto *h_me_raw = static_cast<TH3D *>(h_me_origin->Projection(0, 1, 2));
        h_me_raw->SetDirectory(nullptr);
        auto *h_me_norm = static_cast<TH3D *>(h_me_raw->Clone((group_id + "_ME_norm").c_str()));
        h_me_norm->SetDirectory(nullptr);
        const double int_me = IntegralVisibleRange(h_me_norm, true);
        if (int_me == 0.0) {
          logger.Warn("Zero mixed-event integral for " + group_id + "; skipping group.");
          ++statistics.skipped_zero_mixed_event_groups;
          processed_slices += slices_per_group;
          progress.Update(processed_slices);
          delete h_me_raw;
          delete h_me_norm;
          continue;
        }
        h_me_norm->Scale(1.0 / int_me);

        auto write_slice = [&](TH3D *h_se_raw, TH3D *h_se_norm, const SliceCatalogEntry &entry) {
          auto *h_cf = static_cast<TH3D *>(h_se_norm->Clone("CF3D"));
          h_cf->SetDirectory(nullptr);
          h_cf->Divide(h_me_norm);
          h_cf->GetXaxis()->SetTitle("q_{out} (GeV/c)");
          h_cf->GetYaxis()->SetTitle("q_{side} (GeV/c)");
          h_cf->GetZaxis()->SetTitle("q_{long} (GeV/c)");

          auto *h_me_write = static_cast<TH3D *>(h_me_raw->Clone("ME_raw3d"));
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
            Write1DProjections(h_me_norm, *directory, "ME_norm3d", "Normalized density", true);
          }
          output_file->cd();
          catalog_entries.push_back(entry);
          ++statistics.stored_slices;
          ++processed_slices;
          progress.Update(processed_slices);

          delete h_me_write;
          delete h_cf;
        };

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
          const double raw_phi_low = phi_axis->GetBinLowEdge(1);
          const double raw_phi_high = phi_axis->GetBinUpEdge(phi_axis->GetNbins());
          const auto entry = MakeSliceCatalogEntry(centrality_bin,
                                                   mt_bin,
                                                   static_cast<int>(centrality_index),
                                                   static_cast<int>(mt_index),
                                                   -1,
                                                   raw_phi_low,
                                                   raw_phi_high,
                                                   0.5 * (raw_phi_low + raw_phi_high),
                                                   raw_phi_low,
                                                   raw_phi_high,
                                                   std::numeric_limits<double>::quiet_NaN(),
                                                   true);
          write_slice(h_se_all_raw, h_se_all_norm, entry);
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

          const double raw_phi_low = phi_axis->GetBinLowEdge(phi_index);
          const double raw_phi_high = phi_axis->GetBinUpEdge(phi_index);
          const double raw_phi_center = phi_axis->GetBinCenter(phi_index);
          double display_phi_low = raw_phi_low;
          double display_phi_high = raw_phi_high;
          double display_phi_center = raw_phi_center;
          if (config.build.map_pair_phi_to_symmetric_range && raw_phi_center > TMath::Pi() / 2.0) {
            display_phi_low -= TMath::Pi();
            display_phi_high -= TMath::Pi();
            display_phi_center -= TMath::Pi();
          }

          const auto entry = MakeSliceCatalogEntry(centrality_bin,
                                                   mt_bin,
                                                   static_cast<int>(centrality_index),
                                                   static_cast<int>(mt_index),
                                                   phi_index - 1,
                                                   raw_phi_low,
                                                   raw_phi_high,
                                                   raw_phi_center,
                                                   display_phi_low,
                                                   display_phi_high,
                                                   display_phi_center,
                                                   false);
          write_slice(h_se_raw, h_se_norm, entry);
          delete h_se_raw;
          delete h_se_norm;
        }

        h_se_origin->GetAxis(6)->SetRange(1, n_phi_bins);

        delete h_me_raw;
        delete h_me_norm;
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
    const std::string cf_root_path = input_cf_root_path.has_value()
                                         ? ResolvePath(config.output.output_directory, *input_cf_root_path)
                                         : ResolvePath(config.output.output_directory, config.output.cf_root_name);
    const std::string fit_root_path = ResolvePath(config.output.output_directory, config.output.fit_root_name);
    const std::string fit_summary_path = ResolvePath(config.output.output_directory, config.output.fit_summary_name);

    CreateOrResetRootFile(fit_root_path);
    const bool old_add_directory = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    auto input_file = OpenRootFile(cf_root_path, "READ");
    const std::vector<SliceCatalogEntry> catalog_entries = ReadSliceCatalogTree(*input_file);

    FitRunStatistics statistics;
    statistics.catalog_slices = catalog_entries.size();
    const std::size_t total_selected_slices = static_cast<std::size_t>(std::count_if(
        catalog_entries.begin(), catalog_entries.end(), [&](const SliceCatalogEntry &entry) {
          return MatchSelectedBin(entry, config.fit_centrality_bins, config.fit_mt_bins);
        }));
    std::size_t processed_selected_slices = 0;
    ProgressReporter progress(logger, "fit", total_selected_slices);

    std::unique_ptr<TFile> shared_output_file;
    if (!config.fit.reopen_output_file_per_slice) {
      shared_output_file = OpenRootFile(fit_root_path, "UPDATE");
    }

    const FitModel model = override_model.value_or(config.fit.model);
    logger.Info("Starting fit stage with model=" + ToString(model) + ".");
    progress.Update(0);

    std::vector<LevyFitResult> fit_results;
    for (const SliceCatalogEntry &entry : catalog_entries) {
      if (!MatchSelectedBin(entry, config.fit_centrality_bins, config.fit_mt_bins)) {
        continue;
      }
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
      if (config.fit.options.use_pml) {
        h_se_raw.reset(LoadStoredHistogram3D(*input_file, entry.se_object_path, entry.slice_id + "_se_raw"));
        h_me_raw.reset(LoadStoredHistogram3D(*input_file, entry.me_object_path, entry.slice_id + "_me_raw"));
        if (!h_se_raw || !h_me_raw) {
          logger.Warn("PML requested but raw SE/ME histograms are missing for " + entry.slice_id);
          ++statistics.skipped_missing_raw_histograms;
          ++processed_selected_slices;
          progress.Update(processed_selected_slices);
          continue;
        }
      }

      LevyFitResult fit_result;
      if (FitAndWriteSingleSlice(h_cf.get(),
                                 h_se_raw.get(),
                                 h_me_raw.get(),
                                 entry,
                                 model,
                                 config.fit.options,
                                 fit_root_path,
                                 fit_result,
                                 shared_output_file.get())) {
        fit_results.push_back(fit_result);
        ++statistics.fitted_slices;
      }
      ++processed_selected_slices;
      progress.Update(processed_selected_slices);
    }

    if (shared_output_file) {
      WriteR2Graphs(*shared_output_file, fit_results);
      WriteFitCatalogTree(*shared_output_file, fit_results);
      shared_output_file->Close();
      shared_output_file.reset();
    } else {
      auto output_file = OpenRootFile(fit_root_path, "UPDATE");
      WriteR2Graphs(*output_file, fit_results);
      WriteFitCatalogTree(*output_file, fit_results);
      output_file->Close();
    }
    WriteFitResultsSummaryTsv(fit_summary_path, fit_results);

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
