#include "exp_femto_1d/Workflow.h"

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
#include <utility>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "exp_femto_1d/CatsModel.h"
#include "exp_femto_1d/Config.h"

namespace exp_femto_1d {

  namespace {

    struct RegionDefinition {
      RegionKind kind = RegionKind::kMinBias;
      std::string name;
      double low_1 = 0.0;
      double high_1 = 0.0;
      double low_2 = 0.0;
      double high_2 = 0.0;
      bool has_second_interval = false;
      int region_index = 0;
    };

    struct SliceHistograms {
      std::unique_ptr<TH1D> se_raw;
      std::unique_ptr<TH1D> me_raw;
      std::unique_ptr<TH1D> cf;
    };

    std::string FormatDouble(const double value, const int precision = 2) {
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(precision) << value;
      return stream.str();
    }

    // Keep directory names stable by normalizing near-zero values.
    std::string FormatDirectoryValue(const double value, const int precision = 2) {
      const double stable_value = std::abs(value) < 5.0e-7 ? 0.0 : value;
      return FormatDouble(stable_value, precision);
    }

    std::string BuildSparseObjectPath(const InputConfig &input, const std::string &subtask) {
      return input.task_name + "/" + subtask + "/" + input.sparse_object_name;
    }

    std::string BuildGroupId(const RangeBin &centrality_bin, const RangeBin &mt_bin) {
      return "cent_" + FormatDirectoryValue(centrality_bin.min, 2) + "-" + FormatDirectoryValue(centrality_bin.max, 2)
             + "__mt_" + FormatDirectoryValue(mt_bin.min, 3) + "-" + FormatDirectoryValue(mt_bin.max, 3);
    }

    std::string BuildSliceId(const std::string &group_id, const RegionDefinition &region) {
      return group_id + "__region_" + region.name;
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

    // Build nested ROOT directory trees without leaking the path semantics into callers.
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

    std::array<RegionDefinition, 3> BuildRegionDefinitions() {
      return {
          RegionDefinition{RegionKind::kMinBias, "MinBias", 0.0, kPi, 0.0, 0.0, false, 0},
          RegionDefinition{RegionKind::kInPlane, "InPlane", 0.0, kPi / 4.0, 3.0 * kPi / 4.0, kPi, true, 1},
          RegionDefinition{RegionKind::kOutOfPlane, "OutOfPlane", kPi / 4.0, 3.0 * kPi / 4.0, 0.0, 0.0, false, 2},
      };
    }

    int RegionSummaryBin(const RegionKind region_kind) {
      switch (region_kind) {
        case RegionKind::kMinBias:
          return 1;
        case RegionKind::kInPlane:
          return 2;
        case RegionKind::kOutOfPlane:
          return 3;
      }
      return 1;
    }

    void ResetAxisSelection(TAxis &axis) {
      axis.SetRange(0, 0);
    }

    std::pair<int, int> FindIntegralBinRange(const TH1D &histogram, const double low, const double high) {
      const double epsilon = 1.0e-9;
      const int first_bin = histogram.GetXaxis()->FindFixBin(low + epsilon);
      const int last_bin = histogram.GetXaxis()->FindFixBin(high - epsilon);
      return {first_bin, last_bin};
    }

    double IntegralInRange(const TH1D &histogram, const double low, const double high) {
      const auto [first_bin, last_bin] = FindIntegralBinRange(histogram, low, high);
      return histogram.Integral(first_bin, last_bin);
    }

    // Copy only the requested k* window into a standalone histogram for the public output contract.
    std::unique_ptr<TH1D> BuildWindowedHistogram(const TH1D &source,
                                                 const std::string &name,
                                                 const std::string &title,
                                                 const double low,
                                                 const double high) {
      const double epsilon = 1.0e-9;
      const int first_bin = source.GetXaxis()->FindFixBin(low + epsilon);
      const int last_bin = source.GetXaxis()->FindFixBin(high - epsilon);
      if (last_bin < first_bin) {
        throw std::runtime_error("Invalid k* copy range.");
      }

      std::vector<double> edges;
      edges.reserve(static_cast<std::size_t>(last_bin - first_bin + 2));
      for (int bin = first_bin; bin <= last_bin; ++bin) {
        edges.push_back(source.GetXaxis()->GetBinLowEdge(bin));
      }
      edges.push_back(source.GetXaxis()->GetBinUpEdge(last_bin));

      auto histogram = std::make_unique<TH1D>(name.c_str(), title.c_str(), last_bin - first_bin + 1, edges.data());
      histogram->SetDirectory(nullptr);
      histogram->Sumw2();
      histogram->GetXaxis()->SetTitle(source.GetXaxis()->GetTitle());
      histogram->GetYaxis()->SetTitle(source.GetYaxis()->GetTitle());
      for (int source_bin = first_bin; source_bin <= last_bin; ++source_bin) {
        const int target_bin = source_bin - first_bin + 1;
        histogram->SetBinContent(target_bin, source.GetBinContent(source_bin));
        histogram->SetBinError(target_bin, source.GetBinError(source_bin));
      }
      return histogram;
    }

    // Project one EP region from the 4D sparse after selecting the requested cent/mT slice.
    std::unique_ptr<TH1D> BuildProjection(THnSparseF &sparse, const RangeBin &cent_bin, const RangeBin &mt_bin, const RegionDefinition &region) {
      auto *axis_mt = sparse.GetAxis(1);
      auto *axis_cent = sparse.GetAxis(2);
      auto *axis_ep = sparse.GetAxis(3);

      axis_cent->SetRangeUser(cent_bin.min, cent_bin.max);
      axis_mt->SetRangeUser(mt_bin.min, mt_bin.max);
      axis_ep->SetRangeUser(region.low_1, region.high_1);
      auto *projection = static_cast<TH1D *>(sparse.Projection(0));
      if (projection == nullptr) {
        ResetAxisSelection(*axis_ep);
        ResetAxisSelection(*axis_mt);
        ResetAxisSelection(*axis_cent);
        throw std::runtime_error("Failed to project THnSparseF to TH1D.");
      }

      std::unique_ptr<TH1D> merged_projection(projection);
      if (region.has_second_interval) {
        axis_ep->SetRangeUser(region.low_2, region.high_2);
        std::unique_ptr<TH1D> second_projection(static_cast<TH1D *>(sparse.Projection(0)));
        if (second_projection != nullptr) {
          merged_projection->Add(second_projection.get());
        }
      }

      ResetAxisSelection(*axis_ep);
      ResetAxisSelection(*axis_mt);
      ResetAxisSelection(*axis_cent);

      merged_projection->SetDirectory(nullptr);
      merged_projection->Sumw2();
      merged_projection->GetXaxis()->SetTitle("k* (GeV/c)");
      merged_projection->GetYaxis()->SetTitle("Counts");
      return merged_projection;
    }

    // Build one output slice from raw SE and ME projections using the requested normalization range.
    std::optional<SliceHistograms> BuildCf1D(const TH1D &se_projection,
                                             const TH1D &me_projection,
                                             const BuildCfConfig &build_config) {
      const double se_norm = IntegralInRange(se_projection, build_config.norm_low, build_config.norm_high);
      const double me_norm = IntegralInRange(me_projection, build_config.norm_low, build_config.norm_high);
      if (se_norm <= 0.0 || me_norm <= 0.0 || !std::isfinite(se_norm) || !std::isfinite(me_norm)) {
        return std::nullopt;
      }

      auto se_window = BuildWindowedHistogram(
          se_projection, "SE_raw1d", "SE_raw1d; k* (GeV/c); Counts", build_config.kstar_min, build_config.kstar_max);
      auto me_window = BuildWindowedHistogram(
          me_projection, "ME_raw1d", "ME_raw1d; k* (GeV/c); Counts", build_config.kstar_min, build_config.kstar_max);

      const double normalization_factor = me_norm / se_norm;
      auto me_normalized = std::unique_ptr<TH1D>(static_cast<TH1D *>(me_window->Clone("ME_normalized")));
      me_normalized->SetDirectory(nullptr);
      me_normalized->Scale(1.0 / normalization_factor);

      auto cf_histogram = std::unique_ptr<TH1D>(static_cast<TH1D *>(se_window->Clone("CF1D")));
      cf_histogram->SetDirectory(nullptr);
      cf_histogram->SetTitle("CF1D; k* (GeV/c); C(k*)");
      cf_histogram->Divide(me_normalized.get());

      me_window->SetName("ME_raw1d");
      return SliceHistograms{std::move(se_window), std::move(me_window), std::move(cf_histogram)};
    }

    void WriteSliceHistograms(TFile &output_file, const SliceCatalogEntry &entry, const SliceHistograms &histograms) {
      auto *slice_directory = GetOrCreateDirectoryPath(output_file, entry.slice_directory);
      slice_directory->cd();
      histograms.se_raw->Write("SE_raw1d", TObject::kOverwrite);
      histograms.me_raw->Write("ME_raw1d", TObject::kOverwrite);
      histograms.cf->Write("CF1D", TObject::kOverwrite);
    }

    // Persist the catalog so fit never has to recover metadata from object names.
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
      int centrality_index = -1;
      int mt_index = -1;
      int region_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      std::string region_name;
      int region_kind = 0;
      double ep_low_1 = 0.0;
      double ep_high_1 = 0.0;
      double ep_low_2 = 0.0;
      double ep_high_2 = 0.0;
      int has_second_interval = 0;
      double norm_low = 0.0;
      double norm_high = 0.0;
      double kstar_min = 0.0;
      double kstar_max = 0.0;

      tree->Branch("slice_id", &slice_id);
      tree->Branch("group_id", &group_id);
      tree->Branch("slice_directory", &slice_directory);
      tree->Branch("se_object_path", &se_object_path);
      tree->Branch("me_object_path", &me_object_path);
      tree->Branch("cf_object_path", &cf_object_path);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("region_index", &region_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("region_name", &region_name);
      tree->Branch("region_kind", &region_kind);
      tree->Branch("ep_low_1", &ep_low_1);
      tree->Branch("ep_high_1", &ep_high_1);
      tree->Branch("ep_low_2", &ep_low_2);
      tree->Branch("ep_high_2", &ep_high_2);
      tree->Branch("has_second_interval", &has_second_interval);
      tree->Branch("norm_low", &norm_low);
      tree->Branch("norm_high", &norm_high);
      tree->Branch("kstar_min", &kstar_min);
      tree->Branch("kstar_max", &kstar_max);

      for (const SliceCatalogEntry &entry : entries) {
        slice_id = entry.slice_id;
        group_id = entry.group_id;
        slice_directory = entry.slice_directory;
        se_object_path = entry.se_object_path;
        me_object_path = entry.me_object_path;
        cf_object_path = entry.cf_object_path;
        centrality_index = entry.centrality_index;
        mt_index = entry.mt_index;
        region_index = entry.region_index;
        cent_low = entry.cent_low;
        cent_high = entry.cent_high;
        mt_low = entry.mt_low;
        mt_high = entry.mt_high;
        region_name = entry.region_name;
        region_kind = static_cast<int>(entry.region_kind);
        ep_low_1 = entry.ep_low_1;
        ep_high_1 = entry.ep_high_1;
        ep_low_2 = entry.ep_low_2;
        ep_high_2 = entry.ep_high_2;
        has_second_interval = entry.has_second_interval ? 1 : 0;
        norm_low = entry.norm_low;
        norm_high = entry.norm_high;
        kstar_min = entry.kstar_min;
        kstar_max = entry.kstar_max;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
    }

    std::vector<SliceCatalogEntry> ReadSliceCatalogTree(TTree &tree) {
      std::vector<SliceCatalogEntry> entries;
      TTreeReader reader(&tree);
      TTreeReaderValue<std::string> slice_id(reader, "slice_id");
      TTreeReaderValue<std::string> group_id(reader, "group_id");
      TTreeReaderValue<std::string> slice_directory(reader, "slice_directory");
      TTreeReaderValue<std::string> se_object_path(reader, "se_object_path");
      TTreeReaderValue<std::string> me_object_path(reader, "me_object_path");
      TTreeReaderValue<std::string> cf_object_path(reader, "cf_object_path");
      TTreeReaderValue<int> centrality_index(reader, "centrality_index");
      TTreeReaderValue<int> mt_index(reader, "mt_index");
      TTreeReaderValue<int> region_index(reader, "region_index");
      TTreeReaderValue<double> cent_low(reader, "cent_low");
      TTreeReaderValue<double> cent_high(reader, "cent_high");
      TTreeReaderValue<double> mt_low(reader, "mt_low");
      TTreeReaderValue<double> mt_high(reader, "mt_high");
      TTreeReaderValue<std::string> region_name(reader, "region_name");
      TTreeReaderValue<int> region_kind(reader, "region_kind");
      TTreeReaderValue<double> ep_low_1(reader, "ep_low_1");
      TTreeReaderValue<double> ep_high_1(reader, "ep_high_1");
      TTreeReaderValue<double> ep_low_2(reader, "ep_low_2");
      TTreeReaderValue<double> ep_high_2(reader, "ep_high_2");
      TTreeReaderValue<int> has_second_interval(reader, "has_second_interval");
      TTreeReaderValue<double> norm_low(reader, "norm_low");
      TTreeReaderValue<double> norm_high(reader, "norm_high");
      TTreeReaderValue<double> kstar_min(reader, "kstar_min");
      TTreeReaderValue<double> kstar_max(reader, "kstar_max");

      while (reader.Next()) {
        SliceCatalogEntry entry;
        entry.slice_id = *slice_id;
        entry.group_id = *group_id;
        entry.slice_directory = *slice_directory;
        entry.se_object_path = *se_object_path;
        entry.me_object_path = *me_object_path;
        entry.cf_object_path = *cf_object_path;
        entry.centrality_index = *centrality_index;
        entry.mt_index = *mt_index;
        entry.region_index = *region_index;
        entry.cent_low = *cent_low;
        entry.cent_high = *cent_high;
        entry.mt_low = *mt_low;
        entry.mt_high = *mt_high;
        entry.region_name = *region_name;
        entry.region_kind = static_cast<RegionKind>(*region_kind);
        entry.ep_low_1 = *ep_low_1;
        entry.ep_high_1 = *ep_high_1;
        entry.ep_low_2 = *ep_low_2;
        entry.ep_high_2 = *ep_high_2;
        entry.has_second_interval = (*has_second_interval != 0);
        entry.norm_low = *norm_low;
        entry.norm_high = *norm_high;
        entry.kstar_min = *kstar_min;
        entry.kstar_max = *kstar_max;
        entries.push_back(std::move(entry));
      }
      return entries;
    }

    std::unique_ptr<TH1D> LoadStoredHistogram1D(TFile &file, const std::string &path) {
      auto *histogram = dynamic_cast<TH1D *>(file.Get(path.c_str()));
      if (histogram == nullptr) {
        return nullptr;
      }
      auto clone = std::unique_ptr<TH1D>(static_cast<TH1D *>(histogram->Clone()));
      clone->SetDirectory(nullptr);
      return clone;
    }

    bool MatchSelectedBin(const RangeBin &bin, const std::vector<RangeBin> &selection) {
      return std::any_of(selection.begin(), selection.end(), [&](const RangeBin &candidate) {
        return MatchesRangeBin(bin, candidate);
      });
    }

    bool MatchSelectedSlice(const SliceCatalogEntry &entry,
                            const std::vector<RangeBin> &centrality_selection,
                            const std::vector<RangeBin> &mt_selection) {
      const RangeBin cent_bin{entry.cent_low, entry.cent_high, ""};
      const RangeBin mt_bin{entry.mt_low, entry.mt_high, ""};
      return MatchSelectedBin(cent_bin, centrality_selection) && MatchSelectedBin(mt_bin, mt_selection);
    }

    double EvaluateBaselinePolynomial(const double kstar_gev, const PiPiFitResult &result) {
      return result.baseline_p0
             * (1.0 + result.baseline_p1 * kstar_gev + result.baseline_p2 * kstar_gev * kstar_gev
                + result.baseline_p3 * std::pow(kstar_gev, 3) + result.baseline_p4 * std::pow(kstar_gev, 4));
    }

    // Sample the fitted baseline onto the same binning as the data for portable ROOT output.
    std::unique_ptr<TH1D> BuildBaselineHistogram(const TH1D &data_histogram, const PiPiFitResult &result) {
      auto baseline = std::unique_ptr<TH1D>(static_cast<TH1D *>(data_histogram.Clone("Baseline")));
      baseline->SetDirectory(nullptr);
      baseline->SetTitle("Baseline; k* (GeV/c); Baseline");
      for (int bin = 1; bin <= baseline->GetNbinsX(); ++bin) {
        const double kstar_gev = baseline->GetXaxis()->GetBinCenter(bin);
        baseline->SetBinContent(bin, EvaluateBaselinePolynomial(kstar_gev, result));
        baseline->SetBinError(bin, 0.0);
      }
      return baseline;
    }

    // Sample the pure femtoscopic CATS signal using the recovered source size.
    std::unique_ptr<TH1D> BuildPureFemtoHistogram(const TH1D &data_histogram,
                                                  PiPiCatsModel &model,
                                                  const PiPiFitResult &result) {
      auto pure_femto = std::unique_ptr<TH1D>(static_cast<TH1D *>(data_histogram.Clone("PureFemto")));
      pure_femto->SetDirectory(nullptr);
      pure_femto->SetTitle("PureFemto; k* (GeV/c); CATS");
      for (int bin = 1; bin <= pure_femto->GetNbinsX(); ++bin) {
        const double kstar_mev = 1000.0 * pure_femto->GetXaxis()->GetBinCenter(bin);
        pure_femto->SetBinContent(bin, model.Evaluate(kstar_mev, result.source_size));
        pure_femto->SetBinError(bin, 0.0);
      }
      return pure_femto;
    }

    // Build a simple diagnostic canvas that exposes data, the full fit, and its factorized pieces.
    std::unique_ptr<TCanvas> BuildFitCanvas(TH1D &data_histogram,
                                            TF1 &fit_function,
                                            TH1D &baseline_histogram,
                                            TH1D &pure_femto_histogram,
                                            const PiPiFitResult &result) {
      auto canvas = std::make_unique<TCanvas>("FitCanvas", "FitCanvas", 900, 700);
      data_histogram.SetMarkerStyle(20);
      data_histogram.SetMarkerSize(0.9);
      data_histogram.SetLineColor(kBlack);
      baseline_histogram.SetLineColor(kBlue + 1);
      baseline_histogram.SetLineWidth(2);
      pure_femto_histogram.SetLineColor(kGreen + 2);
      pure_femto_histogram.SetLineWidth(2);
      fit_function.SetLineColor(kRed + 1);
      fit_function.SetLineWidth(2);

      data_histogram.Draw("E1");
      fit_function.Draw("SAME");
      baseline_histogram.Draw("HIST SAME");
      pure_femto_histogram.Draw("HIST SAME");

      auto legend = std::make_unique<TLegend>(0.58, 0.68, 0.88, 0.88);
      legend->SetBorderSize(0);
      legend->AddEntry(&data_histogram, "DataCF", "lep");
      legend->AddEntry(&fit_function, "FitFunction", "l");
      legend->AddEntry(&baseline_histogram, "Baseline", "l");
      legend->AddEntry(&pure_femto_histogram, "PureFemto", "l");
      legend->AddEntry(static_cast<TObject *>(nullptr), ("status=" + std::to_string(result.status)).c_str(), "");
      legend->Draw();
      canvas->Update();
      return canvas;
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

    void WriteFitResultsSummaryTsv(const std::string &path, const std::vector<PiPiFitResult> &results) {
      std::ofstream output(path);
      output << "sliceId\tgroupId\tcentLow\tcentHigh\tmTLow\tmTHigh\tregionName\tfitKstarMax\tusesCoulomb"
                "\tp0\tp0Err\tp1\tp1Err\tp2\tp2Err\tp3\tp3Err\tp4\tp4Err\tsourceSize\tsourceSizeErr"
                "\tchi2\tndf\tfitStatistic\tedm\tstatus\tminuitIstat\tcovarianceQuality\n";
      output << std::fixed << std::setprecision(6);
      for (const PiPiFitResult &result : results) {
        output << result.slice_id << "\t" << result.group_id << "\t" << result.cent_low << "\t" << result.cent_high
               << "\t" << result.mt_low << "\t" << result.mt_high << "\t" << result.region_name << "\t"
               << result.fit_kstar_max << "\t" << (result.uses_coulomb ? 1 : 0) << "\t" << result.baseline_p0 << "\t"
               << result.baseline_p0_err << "\t" << result.baseline_p1 << "\t" << result.baseline_p1_err << "\t"
               << result.baseline_p2 << "\t" << result.baseline_p2_err << "\t" << result.baseline_p3 << "\t"
               << result.baseline_p3_err << "\t" << result.baseline_p4 << "\t" << result.baseline_p4_err << "\t"
               << result.source_size << "\t" << result.source_size_err << "\t" << result.chi2 << "\t" << result.ndf
               << "\t" << result.fit_statistic << "\t" << result.edm << "\t" << result.status << "\t"
               << result.minuit_istat << "\t" << CovarianceQualityToken(result.minuit_istat) << "\n";
      }
    }

    void WriteFitCatalogTree(TFile &output_file, const std::vector<PiPiFitResult> &results) {
      auto *meta_directory = GetOrCreateDirectoryPath(output_file, "meta");
      meta_directory->cd();

      auto tree = std::make_unique<TTree>("FitCatalog", "FitCatalog");
      std::string slice_id;
      std::string group_id;
      std::string slice_directory;
      int centrality_index = -1;
      int mt_index = -1;
      int region_index = -1;
      double cent_low = 0.0;
      double cent_high = 0.0;
      double mt_low = 0.0;
      double mt_high = 0.0;
      std::string region_name;
      double fit_kstar_max = 0.0;
      int uses_coulomb = 0;
      double baseline_p0 = 0.0;
      double baseline_p0_err = 0.0;
      double baseline_p1 = 0.0;
      double baseline_p1_err = 0.0;
      double baseline_p2 = 0.0;
      double baseline_p2_err = 0.0;
      double baseline_p3 = 0.0;
      double baseline_p3_err = 0.0;
      double baseline_p4 = 0.0;
      double baseline_p4_err = 0.0;
      double source_size = 0.0;
      double source_size_err = 0.0;
      double chi2 = 0.0;
      int ndf = 0;
      double fit_statistic = 0.0;
      double edm = 0.0;
      int status = 0;
      int minuit_istat = 0;

      tree->Branch("slice_id", &slice_id);
      tree->Branch("group_id", &group_id);
      tree->Branch("slice_directory", &slice_directory);
      tree->Branch("centrality_index", &centrality_index);
      tree->Branch("mt_index", &mt_index);
      tree->Branch("region_index", &region_index);
      tree->Branch("cent_low", &cent_low);
      tree->Branch("cent_high", &cent_high);
      tree->Branch("mt_low", &mt_low);
      tree->Branch("mt_high", &mt_high);
      tree->Branch("region_name", &region_name);
      tree->Branch("fit_kstar_max", &fit_kstar_max);
      tree->Branch("uses_coulomb", &uses_coulomb);
      tree->Branch("baseline_p0", &baseline_p0);
      tree->Branch("baseline_p0_err", &baseline_p0_err);
      tree->Branch("baseline_p1", &baseline_p1);
      tree->Branch("baseline_p1_err", &baseline_p1_err);
      tree->Branch("baseline_p2", &baseline_p2);
      tree->Branch("baseline_p2_err", &baseline_p2_err);
      tree->Branch("baseline_p3", &baseline_p3);
      tree->Branch("baseline_p3_err", &baseline_p3_err);
      tree->Branch("baseline_p4", &baseline_p4);
      tree->Branch("baseline_p4_err", &baseline_p4_err);
      tree->Branch("source_size", &source_size);
      tree->Branch("source_size_err", &source_size_err);
      tree->Branch("chi2", &chi2);
      tree->Branch("ndf", &ndf);
      tree->Branch("fit_statistic", &fit_statistic);
      tree->Branch("edm", &edm);
      tree->Branch("status", &status);
      tree->Branch("minuit_istat", &minuit_istat);

      for (const PiPiFitResult &result : results) {
        slice_id = result.slice_id;
        group_id = result.group_id;
        slice_directory = result.slice_directory;
        centrality_index = result.centrality_index;
        mt_index = result.mt_index;
        region_index = result.region_index;
        cent_low = result.cent_low;
        cent_high = result.cent_high;
        mt_low = result.mt_low;
        mt_high = result.mt_high;
        region_name = result.region_name;
        fit_kstar_max = result.fit_kstar_max;
        uses_coulomb = result.uses_coulomb ? 1 : 0;
        baseline_p0 = result.baseline_p0;
        baseline_p0_err = result.baseline_p0_err;
        baseline_p1 = result.baseline_p1;
        baseline_p1_err = result.baseline_p1_err;
        baseline_p2 = result.baseline_p2;
        baseline_p2_err = result.baseline_p2_err;
        baseline_p3 = result.baseline_p3;
        baseline_p3_err = result.baseline_p3_err;
        baseline_p4 = result.baseline_p4;
        baseline_p4_err = result.baseline_p4_err;
        source_size = result.source_size;
        source_size_err = result.source_size_err;
        chi2 = result.chi2;
        ndf = result.ndf;
        fit_statistic = result.fit_statistic;
        edm = result.edm;
        status = result.status;
        minuit_istat = result.minuit_istat;
        tree->Fill();
      }

      tree->Write("", TObject::kOverwrite);
    }

    // Collapse the three fixed regions into per-group summary histograms.
    void WriteRegionSummaryHistograms(TFile &output_file, const std::vector<PiPiFitResult> &results) {
      std::map<std::string, std::vector<PiPiFitResult>> grouped_results;
      for (const PiPiFitResult &result : results) {
        grouped_results[result.group_id].push_back(result);
      }

      for (auto &[group_id, group_results] : grouped_results) {
        auto *summary_directory = GetOrCreateDirectoryPath(output_file, "summary/by_region/" + group_id);
        summary_directory->cd();

        TH1D source_size("SourceSize", "SourceSize; Region; Source size (fm)", 3, 0.5, 3.5);
        TH1D norm("Norm", "Norm; Region; p0", 3, 0.5, 3.5);
        TH1D p1("P1", "P1; Region; p1", 3, 0.5, 3.5);
        TH1D p2("P2", "P2; Region; p2", 3, 0.5, 3.5);
        for (TH1D *histogram : std::array<TH1D *, 4>{&source_size, &norm, &p1, &p2}) {
          histogram->GetXaxis()->SetBinLabel(1, "MinBias");
          histogram->GetXaxis()->SetBinLabel(2, "InPlane");
          histogram->GetXaxis()->SetBinLabel(3, "OutOfPlane");
        }

        for (const PiPiFitResult &result : group_results) {
          const int summary_bin = RegionSummaryBin(static_cast<RegionKind>(result.region_index));
          source_size.SetBinContent(summary_bin, result.source_size);
          source_size.SetBinError(summary_bin, result.source_size_err);
          norm.SetBinContent(summary_bin, result.baseline_p0);
          norm.SetBinError(summary_bin, result.baseline_p0_err);
          p1.SetBinContent(summary_bin, result.baseline_p1);
          p1.SetBinError(summary_bin, result.baseline_p1_err);
          p2.SetBinContent(summary_bin, result.baseline_p2);
          p2.SetBinError(summary_bin, result.baseline_p2_err);
        }

        source_size.Write("SourceSize", TObject::kOverwrite);
        norm.Write("Norm", TObject::kOverwrite);
        p1.Write("P1", TObject::kOverwrite);
        p2.Write("P2", TObject::kOverwrite);
      }
    }

    PiPiFitResult BuildFitResult(const SliceCatalogEntry &entry,
                                 const FitConfig &fit_config,
                                 TF1 &fit_function,
                                 const int fit_status,
                                 const int covariance_status,
                                 const double edm) {
      PiPiFitResult result;
      result.slice_id = entry.slice_id;
      result.group_id = entry.group_id;
      result.slice_directory = BuildFitDirectory(entry.slice_id);
      result.centrality_index = entry.centrality_index;
      result.mt_index = entry.mt_index;
      result.region_index = entry.region_index;
      result.cent_low = entry.cent_low;
      result.cent_high = entry.cent_high;
      result.mt_low = entry.mt_low;
      result.mt_high = entry.mt_high;
      result.region_name = entry.region_name;
      result.fit_kstar_max = fit_config.fit_kstar_max;
      result.uses_coulomb = fit_config.use_coulomb;
      result.baseline_p0 = fit_function.GetParameter(0);
      result.baseline_p0_err = fit_function.GetParError(0);
      result.baseline_p1 = fit_function.GetParameter(1);
      result.baseline_p1_err = fit_function.GetParError(1);
      result.baseline_p2 = fit_function.GetParameter(2);
      result.baseline_p2_err = fit_function.GetParError(2);
      result.baseline_p3 = fit_function.GetParameter(3);
      result.baseline_p3_err = fit_function.GetParError(3);
      result.baseline_p4 = fit_function.GetParameter(4);
      result.baseline_p4_err = fit_function.GetParError(4);
      result.source_size = fit_function.GetParameter(5);
      result.source_size_err = fit_function.GetParError(5);
      result.chi2 = fit_function.GetChisquare();
      result.ndf = fit_function.GetNDF();
      result.fit_statistic = fit_function.GetChisquare();
      result.edm = edm;
      result.status = fit_status;
      result.minuit_istat = covariance_status;
      return result;
    }

    void WriteFitArtifacts(TFile &output_file,
                           const SliceCatalogEntry &entry,
                           TH1D &data_histogram,
                           TF1 &fit_function,
                           TH1D &baseline_histogram,
                           TH1D &pure_femto_histogram,
                           TCanvas &fit_canvas) {
      auto *fit_directory = GetOrCreateDirectoryPath(output_file, BuildFitDirectory(entry.slice_id));
      fit_directory->cd();
      data_histogram.Write("DataCF", TObject::kOverwrite);
      fit_function.Write("FitFunction", TObject::kOverwrite);
      baseline_histogram.Write("Baseline", TObject::kOverwrite);
      pure_femto_histogram.Write("PureFemto", TObject::kOverwrite);
      fit_canvas.Write("FitCanvas", TObject::kOverwrite);
    }

  }  // namespace

  BuildCfRunStatistics RunBuildCf(const ApplicationConfig &config, const Logger &logger) {
    const std::string input_root_path = config.input.input_root;
    const std::string output_root_path = ResolvePath(config.output.output_directory, config.output.cf_root_name);
    EnsureDirectoryExists(config.output.output_directory);
    CreateOrResetRootFile(output_root_path);

    auto input_file = OpenRootFile(input_root_path, "READ");
    auto *same_sparse =
        dynamic_cast<THnSparseF *>(input_file->Get(BuildSparseObjectPath(config.input, config.input.same_event_subtask).c_str()));
    auto *mixed_sparse = dynamic_cast<THnSparseF *>(
        input_file->Get(BuildSparseObjectPath(config.input, config.input.mixed_event_subtask).c_str()));
    if (same_sparse == nullptr || mixed_sparse == nullptr) {
      throw std::runtime_error("Failed to load required same-event or mixed-event THnSparseF objects.");
    }

    std::vector<SliceCatalogEntry> catalog_entries;
    const auto regions = BuildRegionDefinitions();
    BuildCfRunStatistics statistics;
    statistics.requested_groups = config.centrality_bins.size() * config.mt_bins.size();
    ProgressReporter progress(logger, "build-cf", statistics.requested_groups, config.build.progress);

    std::unique_ptr<TFile> shared_output_file;
    if (!config.build.reopen_output_file_per_slice) {
      shared_output_file = OpenRootFile(output_root_path, "UPDATE");
    }

    std::size_t completed_groups = 0;
    for (std::size_t cent_index = 0; cent_index < config.centrality_bins.size(); ++cent_index) {
      for (std::size_t mt_index = 0; mt_index < config.mt_bins.size(); ++mt_index) {
        const RangeBin &cent_bin = config.centrality_bins[cent_index];
        const RangeBin &mt_bin = config.mt_bins[mt_index];
        const std::string group_id = BuildGroupId(cent_bin, mt_bin);

        std::unique_ptr<TH1D> me_projection = BuildProjection(*mixed_sparse, cent_bin, mt_bin, regions[0]);
        const double me_norm = IntegralInRange(*me_projection, config.build.norm_low, config.build.norm_high);
        if (me_norm <= 0.0 || !std::isfinite(me_norm)) {
          ++statistics.skipped_zero_mixed_event_groups;
          progress.Update(++completed_groups);
          continue;
        }

        for (const RegionDefinition &region : regions) {
          std::unique_ptr<TH1D> se_projection = BuildProjection(*same_sparse, cent_bin, mt_bin, region);
          const double se_norm = IntegralInRange(*se_projection, config.build.norm_low, config.build.norm_high);
          if (se_norm <= 0.0 || !std::isfinite(se_norm)) {
            ++statistics.skipped_zero_same_event_slices;
            continue;
          }

          const std::optional<SliceHistograms> histograms =
              BuildCf1D(*se_projection, *me_projection, config.build);
          if (!histograms.has_value()) {
            ++statistics.skipped_zero_same_event_slices;
            continue;
          }

          SliceCatalogEntry entry;
          entry.group_id = group_id;
          entry.slice_id = BuildSliceId(group_id, region);
          entry.slice_directory = BuildSliceDirectory(entry.slice_id);
          entry.se_object_path = entry.slice_directory + "/SE_raw1d";
          entry.me_object_path = entry.slice_directory + "/ME_raw1d";
          entry.cf_object_path = entry.slice_directory + "/CF1D";
          entry.centrality_index = static_cast<int>(cent_index);
          entry.mt_index = static_cast<int>(mt_index);
          entry.region_index = region.region_index;
          entry.cent_low = cent_bin.min;
          entry.cent_high = cent_bin.max;
          entry.mt_low = mt_bin.min;
          entry.mt_high = mt_bin.max;
          entry.region_name = region.name;
          entry.region_kind = region.kind;
          entry.ep_low_1 = region.low_1;
          entry.ep_high_1 = region.high_1;
          entry.ep_low_2 = region.low_2;
          entry.ep_high_2 = region.high_2;
          entry.has_second_interval = region.has_second_interval;
          entry.norm_low = config.build.norm_low;
          entry.norm_high = config.build.norm_high;
          entry.kstar_min = config.build.kstar_min;
          entry.kstar_max = config.build.kstar_max;

          if (config.build.reopen_output_file_per_slice) {
            auto output_file = OpenRootFile(output_root_path, "UPDATE");
            WriteSliceHistograms(*output_file, entry, *histograms);
          } else {
            WriteSliceHistograms(*shared_output_file, entry, *histograms);
          }

          catalog_entries.push_back(entry);
          ++statistics.stored_slices;
        }

        progress.Update(++completed_groups);
      }
    }

    if (shared_output_file == nullptr) {
      shared_output_file = OpenRootFile(output_root_path, "UPDATE");
    }
    WriteSliceCatalogTree(*shared_output_file, catalog_entries);
    progress.Finish();
    return statistics;
  }

  FitRunStatistics RunFit(const ApplicationConfig &config,
                          const Logger &logger,
                          std::optional<std::string> input_cf_root_path) {
    const std::string cf_root_path =
        input_cf_root_path.has_value() ? *input_cf_root_path : ResolvePath(config.output.output_directory, config.output.cf_root_name);
    const std::string fit_root_path = ResolvePath(config.output.output_directory, config.output.fit_root_name);
    const std::string summary_tsv_path = ResolvePath(config.output.output_directory, config.output.fit_summary_name);
    EnsureDirectoryExists(config.output.output_directory);
    CreateOrResetRootFile(fit_root_path);

    auto catalog_entries = LoadSliceCatalog(cf_root_path);
    FitRunStatistics statistics;
    statistics.catalog_slices = catalog_entries.size();

    std::vector<SliceCatalogEntry> selected_entries;
    for (const SliceCatalogEntry &entry : catalog_entries) {
      if (MatchSelectedSlice(entry, config.fit_centrality_bins, config.fit_mt_bins)) {
        selected_entries.push_back(entry);
      }
    }
    statistics.selected_slices = selected_entries.size();

    auto cf_file = OpenRootFile(cf_root_path, "READ");
    std::unique_ptr<TFile> shared_output_file;
    if (!config.fit.reopen_output_file_per_slice) {
      shared_output_file = OpenRootFile(fit_root_path, "UPDATE");
    }

    PiPiCatsModel model(config.fit.cats_num_mom_bins,
                        config.fit.cats_kmin_mev,
                        config.fit.cats_kmax_mev,
                        config.fit.use_coulomb);
    if (1000.0 * config.fit.fit_kstar_max > config.fit.cats_kmax_mev + 1.0e-6) {
      throw std::runtime_error("fit.fit_kstar_max exceeds the configured CATS momentum grid.");
    }

    std::vector<PiPiFitResult> fit_results;
    ProgressReporter progress(logger, "fit", selected_entries.size(), config.fit.progress);
    std::size_t completed_slices = 0;

    for (const SliceCatalogEntry &entry : selected_entries) {
      std::unique_ptr<TH1D> data_histogram = LoadStoredHistogram1D(*cf_file, entry.cf_object_path);
      if (data_histogram == nullptr) {
        ++statistics.skipped_missing_objects;
        progress.Update(++completed_slices);
        continue;
      }

      data_histogram->SetName("DataCF");
      std::unique_ptr<TF1> fit_function(model.BuildFitFunction(entry.slice_id + "_fit", config.fit));
      fit_function->SetRange(0.0, config.fit.fit_kstar_max);

      const TFitResultPtr fit_result_ptr = data_histogram->Fit(fit_function.get(), "QRS0");
      const int fit_status = static_cast<int>(fit_result_ptr);
      const TFitResult *fit_result_object = fit_result_ptr.Get();
      const int covariance_status = fit_result_object != nullptr ? fit_result_object->CovMatrixStatus() : -1;
      const double edm = fit_result_object != nullptr ? fit_result_object->Edm() : std::numeric_limits<double>::quiet_NaN();
      if (fit_status == 0) {
        ++statistics.fitted_slices;
      } else {
        ++statistics.skipped_failed_fits;
      }

      PiPiFitResult result = BuildFitResult(entry, config.fit, *fit_function, fit_status, covariance_status, edm);
      fit_results.push_back(result);

      std::unique_ptr<TH1D> baseline_histogram = BuildBaselineHistogram(*data_histogram, result);
      std::unique_ptr<TH1D> pure_femto_histogram = BuildPureFemtoHistogram(*data_histogram, model, result);

      const bool old_batch_state = gROOT->IsBatch();
      gROOT->SetBatch(true);
      std::unique_ptr<TCanvas> fit_canvas =
          BuildFitCanvas(*data_histogram, *fit_function, *baseline_histogram, *pure_femto_histogram, result);
      gROOT->SetBatch(old_batch_state);

      if (config.fit.reopen_output_file_per_slice) {
        auto output_file = OpenRootFile(fit_root_path, "UPDATE");
        WriteFitArtifacts(*output_file,
                          entry,
                          *data_histogram,
                          *fit_function,
                          *baseline_histogram,
                          *pure_femto_histogram,
                          *fit_canvas);
      } else {
        WriteFitArtifacts(*shared_output_file,
                          entry,
                          *data_histogram,
                          *fit_function,
                          *baseline_histogram,
                          *pure_femto_histogram,
                          *fit_canvas);
      }

      progress.Update(++completed_slices);
    }

    if (shared_output_file == nullptr) {
      shared_output_file = OpenRootFile(fit_root_path, "UPDATE");
    }
    WriteFitCatalogTree(*shared_output_file, fit_results);
    WriteRegionSummaryHistograms(*shared_output_file, fit_results);
    WriteFitResultsSummaryTsv(summary_tsv_path, fit_results);
    progress.Finish();
    return statistics;
  }

  std::vector<SliceCatalogEntry> LoadSliceCatalog(const std::string &cf_root_path) {
    auto file = OpenRootFile(cf_root_path, "READ");
    auto *tree = dynamic_cast<TTree *>(file->Get("meta/SliceCatalog"));
    if (tree == nullptr) {
      throw std::runtime_error("SliceCatalog is missing from CF ROOT file: " + cf_root_path);
    }
    return ReadSliceCatalogTree(*tree);
  }

}  // namespace exp_femto_1d
