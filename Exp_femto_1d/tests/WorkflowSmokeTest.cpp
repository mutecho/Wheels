#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "THnSparse.h"
#include "exp_femto_1d/Config.h"
#include "exp_femto_1d/Logging.h"
#include "exp_femto_1d/Workflow.h"

namespace {

  void Expect(const bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void FillSparse(THnSparseF &sparse,
                  const double kstar,
                  const double mt,
                  const double centrality,
                  const double ep,
                  const double weight) {
    double values[4] = {kstar, mt, centrality, ep};
    sparse.Fill(values, weight);
  }

  std::string WriteToyInput(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    auto *task = output.mkdir("task");
    auto *same_dir = task->mkdir("Same");
    auto *mixed_dir = task->mkdir("Mixed");

    const int bins[4] = {48, 8, 8, 16};
    const double min[4] = {0.0, 0.2, 0.0, 0.0};
    const double max[4] = {0.8, 0.5, 10.0, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 4, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 4, bins, min, max);
    for (double kstar : {0.02, 0.04, 0.08, 0.12, 0.20, 0.55, 0.70}) {
      for (double mt : {0.30, 0.42}) {
        for (double ep : {0.10, 0.30, 1.30, 2.60}) {
          const double same_weight = 25.0 + 10.0 * std::exp(-7.0 * kstar) + (ep < 0.8 ? 3.0 : 0.0)
                                     + 2.0 * (mt > 0.35 ? 1.0 : 0.0);
          const double mixed_weight =
              22.0 + 4.0 * std::exp(-2.0 * kstar) + (mt > 0.35 ? 1.5 : 0.0) + 5.0 * ep;
          FillSparse(*same, kstar, mt, 5.0, ep, same_weight);
          FillSparse(*mixed, kstar, mt, 5.0, ep, mixed_weight);
        }
      }
    }

    same_dir->cd();
    same->Write("sparse");
    mixed_dir->cd();
    mixed->Write("sparse");
    output.Close();
    return path.string();
  }

  std::string WriteConfig(const std::filesystem::path &path,
                          const std::string &input_root,
                          const std::string &output_dir,
                          const std::string &cf_root_name = "workflow_cf.root",
                          const std::string &fit_root_name = "workflow_fit.root",
                          const std::string &fit_summary_name = "workflow.tsv",
                          const bool split_mixed_event_by_phi = false,
                          const unsigned cf_rebin_factor = 2U) {
    std::ofstream output(path);
    output << "[input]\n";
    output << "input_root = \"" << input_root << "\"\n";
    output << "task_name = \"task\"\n";
    output << "same_event_subtask = \"Same\"\n";
    output << "mixed_event_subtask = \"Mixed\"\n";
    output << "sparse_object_name = \"sparse\"\n\n";
    output << "[output]\n";
    output << "output_directory = \"" << output_dir << "\"\n";
    output << "cf_root_name = \"" << cf_root_name << "\"\n";
    output << "fit_root_name = \"" << fit_root_name << "\"\n";
    output << "fit_summary_name = \"" << fit_summary_name << "\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "norm_low = 0.5\n";
    output << "norm_high = 0.8\n";
    output << "kstar_min = 0.0\n";
    output << "kstar_max = 0.8\n";
    output << "cf_rebin_factor = " << cf_rebin_factor << "\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "cf_by_mt_show_markers = true\n";
    output << "split_mixed_event_by_phi = " << (split_mixed_event_by_phi ? "true" : "false") << "\n";
    output << "progress = false\n\n";
    output << "[fit]\n";
    output << "fit_kstar_max = 0.20\n";
    output << "use_coulomb = false\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "progress = false\n";
    output << "baseline_p0_init = 1.0\n";
    output << "baseline_p1_init = 0.0\n";
    output << "baseline_p2_init = 0.0\n";
    output << "source_size_init = 6.0\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.35\n\n";
    output << "[[bins.mt]]\nmin = 0.35\nmax = 0.5\n\n";
    output << "[[fit_selection.mt]]\nmin = 0.2\nmax = 0.35\n";
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_1d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_1d_workflow_smoke";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const std::string config_path = WriteConfig(temp_dir / "config.toml", input_root, temp_dir.string());
  const std::string split_config_path = WriteConfig(temp_dir / "config_split.toml",
                                                    input_root,
                                                    temp_dir.string(),
                                                    "workflow_cf_split.root",
                                                    "workflow_fit_split.root",
                                                    "workflow_split.tsv",
                                                    true,
                                                    2U);

  const ApplicationConfig config = LoadApplicationConfig(config_path);
  const Logger logger(LogLevel::kError);
  const BuildCfRunStatistics build_stats = RunBuildCf(config, logger);
  Expect(build_stats.requested_groups == 2, "expected two toy groups");
  Expect(build_stats.stored_slices == 6, "build-cf should produce three region slices for each mT bin");

  TFile cf_file((temp_dir / "workflow_cf.root").string().c_str(), "READ");
  Expect(cf_file.Get("meta/SliceCatalog") != nullptr, "SliceCatalog missing");
  Expect(cf_file.Get("slices") != nullptr, "slices directory missing");

  const auto entries = LoadSliceCatalog((temp_dir / "workflow_cf.root").string());
  Expect(entries.size() == 6, "catalog size mismatch");
  auto *se_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].se_object_path.c_str()));
  auto *me_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].me_object_path.c_str()));
  auto *cf_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].cf_object_path.c_str()));
  auto *min_bias_canvas = dynamic_cast<TCanvas *>(cf_file.Get("cent_slices/cent_0.00-10.00/MinBias/CFByMtCanvas"));
  auto *in_plane_canvas = dynamic_cast<TCanvas *>(cf_file.Get("cent_slices/cent_0.00-10.00/InPlane/CFByMtCanvas"));
  auto *out_of_plane_canvas = dynamic_cast<TCanvas *>(cf_file.Get("cent_slices/cent_0.00-10.00/OutOfPlane/CFByMtCanvas"));
  Expect(se_histogram != nullptr, "SE_raw1d missing");
  Expect(me_histogram != nullptr, "ME_raw1d missing");
  Expect(cf_histogram != nullptr, "CF1D missing");
  Expect(entries[0].cf_rebin_factor == 2, "catalog should record CF rebin factor");
  Expect(cf_histogram->GetNbinsX() * 2 == se_histogram->GetNbinsX(), "CF1D should use rebinned SE/ME binning");
  Expect(cf_histogram->GetNbinsX() * 2 == me_histogram->GetNbinsX(), "CF1D should be rebinned relative to raw ME");
  Expect(cf_histogram->Integral() > 0.0, "CF1D should not be empty");
  Expect(min_bias_canvas != nullptr, "MinBias CF-by-mT canvas missing");
  Expect(in_plane_canvas != nullptr, "InPlane CF-by-mT canvas missing");
  Expect(out_of_plane_canvas != nullptr, "OutOfPlane CF-by-mT canvas missing");
  auto *enabled_marker_graph =
      dynamic_cast<TGraph *>(min_bias_canvas->GetListOfPrimitives()->FindObject("CFTrend__mt_0.200-0.350"));
  Expect(enabled_marker_graph != nullptr, "CF-by-mT trend graph missing");
  Expect(enabled_marker_graph->GetMarkerSize() > 0.0, "CF-by-mT marker switch should enable markers");

  const ApplicationConfig split_config = LoadApplicationConfig(split_config_path);
  const BuildCfRunStatistics split_build_stats = RunBuildCf(split_config, logger);
  Expect(split_build_stats.stored_slices == 6, "split-ME build-cf should keep every toy slice");
  Expect(split_build_stats.skipped_zero_mixed_event_slices == 0, "toy split-ME build should not skip ME slices");
  const auto split_entries = LoadSliceCatalog((temp_dir / "workflow_cf_split.root").string());
  Expect(split_entries.size() == 6, "split-ME catalog size mismatch");
  Expect(split_entries[1].split_mixed_event_by_phi, "catalog should record split-ME mode");

  TFile split_cf_file((temp_dir / "workflow_cf_split.root").string().c_str(), "READ");
  auto *split_in_plane_me = dynamic_cast<TH1D *>(split_cf_file.Get(split_entries[1].me_object_path.c_str()));
  auto *integrated_in_plane_me = dynamic_cast<TH1D *>(cf_file.Get(entries[1].me_object_path.c_str()));
  Expect(split_in_plane_me != nullptr && integrated_in_plane_me != nullptr, "ME comparison histograms missing");
  Expect(split_in_plane_me->Integral() < integrated_in_plane_me->Integral(),
         "split-ME InPlane denominator should use a narrower EP range than integrated ME");

  const FitRunStatistics fit_stats = RunFit(config, logger);
  Expect(fit_stats.catalog_slices == 6, "fit should read six catalog slices");
  Expect(fit_stats.selected_slices == 3, "fit should select the configured mT bin");
  Expect(fit_stats.fitted_slices + fit_stats.skipped_failed_fits == 3, "every selected slice should be attempted");

  TFile fit_file((temp_dir / "workflow_fit.root").string().c_str(), "READ");
  Expect(fit_file.Get("meta/FitCatalog") != nullptr, "FitCatalog missing");
  Expect(fit_file.Get("summary/by_region") != nullptr, "summary/by_region missing");
  Expect(fit_file.Get(("fits/" + entries[0].slice_id + "/DataCF").c_str()) != nullptr, "DataCF missing");
  Expect(fit_file.Get(("fits/" + entries[0].slice_id + "/FitFunction").c_str()) != nullptr, "FitFunction missing");

  std::ifstream summary((temp_dir / "workflow.tsv").string());
  std::string header;
  std::getline(summary, header);
  Expect(header.find("sliceId") != std::string::npos, "summary TSV header missing");
  Expect(header.find("covarianceQuality") != std::string::npos, "summary TSV covariance column missing");

  return 0;
}
