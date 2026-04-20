#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "TFile.h"
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
      for (double ep : {0.10, 0.30, 1.30, 2.60}) {
        const double same_weight = 25.0 + 10.0 * std::exp(-7.0 * kstar) + (ep < 0.8 ? 3.0 : 0.0);
        const double mixed_weight = 22.0 + 4.0 * std::exp(-2.0 * kstar);
        FillSparse(*same, kstar, 0.30, 5.0, ep, same_weight);
        FillSparse(*mixed, kstar, 0.30, 5.0, ep, mixed_weight);
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
                          const std::string &output_dir) {
    std::ofstream output(path);
    output << "[input]\n";
    output << "input_root = \"" << input_root << "\"\n";
    output << "task_name = \"task\"\n";
    output << "same_event_subtask = \"Same\"\n";
    output << "mixed_event_subtask = \"Mixed\"\n";
    output << "sparse_object_name = \"sparse\"\n\n";
    output << "[output]\n";
    output << "output_directory = \"" << output_dir << "\"\n";
    output << "cf_root_name = \"workflow_cf.root\"\n";
    output << "fit_root_name = \"workflow_fit.root\"\n";
    output << "fit_summary_name = \"workflow.tsv\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "norm_low = 0.5\n";
    output << "norm_high = 0.8\n";
    output << "kstar_min = 0.0\n";
    output << "kstar_max = 0.8\n";
    output << "reopen_output_file_per_slice = false\n";
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
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_1d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_1d_workflow_smoke";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const std::string config_path = WriteConfig(temp_dir / "config.toml", input_root, temp_dir.string());

  const ApplicationConfig config = LoadApplicationConfig(config_path);
  const Logger logger(LogLevel::kError);
  const BuildCfRunStatistics build_stats = RunBuildCf(config, logger);
  Expect(build_stats.requested_groups == 1, "expected one toy group");
  Expect(build_stats.stored_slices == 3, "build-cf should produce three region slices");

  TFile cf_file((temp_dir / "workflow_cf.root").string().c_str(), "READ");
  Expect(cf_file.Get("meta/SliceCatalog") != nullptr, "SliceCatalog missing");
  Expect(cf_file.Get("slices") != nullptr, "slices directory missing");

  const auto entries = LoadSliceCatalog((temp_dir / "workflow_cf.root").string());
  Expect(entries.size() == 3, "catalog size mismatch");
  auto *se_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].se_object_path.c_str()));
  auto *me_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].me_object_path.c_str()));
  auto *cf_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].cf_object_path.c_str()));
  Expect(se_histogram != nullptr, "SE_raw1d missing");
  Expect(me_histogram != nullptr, "ME_raw1d missing");
  Expect(cf_histogram != nullptr, "CF1D missing");
  Expect(cf_histogram->Integral() > 0.0, "CF1D should not be empty");

  const FitRunStatistics fit_stats = RunFit(config, logger);
  Expect(fit_stats.catalog_slices == 3, "fit should read three catalog slices");
  Expect(fit_stats.selected_slices == 3, "fit should select every built slice");
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
