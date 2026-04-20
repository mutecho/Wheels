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

    const int bins[4] = {32, 8, 8, 16};
    const double min[4] = {0.0, 0.2, 0.0, 0.0};
    const double max[4] = {0.8, 0.5, 10.0, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 4, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 4, bins, min, max);
    for (double kstar : {0.05, 0.10, 0.55, 0.70}) {
      for (double ep : {0.10, 0.30, 1.30, 2.60}) {
        FillSparse(*same, kstar, 0.30, 5.0, ep, 12.0 + 30.0 * std::exp(-6.0 * kstar));
        FillSparse(*mixed, kstar, 0.30, 5.0, ep, 10.0 + 12.0 * std::exp(-2.0 * kstar));
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
    output << "cf_root_name = \"catalog_test.root\"\n";
    output << "fit_root_name = \"catalog_fit.root\"\n";
    output << "fit_summary_name = \"catalog.tsv\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "norm_low = 0.5\n";
    output << "norm_high = 0.8\n";
    output << "kstar_min = 0.0\n";
    output << "kstar_max = 0.8\n";
    output << "reopen_output_file_per_slice = true\n";
    output << "progress = false\n\n";
    output << "[fit]\n";
    output << "fit_kstar_max = 0.20\n";
    output << "use_coulomb = false\n";
    output << "progress = false\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_1d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_1d_catalog_test";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const std::string config_path = WriteConfig(temp_dir / "config.toml", input_root, temp_dir.string());

  const ApplicationConfig config = LoadApplicationConfig(config_path);
  const Logger logger(LogLevel::kError);
  const BuildCfRunStatistics build_stats = RunBuildCf(config, logger);
  Expect(build_stats.stored_slices == 3, "expected three region slices");
  Expect(build_stats.skipped_zero_mixed_event_groups == 0, "toy input should keep the group");
  Expect(build_stats.skipped_zero_same_event_slices == 0, "toy input should keep every slice");

  const auto entries = LoadSliceCatalog((temp_dir / "catalog_test.root").string());
  Expect(entries.size() == 3, "catalog size mismatch");
  Expect(entries.front().slice_directory.rfind("slices/", 0) == 0, "slice directory should live under slices/");
  Expect(entries.front().norm_low == 0.5, "catalog should persist norm_low");
  Expect(entries.front().kstar_max == 0.8, "catalog should persist kstar_max");
  Expect(entries[0].region_name == "MinBias", "first region should be MinBias");
  Expect(entries[1].has_second_interval, "InPlane should keep the split EP interval metadata");
  Expect(entries[2].region_name == "OutOfPlane", "third region should be OutOfPlane");
  Expect(entries[0].cf_object_path.find("/CF1D") != std::string::npos, "catalog should store CF object path");

  TFile cf_file((temp_dir / "catalog_test.root").string().c_str(), "READ");
  auto *catalog_tree = cf_file.Get("meta/SliceCatalog");
  auto *cf_histogram = dynamic_cast<TH1D *>(cf_file.Get(entries[0].cf_object_path.c_str()));
  Expect(catalog_tree != nullptr, "SliceCatalog tree missing");
  Expect(cf_histogram != nullptr, "CF1D histogram missing");
  Expect(cf_histogram->GetNbinsX() > 0, "CF1D histogram should have bins");

  return 0;
}
