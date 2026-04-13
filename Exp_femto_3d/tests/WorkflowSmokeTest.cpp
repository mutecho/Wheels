#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "TFile.h"
#include "THnSparse.h"
#include "TTree.h"
#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Logging.h"
#include "exp_femto_3d/Workflow.h"

namespace {

  void Expect(bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void FillSparse(THnSparseF &sparse,
                  const double q_out,
                  const double q_side,
                  const double q_long,
                  const double mt,
                  const double centrality,
                  const double phi,
                  const double weight) {
    double values[7] = {q_out, q_side, q_long, mt, centrality, 0.0, phi};
    sparse.Fill(values, weight);
  }

  std::string WriteToyInput(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    auto *task = output.mkdir("task");
    auto *same_dir = task->mkdir("Same");
    auto *mixed_dir = task->mkdir("Mixed");

    const int bins[7] = {6, 6, 6, 2, 2, 1, 3};
    const double min[7] = {-0.15, -0.15, -0.15, 0.2, 0.0, -0.5, 0.0};
    const double max[7] = {0.15, 0.15, 0.15, 0.4, 10.0, 0.5, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    for (double phi : {0.3, 1.3, 2.5}) {
      FillSparse(*same, 0.01, 0.00, 0.01, 0.3, 5.0, phi, 50.0);
      FillSparse(*same, -0.01, 0.01, -0.01, 0.3, 5.0, phi, 42.0);
      FillSparse(*same, 0.02, -0.01, 0.00, 0.3, 5.0, phi, 35.0);
      FillSparse(*mixed, 0.01, 0.00, 0.01, 0.3, 5.0, phi, 55.0);
      FillSparse(*mixed, -0.01, 0.01, -0.01, 0.3, 5.0, phi, 48.0);
      FillSparse(*mixed, 0.02, -0.01, 0.00, 0.3, 5.0, phi, 40.0);
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
    output << "cf_root_name = \"smoke_cf.root\"\n";
    output << "fit_root_name = \"smoke_fit.root\"\n";
    output << "fit_summary_name = \"smoke_fit.tsv\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "map_pair_phi_to_symmetric_range = true\n";
    output << "write_normalized_se_me_1d_projections = true\n";
    output << "reopen_output_file_per_slice = false\n\n";
    output << "[fit]\n";
    output << "model = \"full\"\n";
    output << "use_coulomb = false\n";
    output << "use_core_halo_lambda = true\n";
    output << "use_q2_baseline = false\n";
    output << "use_pml = false\n";
    output << "fit_q_max = 0.12\n";
    output << "reopen_output_file_per_slice = false\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_3d_workflow_smoke";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const std::string config_path = WriteConfig(temp_dir / "config.toml", input_root, temp_dir.string());

  const ApplicationConfig config = LoadApplicationConfig(config_path);
  const Logger logger(LogLevel::kError);

  const BuildCfRunStatistics build_stats = RunBuildCf(config, logger);
  Expect(build_stats.stored_slices == 4, "build-cf should produce 4 slices");

  const FitRunStatistics fit_stats = RunFit(config, logger);
  Expect(fit_stats.selected_slices == 4, "fit should select every built slice");

  TFile cf_file((temp_dir / "smoke_cf.root").string().c_str(), "READ");
  Expect(cf_file.Get("meta/SliceCatalog") != nullptr, "SliceCatalog missing");
  Expect(cf_file.Get("slices") != nullptr, "slices directory missing");

  TFile fit_file((temp_dir / "smoke_fit.root").string().c_str(), "READ");
  Expect(fit_file.Get("meta/FitCatalog") != nullptr, "FitCatalog missing");
  Expect(fit_file.Get("summary/R2_vs_phi") != nullptr, "R2_vs_phi summary missing");

  std::ifstream summary((temp_dir / "smoke_fit.tsv").string());
  std::string header;
  std::getline(summary, header);
  Expect(header.find("sliceId") != std::string::npos, "summary TSV header missing");

  std::cout << "workflow_smoke_test passed\n";
  return 0;
}
