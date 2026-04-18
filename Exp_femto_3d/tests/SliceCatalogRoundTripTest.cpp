#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

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

    const int bins[7] = {4, 4, 4, 2, 2, 1, 3};
    const double min[7] = {-0.2, -0.2, -0.2, 0.2, 0.0, -0.5, 0.0};
    const double max[7] = {0.2, 0.2, 0.2, 0.4, 10.0, 0.5, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    for (double phi : {0.3, 1.3, 2.5}) {
      FillSparse(*same, 0.01, 0.01, 0.01, 0.3, 5.0, phi, 10.0);
      FillSparse(*same, -0.01, 0.01, -0.01, 0.3, 5.0, phi, 8.0);
      FillSparse(*mixed, 0.01, 0.01, 0.01, 0.3, 5.0, phi, 12.0);
      FillSparse(*mixed, -0.01, 0.01, -0.01, 0.3, 5.0, phi, 11.0);
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
    output << "map_pair_phi_to_symmetric_range = true\n";
    output << "write_normalized_se_me_1d_projections = false\n";
    output << "reopen_output_file_per_slice = true\n";
    output << "progress = false\n\n";
    output << "[fit]\n";
    output << "model = \"diag\"\n";
    output << "use_coulomb = false\n";
    output << "use_core_halo_lambda = true\n";
    output << "use_q2_baseline = false\n";
    output << "use_pml = false\n";
    output << "fit_q_max = 0.15\n";
    output << "progress = false\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    return path.string();
  }

  // Write a legacy SliceCatalog tree without the build phi mapping branch so the
  // loader must infer the old build semantics from raw/display phi coordinates.
  std::string WriteLegacyCatalog(const std::filesystem::path &path, const std::vector<exp_femto_3d::SliceCatalogEntry> &entries) {
    TFile output(path.string().c_str(), "RECREATE");
    auto *meta_directory = output.mkdir("meta");
    meta_directory->cd();

    {
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

      for (const exp_femto_3d::SliceCatalogEntry &entry : entries) {
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
    }
    output.Close();
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_3d_catalog_test";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const std::string config_path = WriteConfig(temp_dir / "config.toml", input_root, temp_dir.string());

  const ApplicationConfig config = LoadApplicationConfig(config_path);
  const Logger logger(LogLevel::kError);
  const BuildCfRunStatistics build_stats = RunBuildCf(config, logger);
  Expect(build_stats.stored_slices == 4, "expected phi-all + 3 phi slices");

  const auto entries = LoadSliceCatalog((temp_dir / "catalog_test.root").string());
  Expect(entries.size() == 4, "catalog size mismatch");
  Expect(entries.front().slice_directory.rfind("slices/", 0) == 0, "slice directory should live under slices/");
  Expect(entries[0].cf_object_path.find("/CF3D") != std::string::npos, "catalog should store CF object path");
  Expect(entries[1].build_uses_symmetric_phi_range, "catalog should persist build phi mapping metadata");
  Expect(entries[3].display_phi_center < 0.0, "mapped phi slice should keep a negative display center");

  const auto legacy_entries = LoadSliceCatalog(WriteLegacyCatalog(temp_dir / "legacy_catalog.root", entries));
  Expect(legacy_entries.size() == entries.size(), "legacy catalog size mismatch");
  Expect(legacy_entries[1].build_uses_symmetric_phi_range,
         "legacy catalog reader should infer mapped build phi metadata");
  Expect(legacy_entries[3].display_phi_center < 0.0, "legacy catalog should preserve stored display phi values");

  std::cout << "slice_catalog_roundtrip_test passed\n";
  return 0;
}
