#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "TFile.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
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
                          const std::string &output_dir,
                          const std::string &cf_root_name,
                          const std::string &fit_root_name,
                          const std::string &fit_summary_name,
                          const bool build_map_pair_phi_to_symmetric_range,
                          const std::optional<bool> fit_map_pair_phi_to_symmetric_range) {
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
    output << "map_pair_phi_to_symmetric_range = " << (build_map_pair_phi_to_symmetric_range ? "true" : "false") << "\n";
    output << "write_normalized_se_me_1d_projections = true\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "progress = false\n\n";
    output << "[fit]\n";
    output << "model = \"full\"\n";
    output << "use_coulomb = false\n";
    output << "use_core_halo_lambda = true\n";
    output << "use_q2_baseline = false\n";
    output << "use_pml = false\n";
    output << "fit_q_max = 0.12\n";
    if (fit_map_pair_phi_to_symmetric_range.has_value()) {
      output << "map_pair_phi_to_symmetric_range = "
             << (*fit_map_pair_phi_to_symmetric_range ? "true" : "false") << "\n";
    }
    output << "reopen_output_file_per_slice = false\n";
    output << "progress = false\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    return path.string();
  }

  struct FitCatalogInspection {
    double min_phi = 0.0;
    double max_phi = 0.0;
    bool fit_uses_symmetric_phi_range = false;
    bool saw_non_integrated_slice = false;
  };

  FitCatalogInspection InspectFitCatalog(const std::filesystem::path &fit_root_path) {
    TFile fit_file(fit_root_path.string().c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(fit_file.Get("meta/FitCatalog"));
    Expect(tree != nullptr, "FitCatalog missing");

    TTreeReader reader(tree);
    TTreeReaderValue<double> phi(reader, "phi");
    TTreeReaderValue<int> is_phi_integrated(reader, "is_phi_integrated");
    TTreeReaderValue<int> fit_uses_symmetric_phi_range(reader, "fit_uses_symmetric_phi_range");

    FitCatalogInspection inspection;
    while (reader.Next()) {
      if (*is_phi_integrated != 0) {
        continue;
      }
      if (!inspection.saw_non_integrated_slice) {
        inspection.min_phi = *phi;
        inspection.max_phi = *phi;
        inspection.fit_uses_symmetric_phi_range = (*fit_uses_symmetric_phi_range != 0);
        inspection.saw_non_integrated_slice = true;
        continue;
      }
      inspection.min_phi = std::min(inspection.min_phi, *phi);
      inspection.max_phi = std::max(inspection.max_phi, *phi);
      Expect(inspection.fit_uses_symmetric_phi_range == (*fit_uses_symmetric_phi_range != 0),
             "fit phi mapping flag should be stable across slices");
    }

    Expect(inspection.saw_non_integrated_slice, "FitCatalog should contain phi-differential slices");
    return inspection;
  }

  std::pair<double, double> InspectSummaryPhiFitRange(const std::filesystem::path &fit_root_path) {
    TFile fit_file(fit_root_path.string().c_str(), "READ");
    auto *fit_function =
        dynamic_cast<TF1 *>(fit_file.Get("summary/R2_vs_phi/cent_0.00-10.00__mt_0.20-0.40/Rout2_phi_fit"));
    Expect(fit_function != nullptr, "Rout2_phi_fit summary function missing");
    return {fit_function->GetXmin(), fit_function->GetXmax()};
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_3d_workflow_smoke";
  std::filesystem::create_directories(temp_dir);
  const std::string input_root = WriteToyInput(temp_dir / "input.root");
  const Logger logger(LogLevel::kError);

  const std::string mapped_follow_config_path = WriteConfig(temp_dir / "mapped_follow.toml",
                                                            input_root,
                                                            temp_dir.string(),
                                                            "mapped_cf.root",
                                                            "mapped_follow_fit.root",
                                                            "mapped_follow.tsv",
                                                            true,
                                                            std::nullopt);
  const ApplicationConfig mapped_follow_config = LoadApplicationConfig(mapped_follow_config_path);
  const BuildCfRunStatistics mapped_build_stats = RunBuildCf(mapped_follow_config, logger);
  Expect(mapped_build_stats.stored_slices == 4, "mapped build-cf should produce 4 slices");
  const FitRunStatistics mapped_follow_fit_stats = RunFit(mapped_follow_config, logger);
  Expect(mapped_follow_fit_stats.selected_slices == 4, "follow-input fit should select every built slice");

  TFile mapped_cf_file((temp_dir / "mapped_cf.root").string().c_str(), "READ");
  Expect(mapped_cf_file.Get("meta/SliceCatalog") != nullptr, "mapped SliceCatalog missing");
  Expect(mapped_cf_file.Get("slices") != nullptr, "mapped slices directory missing");

  const FitCatalogInspection mapped_follow_inspection = InspectFitCatalog(temp_dir / "mapped_follow_fit.root");
  Expect(mapped_follow_inspection.fit_uses_symmetric_phi_range,
         "fit should follow the mapped CF metadata when no fit override is given");
  Expect(mapped_follow_inspection.min_phi < 0.0, "mapped follow fit should keep a negative phi slice");
  const auto mapped_follow_phi_range = InspectSummaryPhiFitRange(temp_dir / "mapped_follow_fit.root");
  Expect(std::abs(mapped_follow_phi_range.first + TMath::Pi() / 2.0) < 1.0e-6,
         "mapped follow phi fit minimum should be -pi/2");
  Expect(std::abs(mapped_follow_phi_range.second - TMath::Pi() / 2.0) < 1.0e-6,
         "mapped follow phi fit maximum should be pi/2");

  const std::string mapped_override_raw_config_path = WriteConfig(temp_dir / "mapped_override_raw.toml",
                                                                  input_root,
                                                                  temp_dir.string(),
                                                                  "mapped_cf.root",
                                                                  "mapped_override_raw_fit.root",
                                                                  "mapped_override_raw.tsv",
                                                                  false,
                                                                  false);
  const ApplicationConfig mapped_override_raw_config = LoadApplicationConfig(mapped_override_raw_config_path);
  const FitRunStatistics mapped_override_raw_fit_stats = RunFit(mapped_override_raw_config, logger);
  Expect(mapped_override_raw_fit_stats.selected_slices == 4,
         "raw override fit should select every slice from the mapped CF");
  const FitCatalogInspection mapped_override_raw_inspection =
      InspectFitCatalog(temp_dir / "mapped_override_raw_fit.root");
  Expect(!mapped_override_raw_inspection.fit_uses_symmetric_phi_range,
         "fit override should switch the mapped CF back to raw phi semantics");
  Expect(mapped_override_raw_inspection.min_phi >= 0.0, "raw override fit should not contain negative phi slices");
  Expect(mapped_override_raw_inspection.max_phi > TMath::Pi() / 2.0,
         "raw override fit should expose the original high-phi slice");
  const auto mapped_override_raw_phi_range = InspectSummaryPhiFitRange(temp_dir / "mapped_override_raw_fit.root");
  Expect(std::abs(mapped_override_raw_phi_range.first - 0.0) < 1.0e-6,
         "raw override phi fit minimum should be 0");
  Expect(std::abs(mapped_override_raw_phi_range.second - TMath::Pi()) < 1.0e-6,
         "raw override phi fit maximum should be pi");

  const std::string raw_override_mapped_config_path = WriteConfig(temp_dir / "raw_override_mapped.toml",
                                                                  input_root,
                                                                  temp_dir.string(),
                                                                  "raw_cf.root",
                                                                  "raw_override_mapped_fit.root",
                                                                  "raw_override_mapped.tsv",
                                                                  false,
                                                                  true);
  const ApplicationConfig raw_override_mapped_config = LoadApplicationConfig(raw_override_mapped_config_path);
  const BuildCfRunStatistics raw_build_stats = RunBuildCf(raw_override_mapped_config, logger);
  Expect(raw_build_stats.stored_slices == 4, "raw build-cf should produce 4 slices");
  const FitRunStatistics raw_override_mapped_fit_stats = RunFit(raw_override_mapped_config, logger);
  Expect(raw_override_mapped_fit_stats.selected_slices == 4,
         "mapped override fit should select every slice from the raw CF");
  const FitCatalogInspection raw_override_mapped_inspection =
      InspectFitCatalog(temp_dir / "raw_override_mapped_fit.root");
  Expect(raw_override_mapped_inspection.fit_uses_symmetric_phi_range,
         "fit override should remap raw CF phi coordinates into the symmetric range");
  Expect(raw_override_mapped_inspection.min_phi < 0.0, "mapped override fit should contain a negative phi slice");
  const auto raw_override_mapped_phi_range = InspectSummaryPhiFitRange(temp_dir / "raw_override_mapped_fit.root");
  Expect(std::abs(raw_override_mapped_phi_range.first + TMath::Pi() / 2.0) < 1.0e-6,
         "mapped override phi fit minimum should be -pi/2");
  Expect(std::abs(raw_override_mapped_phi_range.second - TMath::Pi() / 2.0) < 1.0e-6,
         "mapped override phi fit maximum should be pi/2");

  std::ifstream summary((temp_dir / "mapped_follow.tsv").string());
  std::string header;
  std::getline(summary, header);
  Expect(header.find("sliceId") != std::string::npos, "summary TSV header missing");

  std::cout << "workflow_smoke_test passed\n";
  return 0;
}
