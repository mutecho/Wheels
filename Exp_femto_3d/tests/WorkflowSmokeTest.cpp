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

    const int bins[7] = {6, 6, 6, 2, 2, 1, 4};
    const double min[7] = {-0.15, -0.15, -0.15, 0.2, 0.0, -0.5, 0.0};
    const double max[7] = {0.15, 0.15, 0.15, 0.4, 10.0, 0.5, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    for (double phi : {0.3, 1.1, 2.0, 2.8}) {
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
                          const std::string &fit_report_root_name,
                          const bool build_map_pair_phi_to_symmetric_range,
                          const std::optional<bool> fit_map_pair_phi_to_symmetric_range,
                          const std::string &fit_coulomb_lines = "use_coulomb = false\n",
                          const bool use_pml = false,
                          const std::string &fit_parameter_lines = "") {
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
    output << "fit_report_directory = \"" << output_dir << "\"\n";
    output << "fit_report_root_name = \"" << fit_report_root_name << "\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "map_pair_phi_to_symmetric_range = " << (build_map_pair_phi_to_symmetric_range ? "true" : "false") << "\n";
    output << "write_normalized_se_me_1d_projections = true\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "progress = false\n\n";
    output << "[fit]\n";
    output << "model = \"full\"\n";
    output << fit_coulomb_lines;
    output << "use_core_halo_lambda = true\n";
    output << "use_q2_baseline = false\n";
    output << "use_pml = " << (use_pml ? "true" : "false") << "\n";
    output << "fit_q_max = 0.12\n";
    if (fit_map_pair_phi_to_symmetric_range.has_value()) {
      output << "map_pair_phi_to_symmetric_range = "
             << (*fit_map_pair_phi_to_symmetric_range ? "true" : "false") << "\n";
    }
    output << "reopen_output_file_per_slice = false\n";
    output << "progress = false\n\n";
    output << fit_parameter_lines;
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

  FitCatalogInspection InspectFitCatalog(const std::filesystem::path &fit_root_path,
                                         const std::string &expected_coulomb_mode = "none",
                                         const std::optional<double> expected_lambda = std::nullopt,
                                         const std::optional<double> expected_alpha = std::nullopt) {
    TFile fit_file(fit_root_path.string().c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(fit_file.Get("meta/FitCatalog"));
    Expect(tree != nullptr, "FitCatalog missing");
    Expect(tree->GetBranch("coulomb_mode") != nullptr, "FitCatalog coulomb_mode branch missing");
    Expect(tree->GetBranch("coulombMode") != nullptr, "FitCatalog coulombMode branch missing");
    Expect(tree->GetBranch("finite_source_mode") != nullptr, "FitCatalog finite_source_mode branch missing");
    Expect(tree->GetBranch("finiteSourceMode") != nullptr, "FitCatalog finiteSourceMode branch missing");
    Expect(tree->GetBranch("finite_source_radius_fm") != nullptr,
           "FitCatalog finite_source_radius_fm branch missing");
    Expect(tree->GetBranch("finiteSourceRadiusFm") != nullptr,
           "FitCatalog finiteSourceRadiusFm branch missing");
    Expect(tree->GetBranch("qn_index") != nullptr, "FitCatalog qn_index branch missing");
    Expect(tree->GetBranch("qn_label") != nullptr, "FitCatalog qn_label branch missing");
    Expect(tree->GetBranch("is_qn_integrated") != nullptr, "FitCatalog is_qn_integrated branch missing");

    TTreeReader reader(tree);
    TTreeReaderValue<double> phi(reader, "phi");
    TTreeReaderValue<int> uses_coulomb(reader, "uses_coulomb");
    TTreeReaderValue<std::string> coulomb_mode(reader, "coulomb_mode");
    TTreeReaderValue<std::string> finite_source_mode(reader, "finite_source_mode");
    TTreeReaderValue<double> finite_source_radius_fm(reader, "finite_source_radius_fm");
    TTreeReaderValue<int> qn_index(reader, "qn_index");
    TTreeReaderValue<std::string> qn_label(reader, "qn_label");
    TTreeReaderValue<int> is_qn_integrated(reader, "is_qn_integrated");
    TTreeReaderValue<int> is_phi_integrated(reader, "is_phi_integrated");
    TTreeReaderValue<int> fit_uses_symmetric_phi_range(reader, "fit_uses_symmetric_phi_range");
    TTreeReaderValue<double> lambda(reader, "lambda");
    TTreeReaderValue<double> lambda_err(reader, "lambda_err");
    TTreeReaderValue<double> alpha(reader, "alpha");
    TTreeReaderValue<double> alpha_err(reader, "alpha_err");

    FitCatalogInspection inspection;
    while (reader.Next()) {
      Expect(*coulomb_mode == expected_coulomb_mode, "FitCatalog Coulomb mode mismatch");
      Expect((*uses_coulomb != 0) == (expected_coulomb_mode != "none"), "FitCatalog uses_coulomb mismatch");
      Expect(*qn_index == -1 && *qn_label == "qn_all" && *is_qn_integrated != 0,
             "default smoke FitCatalog should remain qn-integrated");
      if (expected_lambda.has_value()) {
        Expect(std::abs(*lambda - *expected_lambda) < 1.0e-9, "fixed lambda value mismatch");
        Expect(std::abs(*lambda_err) < 1.0e-9, "fixed lambda error should be zero");
      }
      if (expected_alpha.has_value()) {
        Expect(std::abs(*alpha - *expected_alpha) < 1.0e-9, "fixed alpha value mismatch");
        Expect(std::abs(*alpha_err) < 1.0e-9, "fixed alpha error should be zero");
      }
      if (expected_coulomb_mode == "finite_source") {
        Expect(!finite_source_mode->empty(), "finite-source mode should be recorded");
        Expect(std::isfinite(*finite_source_radius_fm) && *finite_source_radius_fm > 0.0,
               "finite-source radius should be finite and positive");
      } else {
        Expect(finite_source_mode->empty(), "non-finite Coulomb modes should not record finiteSourceMode");
      }
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

  void ExpectCoulombKernelCatalog(const std::filesystem::path &root_path, const bool expect_rows) {
    TFile file(root_path.string().c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(file.Get("meta/CoulombKernelCatalog"));
    Expect(tree != nullptr, "CoulombKernelCatalog missing");
    Expect(tree->GetBranch("finite_source_mode") != nullptr, "kernel catalog finite_source_mode branch missing");
    Expect(tree->GetBranch("seed_radius_fm") != nullptr, "kernel catalog seed radius branch missing");
    Expect(tree->GetBranch("final_radius_fm") != nullptr, "kernel catalog final radius branch missing");
    Expect(tree->GetBranch("kstar_bin_count") != nullptr, "kernel catalog kstar bin branch missing");
    Expect(tree->GetBranch("qn_index") != nullptr, "kernel catalog qn_index branch missing");
    Expect(tree->GetBranch("qn_label") != nullptr, "kernel catalog qn_label branch missing");
    Expect(tree->GetBranch("is_qn_integrated") != nullptr, "kernel catalog is_qn_integrated branch missing");
    if (expect_rows) {
      Expect(tree->GetEntries() > 0, "finite-source kernel catalog should contain rows");
    }
  }

  void ExpectReportOutputs(const std::filesystem::path &fit_report_root_path) {
    TFile report_file(fit_report_root_path.string().c_str(), "READ");
    Expect(!report_file.IsZombie(), "fit report ROOT file should be readable");
    Expect(report_file.Get("meta/FitCatalog") != nullptr, "fit report should include FitCatalog summary fields");
    Expect(report_file.Get("meta/CoulombKernelCatalog") != nullptr,
           "fit report should include CoulombKernelCatalog summary fields");
    Expect(report_file.Get("summary/R2_vs_phi/cent_0.00-10.00__mt_0.20-0.40/Rside2_vs_phi") != nullptr,
           "fit report should include legacy R2 summary graphs");
    Expect(report_file.Get("source_parameters/cent_0.00-10.00/mt_0.20-0.40/source_parameters_overview_canvas")
               != nullptr,
           "fit report should include source parameter overview canvas");
    Expect(report_file.Get("eps_vs_mt/cent_0.00-10.00/epsf_vs_mt") != nullptr,
           "fit report should include eps vs mt graph");
    Expect(report_file.Get("eps_vs_mt/cent_0.00-10.00/epsf_vs_mt_canvas") != nullptr,
           "fit report should include eps vs mt canvas");
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
                                                            "mapped_follow_report.root",
                                                            true,
                                                            std::nullopt);
  const ApplicationConfig mapped_follow_config = LoadApplicationConfig(mapped_follow_config_path);
  const BuildCfRunStatistics mapped_build_stats = RunBuildCf(mapped_follow_config, logger);
  Expect(mapped_build_stats.stored_slices == 5, "mapped build-cf should produce 5 seam-safe slices");
  const FitRunStatistics mapped_follow_fit_stats = RunFit(mapped_follow_config, logger);
  Expect(mapped_follow_fit_stats.selected_slices == 5, "follow-input fit should select every built slice");

  TFile mapped_cf_file((temp_dir / "mapped_cf.root").string().c_str(), "READ");
  Expect(mapped_cf_file.Get("meta/SliceCatalog") != nullptr, "mapped SliceCatalog missing");
  Expect(mapped_cf_file.Get("slices") != nullptr, "mapped slices directory missing");

  const std::string report_collision_config_path = WriteConfig(temp_dir / "report_collision.toml",
                                                               input_root,
                                                               temp_dir.string(),
                                                               "mapped_cf.root",
                                                               "report_collision_fit.root",
                                                               "report_collision.tsv",
                                                               "report_collision_fit.root",
                                                               true,
                                                               std::nullopt);
  const ApplicationConfig report_collision_config = LoadApplicationConfig(report_collision_config_path);
  bool saw_report_collision = false;
  try {
    (void)RunFit(report_collision_config, logger);
  } catch (const std::runtime_error &error) {
    saw_report_collision = std::string(error.what()).find("distinct from the CF and detailed fit ROOT files")
                           != std::string::npos;
  }
  Expect(saw_report_collision, "fit report path must not overwrite the detailed fit ROOT file");

  const FitCatalogInspection mapped_follow_inspection = InspectFitCatalog(temp_dir / "mapped_follow_fit.root");
  Expect(mapped_follow_inspection.fit_uses_symmetric_phi_range,
         "fit should follow the mapped CF metadata when no fit override is given");
  Expect(mapped_follow_inspection.min_phi < 0.0, "mapped follow fit should keep a negative phi slice");
  const auto mapped_follow_phi_range = InspectSummaryPhiFitRange(temp_dir / "mapped_follow_fit.root");
  Expect(std::abs(mapped_follow_phi_range.first + TMath::Pi() / 2.0) < 1.0e-6,
         "mapped follow phi fit minimum should be -pi/2");
  Expect(std::abs(mapped_follow_phi_range.second - TMath::Pi() / 2.0) < 1.0e-6,
         "mapped follow phi fit maximum should be pi/2");
  ExpectReportOutputs(temp_dir / "mapped_follow_report.root");

  const std::string fixed_parameter_lines = R"toml(
[fit.parameters.lambda]
fixed_value = 0.65

[fit.parameters.alpha]
fixed_value = 1.20

)toml";

  const std::string fixed_non_pml_config_path = WriteConfig(temp_dir / "fixed_non_pml.toml",
                                                            input_root,
                                                            temp_dir.string(),
                                                            "mapped_cf.root",
                                                            "fixed_non_pml_fit.root",
                                                            "fixed_non_pml.tsv",
                                                            "fixed_non_pml_report.root",
                                                            true,
                                                            std::nullopt,
                                                            "use_coulomb = false\n",
                                                            false,
                                                            fixed_parameter_lines);
  const ApplicationConfig fixed_non_pml_config = LoadApplicationConfig(fixed_non_pml_config_path);
  const FitRunStatistics fixed_non_pml_stats = RunFit(fixed_non_pml_config, logger);
  Expect(fixed_non_pml_stats.fitted_slices == 5, "fixed-parameter chi2 fit should fit every selected slice");
  (void)InspectFitCatalog(temp_dir / "fixed_non_pml_fit.root", "none", 0.65, 1.20);

  const std::string fixed_pml_config_path = WriteConfig(temp_dir / "fixed_pml.toml",
                                                        input_root,
                                                        temp_dir.string(),
                                                        "mapped_cf.root",
                                                        "fixed_pml_fit.root",
                                                        "fixed_pml.tsv",
                                                        "fixed_pml_report.root",
                                                        true,
                                                        std::nullopt,
                                                        "use_coulomb = false\n",
                                                        true,
                                                        fixed_parameter_lines);
  const ApplicationConfig fixed_pml_config = LoadApplicationConfig(fixed_pml_config_path);
  const FitRunStatistics fixed_pml_stats = RunFit(fixed_pml_config, logger);
  Expect(fixed_pml_stats.fitted_slices == 5, "fixed-parameter PML fit should fit every selected slice");
  (void)InspectFitCatalog(temp_dir / "fixed_pml_fit.root", "none", 0.65, 1.20);

  const std::string mapped_override_raw_config_path = WriteConfig(temp_dir / "mapped_override_raw.toml",
                                                                  input_root,
                                                                  temp_dir.string(),
                                                                  "mapped_cf.root",
                                                                  "mapped_override_raw_fit.root",
                                                                  "mapped_override_raw.tsv",
                                                                  "mapped_override_raw_report.root",
                                                                  false,
                                                                  false);
  const ApplicationConfig mapped_override_raw_config = LoadApplicationConfig(mapped_override_raw_config_path);
  const FitRunStatistics mapped_override_raw_fit_stats = RunFit(mapped_override_raw_config, logger);
  Expect(mapped_override_raw_fit_stats.selected_slices == 5,
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
                                                                  "raw_override_mapped_report.root",
                                                                  false,
                                                                  true);
  const ApplicationConfig raw_override_mapped_config = LoadApplicationConfig(raw_override_mapped_config_path);
  const BuildCfRunStatistics raw_build_stats = RunBuildCf(raw_override_mapped_config, logger);
  Expect(raw_build_stats.stored_slices == 5, "raw build-cf should produce 5 slices");
  const FitRunStatistics raw_override_mapped_fit_stats = RunFit(raw_override_mapped_config, logger);
  Expect(raw_override_mapped_fit_stats.selected_slices == 5,
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
  Expect(header.find("coulombMode") != std::string::npos, "summary TSV coulombMode header missing");
  Expect(header.find("finiteSourceRadiusFm") != std::string::npos,
         "summary TSV finiteSourceRadiusFm header missing");
  ExpectCoulombKernelCatalog(temp_dir / "mapped_follow_fit.root", false);
  ExpectCoulombKernelCatalog(temp_dir / "mapped_follow_report.root", false);

#ifndef EXP_FEMTO_3D_HAS_CATS
  const std::string finite_unavailable_config_path = WriteConfig(temp_dir / "finite_unavailable.toml",
                                                                 input_root,
                                                                 temp_dir.string(),
                                                                 "mapped_cf.root",
                                                                 "finite_unavailable_fit.root",
                                                                 "finite_unavailable.tsv",
                                                                 "finite_unavailable_report.root",
                                                                 true,
                                                                 std::nullopt,
                                                                 "coulomb_mode = \"finite_source\"\n"
                                                                 "finite_source_mode = \"fixed_1d\"\n");
  const ApplicationConfig finite_unavailable_config = LoadApplicationConfig(finite_unavailable_config_path);
  bool saw_missing_cats_error = false;
  try {
    (void)RunFit(finite_unavailable_config, logger);
  } catch (const std::runtime_error &error) {
    saw_missing_cats_error = std::string(error.what()).find("requires CATS support") != std::string::npos;
  }
  Expect(saw_missing_cats_error, "finite_source should fail clearly when CATS is not compiled in");
#endif

#ifdef EXP_FEMTO_3D_HAS_CATS
  const std::string finite_fixed_config_path = WriteConfig(temp_dir / "finite_fixed.toml",
                                                           input_root,
                                                           temp_dir.string(),
                                                           "mapped_cf.root",
                                                           "finite_fixed_fit.root",
                                                           "finite_fixed.tsv",
                                                           "finite_fixed_report.root",
                                                           true,
                                                           std::nullopt,
                                                           "coulomb_mode = \"finite_source\"\n"
                                                           "finite_source_mode = \"fixed_1d\"\n");
  const ApplicationConfig finite_fixed_config = LoadApplicationConfig(finite_fixed_config_path);
  const FitRunStatistics finite_fixed_stats = RunFit(finite_fixed_config, logger);
  Expect(finite_fixed_stats.selected_slices == 5, "finite-source fit should select every built slice");
  Expect(finite_fixed_stats.fitted_slices == 5, "finite-source fit should fit every selected slice");
  (void)InspectFitCatalog(temp_dir / "finite_fixed_fit.root", "finite_source");
  ExpectCoulombKernelCatalog(temp_dir / "finite_fixed_fit.root", true);
  ExpectCoulombKernelCatalog(temp_dir / "finite_fixed_report.root", true);

  const std::string finite_iterative_config_path = WriteConfig(temp_dir / "finite_iterative.toml",
                                                               input_root,
                                                               temp_dir.string(),
                                                               "mapped_cf.root",
                                                               "finite_iterative_fit.root",
                                                               "finite_iterative.tsv",
                                                               "finite_iterative_report.root",
                                                               true,
                                                               std::nullopt,
                                                               "coulomb_mode = \"finite_source\"\n"
                                                               "finite_source_mode = \"iterative_1d\"\n");
  const ApplicationConfig finite_iterative_config = LoadApplicationConfig(finite_iterative_config_path);
  const FitRunStatistics finite_iterative_stats = RunFit(finite_iterative_config, logger);
  Expect(finite_iterative_stats.selected_slices == 5, "iterative finite-source fit should select every built slice");
  Expect(finite_iterative_stats.fitted_slices == 5, "iterative finite-source fit should fit every selected slice");
  (void)InspectFitCatalog(temp_dir / "finite_iterative_fit.root", "finite_source");
  ExpectCoulombKernelCatalog(temp_dir / "finite_iterative_fit.root", true);
  ExpectCoulombKernelCatalog(temp_dir / "finite_iterative_report.root", true);
#endif

  std::cout << "workflow_smoke_test passed\n";
  return 0;
}
