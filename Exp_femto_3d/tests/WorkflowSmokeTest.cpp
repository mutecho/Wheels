#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2.h"
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
    const double max[7] = {0.15, 0.15, 0.15, 0.4, 20.0, 0.5, 3.14159265358979323846};

    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, min, max);
    for (double phi : {0.3, 1.1, 2.0, 2.8}) {
      for (const double centrality : {5.0, 15.0}) {
        FillSparse(*same, 0.01, 0.00, 0.01, 0.3, centrality, phi, 50.0);
        FillSparse(*same, -0.01, 0.01, -0.01, 0.3, centrality, phi, 42.0);
        FillSparse(*same, 0.02, -0.01, 0.00, 0.3, centrality, phi, 35.0);
        FillSparse(*mixed, 0.01, 0.00, 0.01, 0.3, centrality, phi, 55.0);
        FillSparse(*mixed, -0.01, 0.01, -0.01, 0.3, centrality, phi, 48.0);
        FillSparse(*mixed, 0.02, -0.01, 0.00, 0.3, centrality, phi, 40.0);
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
                          const std::string &cf_root_name,
                          const std::string &fit_root_name,
                          const std::string &fit_summary_name,
                          const std::string &fit_report_root_name,
                          const bool build_map_pair_phi_to_symmetric_range,
                          const std::optional<bool> fit_map_pair_phi_to_symmetric_range,
                          const std::string &fit_coulomb_lines = "use_coulomb = false\n",
                          const bool use_pml = false,
                          const std::string &fit_parameter_lines = "",
                          const std::string &profile_likelihood_lines = "",
                          const std::string &profile_root_name = "profile_likelihood.root",
                          const bool fit_progress = false,
                          const std::string &additional_selection_bins = "") {
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
    output << "profile_root_name = \"" << profile_root_name << "\"\n";
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
    output << "progress = " << (fit_progress ? "true" : "false") << "\n\n";
    output << fit_parameter_lines;
    output << profile_likelihood_lines;
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.mt]]\nmin = 0.2\nmax = 0.4\n";
    output << additional_selection_bins;
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

  void ExpectProfileTreeContract(TTree &tree, const std::string &name) {
    Expect(tree.GetBranch("coordinates") != nullptr, name + " coordinates branch missing");
    Expect(tree.GetBranch("objective") != nullptr, name + " objective branch missing");
    Expect(tree.GetBranch("status") != nullptr, name + " status branch missing");
    Expect(tree.GetBranch("parameter_values") != nullptr, name + " parameter_values branch missing");
    Expect(tree.GetBranch("parameter_errors") != nullptr, name + " parameter_errors branch missing");
    Expect(tree.GetBranch("at_lower_bound") != nullptr, name + " lower-bound branch missing");
    Expect(tree.GetBranch("at_upper_bound") != nullptr, name + " upper-bound branch missing");
    Expect(tree.GetBranch("edm") != nullptr, name + " EDM branch missing");
    Expect(tree.GetBranch("model_domain_valid") != nullptr, name + " model-domain validity branch missing");
    Expect(tree.GetBranch("objective_valid") != nullptr, name + " objective-validity branch missing");
    Expect(tree.GetBranch("fcn_calls") != nullptr, name + " FCN count branch missing");
    Expect(tree.GetBranch("total_wall_ms") != nullptr, name + " total timing branch missing");
    Expect(tree.GetBranch("migrad_wall_ms") != nullptr, name + " MIGRAD timing branch missing");
    Expect(tree.GetBranch("hesse_wall_ms") != nullptr, name + " HESSE timing branch missing");
    Expect(tree.GetBranch("hesse_ran") != nullptr, name + " HESSE execution flag missing");
    Expect(tree.GetBranch("parameter_errors_valid") != nullptr, name + " parameter-error validity flag missing");
  }

  bool EquivalentProfileNumber(const double left, const double right) {
    if (std::isnan(left) && std::isnan(right)) return true;
    return std::abs(left - right) <= 1.0e-10 * (1.0 + std::max(std::abs(left), std::abs(right)));
  }

  void ExpectProfileNumericalEquivalence(const std::filesystem::path &serial_path,
                                         const std::filesystem::path &process_path) {
    TFile serial(serial_path.string().c_str(), "READ");
    TFile process(process_path.string().c_str(), "READ");
    auto *serial_catalog = dynamic_cast<TTree *>(serial.Get("meta/ProfileLikelihoodCatalog"));
    auto *process_catalog = dynamic_cast<TTree *>(process.Get("meta/ProfileLikelihoodCatalog"));
    Expect(!serial.IsZombie() && !process.IsZombie() && serial_catalog != nullptr && process_catalog != nullptr,
           "serial/process profile catalogs must be readable");
    Expect(serial_catalog->GetEntries() == process_catalog->GetEntries(),
           "serial/process profile catalogs must have the same rows");
    TTreeReader serial_catalog_reader(serial_catalog);
    TTreeReader process_catalog_reader(process_catalog);
    TTreeReaderValue<std::string> serial_slice(serial_catalog_reader, "slice_id");
    TTreeReaderValue<std::string> process_slice(process_catalog_reader, "slice_id");
    TTreeReaderValue<std::string> serial_scan(serial_catalog_reader, "scan_id");
    TTreeReaderValue<std::string> process_scan(process_catalog_reader, "scan_id");
    TTreeReaderValue<std::string> serial_reference_source(serial_catalog_reader, "reference_source");
    TTreeReaderValue<std::string> process_reference_source(process_catalog_reader, "reference_source");
    TTreeReaderValue<double> serial_reference(serial_catalog_reader, "reference_objective");
    TTreeReaderValue<double> process_reference(process_catalog_reader, "reference_objective");
    while (serial_catalog_reader.Next()) {
      Expect(process_catalog_reader.Next(), "process profile catalog row missing");
      Expect(*serial_slice == *process_slice && *serial_scan == *process_scan
                 && *serial_reference_source == *process_reference_source
                 && EquivalentProfileNumber(*serial_reference, *process_reference),
             "serial/process catalog reference mismatch");
      const std::string directory = "profiles/" + *serial_slice + "/" + *serial_scan + "/";
      auto *serial_points = dynamic_cast<TTree *>(serial.Get((directory + "ProfilePoints").c_str()));
      auto *process_points = dynamic_cast<TTree *>(process.Get((directory + "ProfilePoints").c_str()));
      auto *serial_attempts = dynamic_cast<TTree *>(serial.Get((directory + "AttemptPoints").c_str()));
      auto *process_attempts = dynamic_cast<TTree *>(process.Get((directory + "AttemptPoints").c_str()));
      Expect(serial_points != nullptr && process_points != nullptr && serial_attempts != nullptr && process_attempts != nullptr,
             "serial/process numerical profile trees missing");
      Expect(serial_points->GetEntries() == process_points->GetEntries()
                 && serial_attempts->GetEntries() == process_attempts->GetEntries(),
             "serial/process point or attempt row count mismatch");

      TTreeReader serial_point_reader(serial_points);
      TTreeReader process_point_reader(process_points);
      TTreeReaderValue<std::vector<double>> serial_coordinates(serial_point_reader, "coordinates");
      TTreeReaderValue<std::vector<double>> process_coordinates(process_point_reader, "coordinates");
      TTreeReaderValue<double> serial_objective(serial_point_reader, "objective");
      TTreeReaderValue<double> process_objective(process_point_reader, "objective");
      TTreeReaderValue<double> serial_delta(serial_point_reader, "delta_neg2logl");
      TTreeReaderValue<double> process_delta(process_point_reader, "delta_neg2logl");
      TTreeReaderValue<std::string> serial_status(serial_point_reader, "status");
      TTreeReaderValue<std::string> process_status(process_point_reader, "status");
      TTreeReaderValue<std::string> serial_winner(serial_point_reader, "winner_seed");
      TTreeReaderValue<std::string> process_winner(process_point_reader, "winner_seed");
      while (serial_point_reader.Next()) {
        Expect(process_point_reader.Next(), "process ProfilePoints row missing");
        Expect(serial_coordinates->size() == process_coordinates->size(),
               "serial/process profile coordinate dimensionality mismatch");
        for (std::size_t coordinate = 0; coordinate < serial_coordinates->size(); ++coordinate) {
          Expect(EquivalentProfileNumber(serial_coordinates->at(coordinate), process_coordinates->at(coordinate)),
                 "serial/process profile coordinate mismatch");
        }
        Expect(*serial_status == *process_status && *serial_winner == *process_winner
                   && EquivalentProfileNumber(*serial_objective, *process_objective)
                   && EquivalentProfileNumber(*serial_delta, *process_delta),
               "serial/process profile winner, status, objective, or delta mismatch");
      }

      TTreeReader serial_attempt_reader(serial_attempts);
      TTreeReader process_attempt_reader(process_attempts);
      TTreeReaderValue<int> serial_attempt_index(serial_attempt_reader, "attempt_index");
      TTreeReaderValue<int> process_attempt_index(process_attempt_reader, "attempt_index");
      TTreeReaderValue<int> serial_point_index(serial_attempt_reader, "point_index");
      TTreeReaderValue<int> process_point_index(process_attempt_reader, "point_index");
      TTreeReaderValue<std::string> serial_seed(serial_attempt_reader, "seed_origin");
      TTreeReaderValue<std::string> process_seed(process_attempt_reader, "seed_origin");
      TTreeReaderValue<std::string> serial_attempt_status(serial_attempt_reader, "status");
      TTreeReaderValue<std::string> process_attempt_status(process_attempt_reader, "status");
      while (serial_attempt_reader.Next()) {
        Expect(process_attempt_reader.Next(), "process AttemptPoints row missing");
        Expect(*serial_attempt_index == *process_attempt_index
                   && *serial_point_index == *process_point_index
                   && *serial_seed == *process_seed
                   && *serial_attempt_status == *process_attempt_status,
               "serial/process attempt order, seed, or status mismatch");
      }
    }
    Expect(!process_catalog_reader.Next(), "process profile catalog has extra rows");
  }

  void ExpectProfileOutput(const std::filesystem::path &profile_root_path,
                           const std::string &slice_id) {
    TFile file(profile_root_path.string().c_str(), "READ");
    Expect(!file.IsZombie(), "profile ROOT file must be readable");
    auto *catalog = dynamic_cast<TTree *>(file.Get("meta/ProfileLikelihoodCatalog"));
    auto *parameters = dynamic_cast<TTree *>(file.Get("meta/ProfileParameterCatalog"));
    Expect(catalog != nullptr, "ProfileLikelihoodCatalog missing");
    Expect(parameters != nullptr, "ProfileParameterCatalog missing");
    Expect(catalog->GetBranch("objective_kind") != nullptr, "profile catalog objective kind missing");
    Expect(catalog->GetBranch("finite_kernel_frozen") != nullptr, "profile catalog frozen-kernel metadata missing");
    Expect(catalog->GetBranch("reference_objective") != nullptr, "profile catalog reference objective missing");
    Expect(catalog->GetEntries() == 3, "profile catalog should contain the requested 1D, 2D, and PSD scans");
    Expect(parameters->GetEntries() == 10, "full-model parameter catalog should contain all ten physical parameters");

    TTreeReader catalog_reader(catalog);
    TTreeReaderValue<std::string> objective_kind(catalog_reader, "objective_kind");
    while (catalog_reader.Next()) {
      Expect(*objective_kind == "neg2logl_pml", "profile objective kind must identify the PML -2 ln L statistic");
    }

    const std::string base = "profiles/" + slice_id + "/";
    auto *one_d_points = dynamic_cast<TTree *>(file.Get((base + "rout2/ProfilePoints").c_str()));
    auto *one_d_attempts = dynamic_cast<TTree *>(file.Get((base + "rout2/AttemptPoints").c_str()));
    Expect(one_d_points != nullptr && one_d_attempts != nullptr, "1D profile point/attempt trees missing");
    ExpectProfileTreeContract(*one_d_points, "1D ProfilePoints");
    ExpectProfileTreeContract(*one_d_attempts, "1D AttemptPoints");
    Expect(one_d_points->GetEntries() == 3, "1D ProfilePoints must retain every requested coordinate exactly once");
    Expect(one_d_attempts->GetEntries() >= one_d_points->GetEntries(), "1D attempts must retain nominal/retry records");
    TTreeReader one_d_reader(one_d_points);
    TTreeReaderValue<int> one_d_stage(one_d_reader, "stage");
    TTreeReaderValue<std::vector<double>> one_d_coordinates(one_d_reader, "coordinates");
    TTreeReaderValue<std::string> one_d_status(one_d_reader, "status");
    TTreeReaderValue<double> one_d_objective(one_d_reader, "objective");
    TTreeReaderValue<double> one_d_slice_objective(one_d_reader, "slice_objective");
    int one_d_index = 0;
    int valid_one_d_points = 0;
    for (const double expected : {0.01, 1.005, 2.0}) {
      Expect(one_d_reader.Next(), "1D ProfilePoints row missing");
      Expect(*one_d_stage == 0 && one_d_coordinates->size() == 1U,
             "1D ProfilePoints must preserve coarse stage and axis dimensionality");
      Expect(std::abs(one_d_coordinates->at(0) - expected) < 1.0e-12,
             "1D ProfilePoints must preserve exact configured grid coordinates");
      if (*one_d_status == "valid") {
        ++valid_one_d_points;
        const double tolerance = 1.0e-7 * (1.0 + std::abs(*one_d_slice_objective));
        Expect(std::isfinite(*one_d_objective) && std::isfinite(*one_d_slice_objective)
                   && *one_d_objective <= *one_d_slice_objective + tolerance,
               "valid profile objective must not exceed its fixed-nuisance likelihood slice");
      }
      ++one_d_index;
    }
    Expect(!one_d_reader.Next() && one_d_index == 3, "1D ProfilePoints must have no implicit under/overflow rows");
    Expect(valid_one_d_points > 0, "toy profile must retain at least one valid profiled point for slice comparison");
    Expect(dynamic_cast<TGraph *>(file.Get((base + "rout2/Profile1D").c_str())) != nullptr,
           "1D profiled likelihood graph missing");
    Expect(dynamic_cast<TGraph *>(file.Get((base + "rout2/Slice1D").c_str())) != nullptr,
           "1D fixed-nuisance likelihood slice graph missing");
    Expect(dynamic_cast<TGraph *>(file.Get((base + "rout2/Nuisance_norm").c_str())) != nullptr,
           "named 1D nuisance trajectory graph missing");
    for (int parameter = 0; parameter < 16; ++parameter) {
      Expect(file.Get((base + "rout2/Nuisance_p" + std::to_string(parameter)).c_str()) == nullptr,
             "positional nuisance object must not be written");
    }
    Expect(file.Get((base + "rout2/Canvas_1D").c_str()) != nullptr, "1D profile canvas missing");
    Expect(file.Get((base + "rout2/NominalPoint").c_str()) != nullptr, "nominal marker missing");
    for (const std::string &removed : {"Profile1D_Coarse",
                                       "Profile1D_Refined",
                                       "FailurePoints1D",
                                       "ProfileMinimum",
                                       "LowerBound",
                                       "UpperBound",
                                       "Canvas_Nuisance"}) {
      Expect(file.Get((base + "rout2/" + removed).c_str()) == nullptr,
             "redundant 1D QA object must not be written: " + removed);
    }

    auto *two_d_points = dynamic_cast<TTree *>(file.Get((base + "rout2_rside2/ProfilePoints").c_str()));
    auto *two_d_attempts = dynamic_cast<TTree *>(file.Get((base + "rout2_rside2/AttemptPoints").c_str()));
    auto *delta = dynamic_cast<TH2 *>(file.Get((base + "rout2_rside2/DeltaNeg2LogL2D").c_str()));
    auto *status = dynamic_cast<TH2 *>(file.Get((base + "rout2_rside2/PointStatus2D").c_str()));
    Expect(two_d_points != nullptr && two_d_attempts != nullptr, "2D profile point/attempt trees missing");
    Expect(delta != nullptr && status != nullptr, "2D profile heatmap/status mask missing");
    ExpectProfileTreeContract(*two_d_points, "2D ProfilePoints");
    ExpectProfileTreeContract(*two_d_attempts, "2D AttemptPoints");
    Expect(two_d_points->GetEntries() == 9, "2D ProfilePoints must retain every requested coordinate exactly once");
    Expect(two_d_attempts->GetEntries() >= two_d_points->GetEntries(), "2D attempts must retain every requested point");
    for (int xbin = 0; xbin <= delta->GetNbinsX() + 1; ++xbin) {
      Expect(delta->GetBinContent(xbin, 0) == 0.0
                 && delta->GetBinContent(xbin, delta->GetNbinsY() + 1) == 0.0,
             "2D likelihood heatmap must leave Y underflow/overflow empty");
    }
    for (int ybin = 0; ybin <= delta->GetNbinsY() + 1; ++ybin) {
      Expect(delta->GetBinContent(0, ybin) == 0.0
                 && delta->GetBinContent(delta->GetNbinsX() + 1, ybin) == 0.0,
             "2D likelihood heatmap must leave X underflow/overflow empty");
    }
    Expect(file.Get((base + "rout2_rside2/Canvas_2D").c_str()) != nullptr, "2D diagnostic canvas missing");
    Expect(file.Get((base + "rout2_rside2/NominalPoint").c_str()) != nullptr, "2D nominal marker missing");
    for (const std::string &removed : {"FailurePoints2D",
                                       "RefinedPoints2D",
                                       "ProfileMinimum",
                                       "LowerBound",
                                       "UpperBound"}) {
      Expect(file.Get((base + "rout2_rside2/" + removed).c_str()) == nullptr,
             "redundant 2D QA object must not be written: " + removed);
    }

    auto *psd_points = dynamic_cast<TTree *>(file.Get((base + "psd_invalid/ProfilePoints").c_str()));
    Expect(psd_points != nullptr, "full-model PSD-invalid profile tree missing");
    TTreeReader psd_reader(psd_points);
    TTreeReaderValue<std::string> psd_status(psd_reader, "status");
    TTreeReaderValue<int> psd_model_domain_valid(psd_reader, "model_domain_valid");
    TTreeReaderValue<double> psd_delta(psd_reader, "delta_neg2logl");
    int invalid_psd_points = 0;
    while (psd_reader.Next()) {
      if (*psd_status == "model_domain_invalid") {
        ++invalid_psd_points;
        Expect(*psd_model_domain_valid == 0, "PSD-invalid point must carry a false model-domain flag");
        Expect(!std::isfinite(*psd_delta), "PSD-invalid point must not be displayed as a valid likelihood delta");
      }
    }
    Expect(invalid_psd_points >= 1, "an unreachable full-model cross term must classify as model_domain_invalid");
    Expect(file.Get((base + "psd_invalid/Failure_model_domain_invalid1D").c_str()) == nullptr,
           "failure category is numerical tree data and must not be duplicated as a QA graph");
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
                                                        fixed_parameter_lines,
                                                        "",
                                                        "pml_profile_disabled.root");
  std::filesystem::remove(temp_dir / "pml_profile_disabled.root");
  const ApplicationConfig fixed_pml_config = LoadApplicationConfig(fixed_pml_config_path);
  const FitRunStatistics fixed_pml_stats = RunFit(fixed_pml_config, logger);
  Expect(fixed_pml_stats.fitted_slices == 5, "fixed-parameter PML fit should fit every selected slice");
  (void)InspectFitCatalog(temp_dir / "fixed_pml_fit.root", "none", 0.65, 1.20);
  Expect(!std::filesystem::exists(temp_dir / "pml_profile_disabled.root"),
         "disabled profile mode must not create or reset its separate ROOT output");

  const std::vector<SliceCatalogEntry> profile_catalog = LoadSliceCatalog((temp_dir / "mapped_cf.root").string());
  Expect(!profile_catalog.empty(), "toy CF must provide a slice for profile diagnostics");
  const std::string profile_slice_id = profile_catalog.front().slice_id;
  const std::string profile_likelihood_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]
retry_strategy = "reference_only"
write_likelihood_slice = true
contour_levels = [1.0, 2.0, 4.0]

[[fit.profile_likelihood.scans]]
id = "rout2"
parameters = ["rout2"]
points = [3]
min = [0.01]
max = [2.0]
refine = false

[[fit.profile_likelihood.scans]]
id = "rout2_rside2"
parameters = ["rout2", "rside2"]
points = [3, 3]
min = [0.01, 0.01]
max = [2.0, 2.0]
refine = false

[[fit.profile_likelihood.scans]]
id = "psd_invalid"
parameters = ["routside2"]
points = [3]
min = [-1000.0]
max = [1000.0]
refine = false

)toml";
  const std::string profile_config_path = WriteConfig(temp_dir / "profile_enabled.toml",
                                                       input_root,
                                                       temp_dir.string(),
                                                       "mapped_cf.root",
                                                       "profile_enabled_fit.root",
                                                       "profile_enabled.tsv",
                                                       "profile_enabled_report.root",
                                                       true,
                                                       std::nullopt,
                                                       "use_coulomb = false\n",
                                                       true,
                                                       fixed_parameter_lines,
                                                       profile_likelihood_lines,
                                                       "profile_enabled.root");
  const ApplicationConfig profile_config = LoadApplicationConfig(profile_config_path);
  const FitRunStatistics profile_stats = RunFit(profile_config, logger);
  Expect(profile_stats.profile_selected_slices == 1, "listed profile mode must select only its exact slice ID");
  Expect(profile_stats.profile_completed_slices == 1, "listed toy profile must complete its nominated slice");
  ExpectProfileOutput(temp_dir / "profile_enabled.root", profile_slice_id);

  const std::string process_profile_lines = R"toml(
[fit.profile_likelihood]
enabled = true
execution_mode = "profile_only"
parallel_backend = "process"
workers = 2
minimizer_backend = "legacy_tminuit"
hesse_strategy = "none"
slice_scope = "fit_selection"
retry_strategy = "reference_only"
write_likelihood_slice = false

[fit.profile_likelihood.checkpoint]
enabled = true
resume = true
run_id = "toy_process_v1"
directory = "toy_process.work"

[[fit.profile_likelihood.scans]]
id = "rout2"
parameters = ["rout2"]
points = [3]
min = [0.01]
max = [2.0]
)toml";
  const std::string two_group_bins = R"toml(

[[bins.centrality]]
min = 10
max = 20
)toml";
  const std::string process_config_path = WriteConfig(
      temp_dir / "profile_process.toml", input_root, temp_dir.string(), "two_group_cf.root",
      "profile_process_production.root", "profile_process_production.tsv", "profile_process_report.root",
      true, std::nullopt, "use_coulomb = false\n", true, fixed_parameter_lines, process_profile_lines,
      "profile_process.root", true, two_group_bins);
  const ApplicationConfig process_config = LoadApplicationConfig(process_config_path);
  const BuildCfRunStatistics two_group_build_stats = RunBuildCf(process_config, logger);
  Expect(two_group_build_stats.stored_slices == 10,
         "two-group build must provide every phi-integrated and differential toy slice");
  std::filesystem::remove(temp_dir / "profile_process.root");
  std::filesystem::remove_all(temp_dir / "toy_process.work" / "toy_process_v1");
  const FitRunStatistics estimate_stats =
      RunFit(process_config, logger, std::nullopt, std::nullopt, true, process_config_path);
  Expect(estimate_stats.profile_estimate_only && estimate_stats.profile_estimated_attempts == 30,
         "estimate-only must report the exact reference-only coarse upper bound");
  Expect(estimate_stats.profile_selected_slices == 10
             && estimate_stats.profile_estimated_groups == 2
             && estimate_stats.profile_estimated_coarse_points_per_slice == 3
             && estimate_stats.profile_estimated_refined_points_per_slice == 0
             && estimate_stats.profile_effective_workers == 2,
         "estimate-only must report group, grid, refinement, and effective concurrency bounds");
  Expect(!std::filesystem::exists(temp_dir / "profile_process.root"),
         "estimate-only must not create the profile output");
  {
    std::ofstream fit_sentinel(temp_dir / "profile_process_production.root");
    fit_sentinel << "fit sentinel";
    std::ofstream report_sentinel(temp_dir / "profile_process_report.root");
    report_sentinel << "report sentinel";
    std::ofstream tsv_sentinel(temp_dir / "profile_process_production.tsv");
    tsv_sentinel << "tsv sentinel";
  }
  const FitRunStatistics process_stats =
      RunFit(process_config, logger, std::nullopt, std::nullopt, false, process_config_path);
  Expect(process_stats.profile_completed_slices == 10,
         "process profile-only run must complete every fit_selection slice");
  Expect(process_stats.profile_output_path == (temp_dir / "profile_process.root").string(),
         "profile-only statistics must expose the resolved diagnostic output path");
  {
    for (const std::string &group : {"cent0_mt0_qn-1", "cent1_mt0_qn-1"}) {
      std::ifstream worker_log(temp_dir / "toy_process.work" / "toy_process_v1" / (group + ".worker.log"));
      Expect(worker_log.is_open(), "process worker log must be retained for every assigned group");
      const std::string worker_output((std::istreambuf_iterator<char>(worker_log)), std::istreambuf_iterator<char>());
      Expect(worker_output.find("fit [") == std::string::npos,
             "process workers must not render independent terminal progress lines");
    }
  }
  {
    TFile process_output((temp_dir / "profile_process.root").string().c_str(), "READ");
    Expect(!process_output.IsZombie(), "merged process profile output must be readable");
    auto *execution = dynamic_cast<TTree *>(process_output.Get("meta/ProfileExecution"));
    auto *profile_catalog_tree = dynamic_cast<TTree *>(process_output.Get("meta/ProfileLikelihoodCatalog"));
    auto *attempts = dynamic_cast<TTree *>(
        process_output.Get(("profiles/" + profile_slice_id + "/rout2/AttemptPoints").c_str()));
    Expect(execution != nullptr && execution->GetEntries() == 1, "process execution metadata missing");
    Expect(profile_catalog_tree != nullptr && profile_catalog_tree->GetEntries() == 10,
           "merged fit_selection catalog must contain exactly selected slices times scans");
    std::string *slice_scope = nullptr;
    int selected_slices = 0;
    int selected_groups = 0;
    int configured_workers = 0;
    int effective_workers = 0;
    execution->SetBranchAddress("slice_scope", &slice_scope);
    execution->SetBranchAddress("selected_slices", &selected_slices);
    execution->SetBranchAddress("selected_groups", &selected_groups);
    execution->SetBranchAddress("configured_workers", &configured_workers);
    execution->SetBranchAddress("effective_workers", &effective_workers);
    execution->GetEntry(0);
    Expect(slice_scope != nullptr && *slice_scope == "fit_selection"
               && selected_slices == 10 && selected_groups == 2
               && configured_workers == 2 && effective_workers == 2,
           "process metadata must preserve parent scope and configured/effective concurrency");
    Expect(attempts != nullptr && attempts->GetBranch("fcn_calls") != nullptr
               && attempts->GetBranch("hesse_ran") != nullptr
               && attempts->GetBranch("parameter_errors_valid") != nullptr,
           "process attempt timing/HESSE metadata missing");
  }
  for (const auto &sentinel : {std::make_pair("profile_process_production.root", "fit sentinel"),
                               std::make_pair("profile_process_report.root", "report sentinel"),
                               std::make_pair("profile_process_production.tsv", "tsv sentinel")}) {
    std::ifstream input(temp_dir / sentinel.first);
    std::string contents;
    std::getline(input, contents);
    Expect(contents == sentinel.second,
           std::string("profile_only must not mutate production output ") + sentinel.first);
  }
  for (const auto &[group, expected_centrality] :
       {std::pair<std::string, int>{"cent0_mt0_qn-1", 0},
        std::pair<std::string, int>{"cent1_mt0_qn-1", 1}}) {
    TFile chunk((temp_dir / "toy_process.work" / "toy_process_v1" / (group + ".root")).string().c_str(), "READ");
    auto *chunk_catalog = dynamic_cast<TTree *>(chunk.Get("meta/ProfileLikelihoodCatalog"));
    Expect(!chunk.IsZombie() && chunk_catalog != nullptr && chunk_catalog->GetEntries() == 5,
           "each process chunk must contain only the five slices in its assigned group");
    int centrality_index = -1;
    chunk_catalog->SetBranchAddress("centrality_index", &centrality_index);
    for (Long64_t row = 0; row < chunk_catalog->GetEntries(); ++row) {
      chunk_catalog->GetEntry(row);
      Expect(centrality_index == expected_centrality,
             "worker chunk must not re-expand fit_selection outside its assigned group");
    }
  }
  const std::string serial_profile_lines = R"toml(
[fit.profile_likelihood]
enabled = true
execution_mode = "profile_only"
parallel_backend = "serial"
workers = 1
minimizer_backend = "legacy_tminuit"
hesse_strategy = "none"
slice_scope = "fit_selection"
retry_strategy = "reference_only"
write_likelihood_slice = false

[[fit.profile_likelihood.scans]]
id = "rout2"
parameters = ["rout2"]
points = [3]
min = [0.01]
max = [2.0]
)toml";
  const std::string serial_profile_config_path = WriteConfig(
      temp_dir / "profile_serial_fit_selection.toml", input_root, temp_dir.string(), "two_group_cf.root",
      "profile_serial_production.root", "profile_serial_production.tsv", "profile_serial_report.root",
      true, std::nullopt, "use_coulomb = false\n", true, fixed_parameter_lines, serial_profile_lines,
      "profile_serial.root", false, two_group_bins);
  const ApplicationConfig serial_profile_config = LoadApplicationConfig(serial_profile_config_path);
  const FitRunStatistics serial_profile_stats = RunFit(serial_profile_config, logger);
  Expect(serial_profile_stats.profile_completed_slices == 10,
         "serial profile-only comparison must complete the same fit_selection slices");
  ExpectProfileNumericalEquivalence(temp_dir / "profile_serial.root", temp_dir / "profile_process.root");
  {
    TFile checkpoint(
        (temp_dir / "toy_process.work" / "toy_process_v1" / "cent0_mt0_qn-1.root").string().c_str(),
        "UPDATE");
    Expect(!checkpoint.IsZombie(), "process checkpoint must be readable before resume pruning test");
    auto *scan_directory = checkpoint.GetDirectory(("profiles/" + profile_slice_id + "/rout2").c_str());
    Expect(scan_directory != nullptr, "process checkpoint scan directory missing");
    scan_directory->cd();
    TGraph legacy_alias(1);
    legacy_alias.SetPoint(0, 0.01, 1.0);
    legacy_alias.Write("Nuisance_p0", TObject::kOverwrite);
    checkpoint.Close();
  }
  const FitRunStatistics resumed_process_stats =
      RunFit(process_config, logger, std::nullopt, std::nullopt, false, process_config_path);
  Expect(resumed_process_stats.profile_completed_slices == 10,
         "matching process checkpoint must resume and still produce a complete result");
  {
    TFile resumed_output((temp_dir / "profile_process.root").string().c_str(), "READ");
    Expect(!resumed_output.IsZombie(), "resumed process output must be readable");
    Expect(resumed_output.Get(("profiles/" + profile_slice_id + "/rout2/Nuisance_p0").c_str()) == nullptr,
           "final merge must prune duplicate nuisance aliases from reused checkpoint chunks");
  }
  const auto complete_profile_size = std::filesystem::file_size(temp_dir / "profile_process.root");

  const std::string one_group_selection = two_group_bins + R"toml(

[[fit_selection.centrality]]
min = 0
max = 10

[[fit_selection.mt]]
min = 0.2
max = 0.4
)toml";
  const std::string changed_selection_config_path = WriteConfig(
      temp_dir / "profile_process_changed_selection.toml", input_root, temp_dir.string(), "two_group_cf.root",
      "profile_process_production.root", "profile_process_production.tsv", "profile_process_report.root",
      true, std::nullopt, "use_coulomb = false\n", true, fixed_parameter_lines, process_profile_lines,
      "profile_process.root", true, one_group_selection);
  const ApplicationConfig changed_selection_config = LoadApplicationConfig(changed_selection_config_path);
  bool saw_selection_checkpoint_mismatch = false;
  try {
    (void)RunFit(changed_selection_config, logger, std::nullopt, std::nullopt, false,
                 changed_selection_config_path);
  } catch (const std::runtime_error &error) {
    saw_selection_checkpoint_mismatch =
        std::string(error.what()).find("Checkpoint contract mismatch") != std::string::npos;
  }
  Expect(saw_selection_checkpoint_mismatch,
         "changing fit_selection must reject chunks built from the previous expanded slice list");
  Expect(std::filesystem::file_size(temp_dir / "profile_process.root") == complete_profile_size,
         "fit_selection checkpoint rejection must preserve the previous complete profile output");

  {
    TFile incomplete_chunk(
        (temp_dir / "toy_process.work" / "toy_process_v1" / "cent0_mt0_qn-1.root").string().c_str(),
        "UPDATE");
    auto *scan_directory = incomplete_chunk.GetDirectory(("profiles/" + profile_slice_id + "/rout2").c_str());
    Expect(scan_directory != nullptr, "checkpoint scan directory missing before completeness test");
    scan_directory->Delete("ProfilePoints;*");
    incomplete_chunk.Close();
  }
  bool saw_checkpoint_mismatch = false;
  try {
    (void)RunFit(process_config, logger, std::nullopt, std::nullopt, false, process_config_path);
  } catch (const std::runtime_error &error) {
    saw_checkpoint_mismatch = std::string(error.what()).find("Checkpoint contract mismatch") != std::string::npos;
  }
  Expect(saw_checkpoint_mismatch, "resume must reject a chunk missing ProfilePoints");
  Expect(std::filesystem::file_size(temp_dir / "profile_process.root") == complete_profile_size,
         "checkpoint rejection must preserve the previous complete profile output");

  std::filesystem::remove(temp_dir / "toy_process.work" / "toy_process_v1" / "cent0_mt0_qn-1.root");
  std::filesystem::remove(temp_dir / "toy_process.work" / "toy_process_v1" / "cent0_mt0_qn-1.digest");
  (void)RunFit(process_config, logger, std::nullopt, std::nullopt, false, process_config_path);
  const auto repaired_profile_size = std::filesystem::file_size(temp_dir / "profile_process.root");
  {
    TFile incomplete_chunk(
        (temp_dir / "toy_process.work" / "toy_process_v1" / "cent0_mt0_qn-1.root").string().c_str(),
        "UPDATE");
    auto *scan_directory = incomplete_chunk.GetDirectory(("profiles/" + profile_slice_id + "/rout2").c_str());
    Expect(scan_directory != nullptr, "checkpoint scan directory missing before attempt-tree test");
    scan_directory->Delete("AttemptPoints;*");
    incomplete_chunk.Close();
  }
  bool saw_missing_attempts = false;
  try {
    (void)RunFit(process_config, logger, std::nullopt, std::nullopt, false, process_config_path);
  } catch (const std::runtime_error &error) {
    saw_missing_attempts = std::string(error.what()).find("Checkpoint contract mismatch") != std::string::npos;
  }
  Expect(saw_missing_attempts, "resume must reject a chunk missing AttemptPoints");
  Expect(std::filesystem::file_size(temp_dir / "profile_process.root") == repaired_profile_size,
         "missing-attempt checkpoint rejection must preserve the previous complete profile output");

  const std::string empty_selection_bins = R"toml(

[[bins.centrality]]
min = 10
max = 20

[[fit_selection.centrality]]
min = 10
max = 20

[[fit_selection.mt]]
min = 0.2
max = 0.4
)toml";
  const std::string empty_selection_config_path = WriteConfig(
      temp_dir / "profile_process_empty_selection.toml", input_root, temp_dir.string(), "mapped_cf.root",
      "empty_selection_fit.root", "empty_selection.tsv", "empty_selection_report.root",
      true, std::nullopt, "use_coulomb = false\n", true, fixed_parameter_lines, process_profile_lines,
      "empty_selection_profile.root", false, empty_selection_bins);
  std::filesystem::remove(temp_dir / "empty_selection_profile.root");
  const ApplicationConfig empty_selection_config = LoadApplicationConfig(empty_selection_config_path);
  bool saw_empty_selection = false;
  try {
    (void)RunFit(empty_selection_config, logger, std::nullopt, std::nullopt, false,
                 empty_selection_config_path);
  } catch (const std::runtime_error &error) {
    saw_empty_selection = std::string(error.what()).find("selected zero slices") != std::string::npos;
  }
  Expect(saw_empty_selection && !std::filesystem::exists(temp_dir / "empty_selection_profile.root"),
         "an empty materialized fit_selection must fail before creating profile output");

  const std::string dynamic_override_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]

[[fit.profile_likelihood.scans]]
id = "full_only_cross_term"
parameters = ["routside2"]
points = [3]
min = [-1.0]
max = [1.0]

)toml";
  const std::string dynamic_override_config_path = WriteConfig(temp_dir / "profile_dynamic_override.toml",
                                                                input_root,
                                                                temp_dir.string(),
                                                                "mapped_cf.root",
                                                                "profile_dynamic_override_fit.root",
                                                                "profile_dynamic_override.tsv",
                                                                "profile_dynamic_override_report.root",
                                                                true,
                                                                std::nullopt,
                                                                "use_coulomb = false\n",
                                                                true,
                                                                fixed_parameter_lines,
                                                                dynamic_override_lines,
                                                                "profile_dynamic_override.root");
  const ApplicationConfig dynamic_override_config = LoadApplicationConfig(dynamic_override_config_path);
  {
    std::ofstream sentinel(temp_dir / "profile_dynamic_override.root");
    sentinel << "profile sentinel";
  }
  bool saw_dynamic_target_error = false;
  try {
    (void)RunFit(dynamic_override_config, logger, FitModel::kDiag);
  } catch (const std::runtime_error &error) {
    saw_dynamic_target_error = std::string(error.what()).find("unknown, inactive, or fixed") != std::string::npos;
  }
  Expect(saw_dynamic_target_error, "CLI model override must revalidate inactive full-only profile targets");
  std::ifstream sentinel_input(temp_dir / "profile_dynamic_override.root");
  std::string sentinel_contents;
  std::getline(sentinel_input, sentinel_contents);
  Expect(sentinel_contents == "profile sentinel", "dynamic profile validation must precede profile output reset");

  const std::string unbounded_target_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]

[[fit.profile_likelihood.scans]]
id = "unbounded_cross_term"
parameters = ["routside2"]
points = [3]

)toml";
  const ApplicationConfig unbounded_target_config = LoadApplicationConfig(WriteConfig(
      temp_dir / "profile_unbounded_target.toml",
      input_root,
      temp_dir.string(),
      "mapped_cf.root",
      "profile_unbounded_target_fit.root",
      "profile_unbounded_target.tsv",
      "profile_unbounded_target_report.root",
      true,
      std::nullopt,
      "use_coulomb = false\n",
      true,
      fixed_parameter_lines,
      unbounded_target_lines,
      "profile_unbounded_target.root"));
  bool saw_unbounded_target_error = false;
  try {
    (void)RunFit(unbounded_target_config, logger);
  } catch (const std::runtime_error &error) {
    saw_unbounded_target_error = std::string(error.what()).find("no finite effective bounds") != std::string::npos;
  }
  Expect(saw_unbounded_target_error,
         "default-unbounded full-model profile targets must require explicit scan bounds");

  const std::string out_of_bounds_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]
[[fit.profile_likelihood.scans]]
id = "rout2_out_of_bounds"
parameters = ["rout2"]
points = [3]
min = [0.01]
max = [401.0]
)toml";
  const ApplicationConfig out_of_bounds_config = LoadApplicationConfig(WriteConfig(
      temp_dir / "profile_out_of_bounds.toml", input_root, temp_dir.string(), "mapped_cf.root",
      "profile_out_of_bounds_fit.root", "profile_out_of_bounds.tsv", "profile_out_of_bounds_report.root",
      true, std::nullopt, "use_coulomb = false\n", true, fixed_parameter_lines, out_of_bounds_lines,
      "profile_out_of_bounds.root"));
  bool saw_out_of_bounds_error = false;
  try {
    (void)RunFit(out_of_bounds_config, logger);
  } catch (const std::runtime_error &error) {
    saw_out_of_bounds_error = std::string(error.what()).find("outside its effective fit bounds") != std::string::npos;
  }
  Expect(saw_out_of_bounds_error, "bounded profile targets must reject scan ranges outside the fit limits");

  const std::string fixed_target_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]

[[fit.profile_likelihood.scans]]
id = "fixed_lambda"
parameters = ["lambda"]
points = [3]
min = [0.2]
max = [0.9]

)toml";
  const ApplicationConfig fixed_target_config = LoadApplicationConfig(WriteConfig(temp_dir / "profile_fixed_target.toml",
                                                                                    input_root,
                                                                                    temp_dir.string(),
                                                                                    "mapped_cf.root",
                                                                                    "profile_fixed_target_fit.root",
                                                                                    "profile_fixed_target.tsv",
                                                                                    "profile_fixed_target_report.root",
                                                                                    true,
                                                                                    std::nullopt,
                                                                                    "use_coulomb = false\n",
                                                                                    true,
                                                                                    fixed_parameter_lines,
                                                                                    fixed_target_lines,
                                                                                    "profile_fixed_target.root"));
  bool saw_fixed_target_error = false;
  try {
    (void)RunFit(fixed_target_config, logger);
  } catch (const std::runtime_error &error) {
    saw_fixed_target_error = std::string(error.what()).find("unknown, inactive, or fixed") != std::string::npos;
  }
  Expect(saw_fixed_target_error, "profile target validation must reject user-fixed parameters");

  const std::string unknown_target_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]

[[fit.profile_likelihood.scans]]
id = "unknown_target"
parameters = ["not_a_parameter"]
points = [3]
min = [0.0]
max = [1.0]

)toml";
  const ApplicationConfig unknown_target_config = LoadApplicationConfig(WriteConfig(temp_dir / "profile_unknown_target.toml",
                                                                                      input_root,
                                                                                      temp_dir.string(),
                                                                                      "mapped_cf.root",
                                                                                      "profile_unknown_target_fit.root",
                                                                                      "profile_unknown_target.tsv",
                                                                                      "profile_unknown_target_report.root",
                                                                                      true,
                                                                                      std::nullopt,
                                                                                      "use_coulomb = false\n",
                                                                                      true,
                                                                                      fixed_parameter_lines,
                                                                                      unknown_target_lines,
                                                                                      "profile_unknown_target.root"));
  bool saw_unknown_target_error = false;
  try {
    (void)RunFit(unknown_target_config, logger);
  } catch (const std::runtime_error &error) {
    saw_unknown_target_error = std::string(error.what()).find("unknown, inactive, or fixed") != std::string::npos;
  }
  Expect(saw_unknown_target_error, "profile target validation must reject unknown parameter names");

  const std::string collision_lines = R"toml(
[fit.profile_likelihood]
enabled = true
slice_scope = "listed"
slice_ids = [")toml" + profile_slice_id + R"toml("]
[[fit.profile_likelihood.scans]]
id = "rout2_collision"
parameters = ["rout2"]
points = [3]
min = [0.01]
max = [2.0]

)toml";
  const ApplicationConfig collision_config = LoadApplicationConfig(WriteConfig(temp_dir / "profile_collision.toml",
                                                                                 input_root,
                                                                                 temp_dir.string(),
                                                                                 "mapped_cf.root",
                                                                                 "profile_collision_fit.root",
                                                                                 "profile_collision.tsv",
                                                                                 "profile_collision_report.root",
                                                                                 true,
                                                                                 std::nullopt,
                                                                                 "use_coulomb = false\n",
                                                                                 true,
                                                                                 fixed_parameter_lines,
                                                                                 collision_lines,
                                                                                 "mapped_cf.root"));
  bool saw_profile_collision = false;
  try {
    (void)RunFit(collision_config, logger);
  } catch (const std::runtime_error &error) {
    saw_profile_collision = std::string(error.what()).find("profile_root_name") != std::string::npos;
  }
  Expect(saw_profile_collision, "profile ROOT output must not collide with the CF ROOT file");

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
