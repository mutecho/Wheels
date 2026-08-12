#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "TFile.h"
#include "TNamed.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TTree.h"
#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Logging.h"
#include "exp_femto_3d/Workflow.h"

namespace {

  void Expect(const bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void ExpectNear(const double actual, const double expected, const std::string &message) {
    if (std::abs(actual - expected) > 1.0e-9) {
      throw std::runtime_error(message + ": expected " + std::to_string(expected) + ", got " + std::to_string(actual));
    }
  }

  void FillSparse(THnSparseF &sparse,
                  const double mt,
                  const double centrality,
                  const double phi,
                  const double weight,
                  const double qn = 0.0) {
    double values[7] = {0.01, -0.01, 0.02, mt, centrality, qn, phi};
    sparse.Fill(values, weight);
  }

  void WriteSentinelRootFile(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    Expect(!output.IsZombie(), "failed to create sentinel ROOT file");
    TNamed sentinel("sentinel_preserve", "existing output must survive failed validation");
    sentinel.Write();
    output.Close();
  }

  bool HasSentinelObject(const std::filesystem::path &path) {
    TFile input(path.string().c_str(), "READ");
    return !input.IsZombie() && input.Get("sentinel_preserve") != nullptr;
  }

  std::string ReadTextFile(const std::filesystem::path &path) {
    std::ifstream input(path);
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
  }

  // Each native mT/phi cell receives a distinct integral so merged selections are observable.
  std::string WriteRebinInput(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    Expect(!output.IsZombie(), "failed to create synthetic ROOT input");
    auto *task = output.mkdir("task");
    auto *same_directory = task->mkdir("Same");
    auto *mixed_directory = task->mkdir("Mixed");

    const int bins[7] = {2, 2, 2, 4, 1, 1, 4};
    const double minimum[7] = {-0.1, -0.1, -0.1, 0.2, 0.0, -0.5, 0.0};
    const double maximum[7] = {0.1, 0.1, 0.1, 0.6, 10.0, 0.5, exp_femto_3d::kPi};
    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    same->Sumw2();
    mixed->Sumw2();

    for (int mt_index = 0; mt_index < 4; ++mt_index) {
      const double mt = 0.25 + 0.1 * static_cast<double>(mt_index);
      for (int phi_index = 0; phi_index < 4; ++phi_index) {
        const double phi = (static_cast<double>(phi_index) + 0.5) * exp_femto_3d::kPi / 4.0;
        const double same_weight = 10.0 * static_cast<double>(mt_index + 1) + static_cast<double>(phi_index + 1);
        FillSparse(*same, mt, 5.0, phi, same_weight);
        FillSparse(*mixed, mt, 5.0, phi, 10.0 * same_weight);
      }
    }

    same_directory->cd();
    same->Write("sparse");
    mixed_directory->cd();
    mixed->Write("sparse");
    output.Close();
    return path.string();
  }

  // Populate two qn intervals with proportional SE/ME yields so qn splitting narrows ME without changing CF closure.
  std::string WriteCombinedSplitRebinInput(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    Expect(!output.IsZombie(), "failed to create combined split/rebin ROOT input");
    auto *task = output.mkdir("task");
    auto *same_directory = task->mkdir("Same");
    auto *mixed_directory = task->mkdir("Mixed");

    const int bins[7] = {2, 2, 2, 4, 1, 2, 4};
    const double minimum[7] = {-0.1, -0.1, -0.1, 0.2, 0.0, 0.0, 0.0};
    const double maximum[7] = {0.1, 0.1, 0.1, 0.6, 10.0, 2.0, exp_femto_3d::kPi};
    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    same->Sumw2();
    mixed->Sumw2();

    for (int mt_index = 0; mt_index < 4; ++mt_index) {
      const double mt = 0.25 + 0.1 * static_cast<double>(mt_index);
      for (int phi_index = 0; phi_index < 4; ++phi_index) {
        const double phi = (static_cast<double>(phi_index) + 0.5) * exp_femto_3d::kPi / 4.0;
        for (int qn_index = 0; qn_index < 2; ++qn_index) {
          const double qn = 0.5 + static_cast<double>(qn_index);
          const double same_weight = static_cast<double>((mt_index + 1) * (phi_index + 1) * (qn_index + 1));
          FillSparse(*same, mt, 5.0, phi, same_weight, qn);
          FillSparse(*mixed, mt, 5.0, phi, 10.0 * same_weight, qn);
        }
      }
    }

    same_directory->cd();
    same->Write("sparse");
    mixed_directory->cd();
    mixed->Write("sparse");
    output.Close();
    return path.string();
  }

  std::string WriteThreePhiInput(const std::filesystem::path &path) {
    TFile output(path.string().c_str(), "RECREATE");
    Expect(!output.IsZombie(), "failed to create synthetic ROOT input");
    auto *task = output.mkdir("task");
    auto *same_directory = task->mkdir("Same");
    auto *mixed_directory = task->mkdir("Mixed");

    const int bins[7] = {2, 2, 2, 4, 1, 1, 3};
    const double minimum[7] = {-0.1, -0.1, -0.1, 0.2, 0.0, -0.5, 0.0};
    const double maximum[7] = {0.1, 0.1, 0.1, 0.6, 10.0, 0.5, exp_femto_3d::kPi};
    auto same = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    auto mixed = std::make_unique<THnSparseF>("sparse", "sparse", 7, bins, minimum, maximum);
    same->Sumw2();
    mixed->Sumw2();

    for (int mt_index = 0; mt_index < 4; ++mt_index) {
      const double mt = 0.25 + 0.1 * static_cast<double>(mt_index);
      for (int phi_index = 0; phi_index < 3; ++phi_index) {
        const double phi = (static_cast<double>(phi_index) + 0.5) * exp_femto_3d::kPi / 3.0;
        FillSparse(*same, mt, 5.0, phi, 1.0);
        FillSparse(*mixed, mt, 5.0, phi, 10.0);
      }
    }

    same_directory->cd();
    same->Write("sparse");
    mixed_directory->cd();
    mixed->Write("sparse");
    output.Close();
    return path.string();
  }

  std::string WriteFactorConfig(const std::filesystem::path &path,
                                const std::string &input_root,
                                const std::string &output_directory,
                                const std::string &cf_root_name = "factor_cf.root",
                                const bool split_same_event_by_qn = false,
                                const bool split_mixed_event_by_qn = false) {
    std::ofstream output(path);
    output << "[input]\n";
    output << "input_root = \"" << input_root << "\"\n";
    output << "task_name = \"task\"\n";
    output << "same_event_subtask = \"Same\"\n";
    output << "mixed_event_subtask = \"Mixed\"\n";
    output << "sparse_object_name = \"sparse\"\n\n";
    output << "[output]\n";
    output << "output_directory = \"" << output_directory << "\"\n";
    output << "cf_root_name = \"" << cf_root_name << "\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "map_pair_phi_to_symmetric_range = false\n";
    output << "write_normalized_se_me_1d_projections = false\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "split_mixed_event_by_phi = true\n";
    output << "split_same_event_by_qn = " << (split_same_event_by_qn ? "true" : "false") << "\n";
    output << "split_mixed_event_by_qn = " << (split_mixed_event_by_qn ? "true" : "false") << "\n";
    output << "progress = false\n\n";
    output << "[build.rebin.mt]\n";
    output << "enabled = true\n";
    output << "factor = 2\n\n";
    output << "[build.rebin.phi]\n";
    output << "enabled = true\n";
    output << "factor = 2\n\n";
    output << "[fit]\n";
    output << "model = \"diag\"\n";
    output << "use_coulomb = false\n";
    output << "use_core_halo_lambda = false\n";
    output << "use_q2_baseline = false\n";
    output << "use_pml = false\n";
    output << "fit_q_max = 0.15\n\n";
    output << "[[bins.centrality]]\n";
    output << "min = 0\n";
    output << "max = 10\n";
    if (split_same_event_by_qn) {
      output << "\n[[bins.qn]]\n";
      output << "label = \"qn_low\"\n";
      output << "min = 0\n";
      output << "max = 1\n";
    }
    return path.string();
  }

  std::string WritePhiRangesMixedConfig(const std::filesystem::path &path,
                                        const std::string &input_root,
                                        const std::string &output_directory) {
    std::ofstream output(path);
    output << std::setprecision(17);
    output << "[input]\n";
    output << "input_root = \"" << input_root << "\"\n";
    output << "task_name = \"task\"\n";
    output << "same_event_subtask = \"Same\"\n";
    output << "mixed_event_subtask = \"Mixed\"\n";
    output << "sparse_object_name = \"sparse\"\n\n";
    output << "[output]\n";
    output << "output_directory = \"" << output_directory << "\"\n";
    output << "cf_root_name = \"ranges_cf.root\"\n";
    output << "log_level = \"error\"\n\n";
    output << "[build]\n";
    output << "map_pair_phi_to_symmetric_range = false\n";
    output << "write_normalized_se_me_1d_projections = false\n";
    output << "reopen_output_file_per_slice = false\n";
    output << "split_mixed_event_by_phi = true\n";
    output << "progress = false\n\n";
    output << "[build.rebin.mt]\n";
    output << "enabled = true\n";
    output << "factor = 2\n\n";
    output << "[build.rebin.phi]\n";
    output << "enabled = true\n\n";
    output << "[fit]\n";
    output << "model = \"diag\"\n";
    output << "fit_q_max = 0.15\n\n";
    output << "[[bins.centrality]]\nmin = 0\nmax = 10\n\n";
    output << "[[bins.phi]]\nmin = 0\nmax = " << exp_femto_3d::kPi / 2.0 << "\n";
    output << "\n[[bins.phi]]\nmin = " << exp_femto_3d::kPi / 2.0 << "\nmax = " << exp_femto_3d::kPi << "\n";
    return path.string();
  }

  double VisibleIntegral(const TH3D &histogram) {
    return histogram.Integral(1, histogram.GetNbinsX(), 1, histogram.GetNbinsY(), 1, histogram.GetNbinsZ());
  }

  double ExpectedSameEventIntegral(const exp_femto_3d::SliceCatalogEntry &entry) {
    const int first_mt = 2 * entry.mt_index;
    const int last_mt = first_mt + 1;
    const int first_phi = entry.is_phi_integrated ? 0 : 2 * entry.phi_index;
    const int last_phi = entry.is_phi_integrated ? 3 : first_phi + 1;
    double integral = 0.0;
    for (int mt_index = first_mt; mt_index <= last_mt; ++mt_index) {
      for (int phi_index = first_phi; phi_index <= last_phi; ++phi_index) {
        integral += 10.0 * static_cast<double>(mt_index + 1) + static_cast<double>(phi_index + 1);
      }
    }
    return integral;
  }

  void ValidateQAxesUnchanged(const TH3D &histogram, const std::string &label) {
    Expect(histogram.GetNbinsX() == 2, label + " q_out bin count changed");
    Expect(histogram.GetNbinsY() == 2, label + " q_side bin count changed");
    Expect(histogram.GetNbinsZ() == 2, label + " q_long bin count changed");
    ExpectNear(histogram.GetXaxis()->GetBinLowEdge(1), -0.1, label + " q_out low edge changed");
    ExpectNear(histogram.GetXaxis()->GetBinUpEdge(2), 0.1, label + " q_out high edge changed");
    ExpectNear(histogram.GetYaxis()->GetBinLowEdge(1), -0.1, label + " q_side low edge changed");
    ExpectNear(histogram.GetYaxis()->GetBinUpEdge(2), 0.1, label + " q_side high edge changed");
    ExpectNear(histogram.GetZaxis()->GetBinLowEdge(1), -0.1, label + " q_long low edge changed");
    ExpectNear(histogram.GetZaxis()->GetBinUpEdge(2), 0.1, label + " q_long high edge changed");
  }

  void ValidateCfHistogram(const TH3D &histogram, const std::string &label) {
    ValidateQAxesUnchanged(histogram, label);
    Expect(histogram.GetSumw2N() > 0, label + " should retain Sumw2 after normalized division");
    const int q_out_bin = histogram.GetXaxis()->FindFixBin(0.01);
    const int q_side_bin = histogram.GetYaxis()->FindFixBin(-0.01);
    const int q_long_bin = histogram.GetZaxis()->FindFixBin(0.02);
    ExpectNear(histogram.GetBinContent(q_out_bin, q_side_bin, q_long_bin), 1.0,
               label + " proportional SE/ME CF value mismatch");
    ExpectNear(VisibleIntegral(histogram), 1.0, label + " should contain one finite proportional CF bin");
    for (int x_bin = 1; x_bin <= histogram.GetNbinsX(); ++x_bin) {
      for (int y_bin = 1; y_bin <= histogram.GetNbinsY(); ++y_bin) {
        for (int z_bin = 1; z_bin <= histogram.GetNbinsZ(); ++z_bin) {
          const double error = histogram.GetBinError(x_bin, y_bin, z_bin);
          Expect(std::isfinite(error) && error >= 0.0, label + " has invalid stored bin error");
        }
      }
    }
  }

  void ExpectFitCatalogMetadata(const std::filesystem::path &path, const long long expected_entries) {
    TFile input(path.string().c_str(), "READ");
    Expect(!input.IsZombie(), "fit ROOT file should be readable: " + path.string());
    auto *tree = dynamic_cast<TTree *>(input.Get("meta/FitCatalog"));
    Expect(tree != nullptr, "FitCatalog should exist in " + path.string());
    Expect(tree->GetEntries() == expected_entries, "FitCatalog entry count mismatch in " + path.string());
    Expect(tree->GetBranch("mt_rebin_enabled") != nullptr, "FitCatalog mT rebin enabled branch missing");
    Expect(tree->GetBranch("mt_rebin_mode") != nullptr, "FitCatalog mT rebin mode branch missing");
    Expect(tree->GetBranch("phi_rebin_enabled") != nullptr, "FitCatalog phi rebin enabled branch missing");
    Expect(tree->GetBranch("phi_rebin_mode") != nullptr, "FitCatalog phi rebin mode branch missing");
    Expect(tree->GetBranch("raw_phi_low") != nullptr, "FitCatalog raw phi low branch missing");
    Expect(tree->GetBranch("display_phi_low") != nullptr, "FitCatalog display phi low branch missing");
  }

  // Build-time checks depend on the concrete ROOT axis and must throw before any slice is accepted.
  void ExpectBuildFailure(exp_femto_3d::ApplicationConfig config,
                          const exp_femto_3d::Logger &logger,
                          const std::string &output_name,
                          const std::string &message_fragment) {
    config.output.cf_root_name = output_name;
    const std::filesystem::path sentinel_path = std::filesystem::path(config.output.output_directory) / output_name;
    WriteSentinelRootFile(sentinel_path);
    bool saw_expected_error = false;
    try {
      (void)exp_femto_3d::RunBuildCf(config, logger);
    } catch (const std::runtime_error &error) {
      saw_expected_error = std::string(error.what()).find(message_fragment) != std::string::npos;
    }
    Expect(saw_expected_error, output_name + " should fail with message containing '" + message_fragment + "'");
    Expect(HasSentinelObject(sentinel_path), output_name + " should not be reset before validation failure");
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_directory =
      std::filesystem::temp_directory_path() / "exp_femto_3d_build_cf_rebin_test";
  std::filesystem::create_directories(temp_directory);
  const std::string input_root = WriteRebinInput(temp_directory / "input.root");
  const ApplicationConfig config =
      LoadApplicationConfig(WriteFactorConfig(temp_directory / "factor.toml", input_root, temp_directory.string()));
  const Logger logger(LogLevel::kError);

  const BuildCfRunStatistics statistics = RunBuildCf(config, logger);
  Expect(statistics.requested_groups == 2, "factor build should request two mT groups");
  Expect(statistics.stored_slices == 6, "two mT groups should each store phi_all plus two phi slices");
  Expect(statistics.skipped_zero_mixed_event_groups == 0 && statistics.skipped_zero_mixed_event_slices == 0
             && statistics.skipped_zero_same_event_slices == 0,
         "fully populated synthetic input should not skip slices");
  Expect(statistics.mt_input_bins == 4 && statistics.mt_output_bins == 2,
         "mT factor=2 should merge four native bins into two outputs");
  Expect(statistics.phi_input_bins == 4 && statistics.phi_output_bins == 2,
         "phi factor=2 should merge four native bins into two outputs");
  Expect(statistics.mt_rebin_enabled && statistics.mt_rebin_mode == RebinMode::kFactor,
         "mT statistics should record enabled factor rebin");
  Expect(statistics.phi_rebin_enabled && statistics.phi_rebin_mode == RebinMode::kFactor,
         "phi statistics should record enabled factor rebin");

  const std::filesystem::path output_path = temp_directory / "factor_cf.root";
  const std::vector<SliceCatalogEntry> catalog = LoadSliceCatalog(output_path.string());
  Expect(catalog.size() == 6, "SliceCatalog should contain all six stored slices");
  std::set<std::string> slice_ids;
  std::set<std::string> slice_directories;
  TFile output_file(output_path.string().c_str(), "READ");
  Expect(!output_file.IsZombie(), "factor output ROOT file should be readable");

  // Catalog rows must preserve final physical ranges, unique indexed paths, and raw projection integrals.
  for (const SliceCatalogEntry &entry : catalog) {
    Expect(entry.mt_index == 0 || entry.mt_index == 1, "catalog mT index should be an output-bin index");
    Expect(entry.phi_index == -1 || entry.phi_index == 0 || entry.phi_index == 1,
           "catalog phi index should be integrated or an output-bin index");
    Expect(entry.mt_rebin_enabled && entry.mt_rebin_mode == "factor", "catalog should record mT factor metadata");
    Expect(entry.phi_rebin_enabled && entry.phi_rebin_mode == "factor", "catalog should record phi factor metadata");
    Expect(entry.split_mixed_event_by_phi, "catalog should record split-ME mode");
    ExpectNear(entry.mt_low, entry.mt_index == 0 ? 0.2 : 0.4, "catalog mT low edge mismatch");
    ExpectNear(entry.mt_high, entry.mt_index == 0 ? 0.4 : 0.6, "catalog mT high edge mismatch");
    if (entry.is_phi_integrated) {
      Expect(entry.phi_index == -1, "phi_all slice should retain integrated index");
      ExpectNear(entry.raw_phi_low, 0.0, "phi_all low edge should cover the full normal axis");
      ExpectNear(entry.raw_phi_high, kPi, "phi_all high edge should cover the full normal axis");
      Expect(entry.slice_id.find("__phi_all") != std::string::npos, "integrated slice ID should remain phi_all");
    } else {
      const double expected_low = static_cast<double>(entry.phi_index) * kPi / 2.0;
      ExpectNear(entry.raw_phi_low, expected_low, "catalog phi low edge mismatch");
      ExpectNear(entry.raw_phi_high, expected_low + kPi / 2.0, "catalog phi high edge mismatch");
      Expect(entry.slice_id.find("__phibin_000" + std::to_string(entry.phi_index)) != std::string::npos,
             "rebinned phi slice ID should include its zero-padded output index");
    }
    Expect(entry.group_id.find("__mtbin_000" + std::to_string(entry.mt_index)) != std::string::npos,
           "rebinned group ID should include its zero-padded mT output index");
    Expect(slice_ids.insert(entry.slice_id).second, "slice IDs must be unique after rebinning");
    Expect(slice_directories.insert(entry.slice_directory).second, "slice directories must be unique after rebinning");

    auto *same_histogram = dynamic_cast<TH3D *>(output_file.Get(entry.se_object_path.c_str()));
    auto *mixed_histogram = dynamic_cast<TH3D *>(output_file.Get(entry.me_object_path.c_str()));
    auto *cf_histogram = dynamic_cast<TH3D *>(output_file.Get(entry.cf_object_path.c_str()));
    Expect(same_histogram != nullptr && mixed_histogram != nullptr,
           "catalog raw SE/ME paths should resolve to TH3D objects");
    Expect(cf_histogram != nullptr, "catalog CF path should resolve to a TH3D object");
    ValidateQAxesUnchanged(*same_histogram, entry.slice_id + " SE");
    ValidateQAxesUnchanged(*mixed_histogram, entry.slice_id + " ME");
    ValidateCfHistogram(*cf_histogram, entry.slice_id + " CF");
    const double expected_same = ExpectedSameEventIntegral(entry);
    ExpectNear(VisibleIntegral(*same_histogram), expected_same, "merged SE integral mismatch");
    ExpectNear(VisibleIntegral(*mixed_histogram), 10.0 * expected_same, "split-ME integral mismatch");
  }

  ApplicationConfig selected_fit = config;
  selected_fit.output.fit_root_name = "selected_fit.root";
  selected_fit.output.fit_summary_name = "selected_fit.tsv";
  selected_fit.output.fit_report_root_name = "selected_report.root";
  selected_fit.fit_mt_bins = {{0.2, 0.4, "mt_020_040"}};
  const FitRunStatistics selected_fit_statistics = RunFit(selected_fit, logger);
  Expect(selected_fit_statistics.catalog_slices == 6, "fit should read all catalog slices before selection");
  Expect(selected_fit_statistics.selected_slices == 3, "dynamic mT fit selection should choose one rebinned group");
  Expect(selected_fit_statistics.fitted_slices == 3, "dynamic mT fit selection should fit one phi-all plus two phi slices");
  ExpectFitCatalogMetadata(temp_directory / "selected_fit.root", 3);
  ExpectFitCatalogMetadata(temp_directory / "selected_report.root", 3);
  const std::string selected_tsv = ReadTextFile(temp_directory / "selected_fit.tsv");
  Expect(selected_tsv.find("mTRebinEnabled\tmTRebinMode\tphiRebinEnabled\tphiRebinMode") != std::string::npos,
         "fit TSV should expose rebin metadata columns");
  Expect(selected_tsv.find("\t1\tfactor\t1\tfactor\t") != std::string::npos,
         "fit TSV should record factor rebin metadata for selected slices");

  ApplicationConfig invalid_fit_selection = config;
  invalid_fit_selection.output.fit_root_name = "invalid_fit.root";
  invalid_fit_selection.fit_mt_bins = {{0.25, 0.45, "not_a_catalog_group"}};
  WriteSentinelRootFile(temp_directory / "invalid_fit.root");
  bool saw_invalid_fit_selection = false;
  try {
    (void)RunFit(invalid_fit_selection, logger);
  } catch (const std::runtime_error &error) {
    saw_invalid_fit_selection = std::string(error.what()).find("does not match any mT group") != std::string::npos;
  }
  Expect(saw_invalid_fit_selection, "invalid dynamic fit_selection.mt should fail against SliceCatalog");
  Expect(HasSentinelObject(temp_directory / "invalid_fit.root"),
         "RunFit should not reset fit output before dynamic fit_selection.mt validation");

  const ApplicationConfig mixed_ranges_config = LoadApplicationConfig(
      WritePhiRangesMixedConfig(temp_directory / "ranges.toml", input_root, temp_directory.string()));
  const BuildCfRunStatistics mixed_ranges_statistics = RunBuildCf(mixed_ranges_config, logger);
  Expect(mixed_ranges_statistics.stored_slices == 6, "factor/ranges mixed rebin should store six slices");
  Expect(mixed_ranges_statistics.mt_rebin_mode == RebinMode::kFactor
             && mixed_ranges_statistics.phi_rebin_mode == RebinMode::kRanges,
         "mixed-mode statistics should record mT factor and phi ranges");
  const auto mixed_ranges_catalog = LoadSliceCatalog((temp_directory / "ranges_cf.root").string());
  Expect(mixed_ranges_catalog[1].mt_rebin_mode == "factor" && mixed_ranges_catalog[1].phi_rebin_mode == "ranges",
         "SliceCatalog should preserve factor/ranges mixed-mode metadata");

  // Exercise the full split/rebin cross-product and verify that both denominator narrowing and metadata survive it.
  const std::string combined_input = WriteCombinedSplitRebinInput(temp_directory / "combined_input.root");
  const ApplicationConfig combined_config = LoadApplicationConfig(
      WriteFactorConfig(temp_directory / "combined.toml",
                        combined_input,
                        temp_directory.string(),
                        "combined_cf.root",
                        true,
                        true));
  const BuildCfRunStatistics combined_statistics = RunBuildCf(combined_config, logger);
  Expect(combined_statistics.requested_groups == 4 && combined_statistics.stored_slices == 12,
         "combined split/rebin build should store two qn states across two mT and three phi outputs");

  const auto combined_catalog = LoadSliceCatalog((temp_directory / "combined_cf.root").string());
  Expect(combined_catalog.size() == 12, "combined SliceCatalog should contain all stored slices");
  const SliceCatalogEntry *qn_all_entry = nullptr;
  const SliceCatalogEntry *qn_low_entry = nullptr;
  for (const SliceCatalogEntry &entry : combined_catalog) {
    Expect(entry.split_mixed_event_by_phi && entry.split_mixed_event_by_qn,
           "combined catalog should record both ME split dimensions");
    Expect(entry.mt_rebin_enabled && entry.mt_rebin_mode == "factor" && entry.phi_rebin_enabled
               && entry.phi_rebin_mode == "factor",
           "combined catalog should record both factor-rebin dimensions");
    if (entry.mt_index == 0 && entry.is_phi_integrated && entry.is_qn_integrated) {
      qn_all_entry = &entry;
    }
    if (entry.mt_index == 0 && entry.is_phi_integrated && entry.qn_label == "qn_low") {
      qn_low_entry = &entry;
    }
  }
  Expect(qn_all_entry != nullptr && qn_low_entry != nullptr,
         "combined catalog should contain matching qn-all and qn-low integrated-phi slices");
  Expect(!qn_low_entry->is_qn_integrated && qn_low_entry->qn_index == 0,
         "combined qn-specific catalog metadata mismatch");

  TFile combined_file((temp_directory / "combined_cf.root").string().c_str(), "READ");
  auto *qn_all_me = dynamic_cast<TH3D *>(combined_file.Get(qn_all_entry->me_object_path.c_str()));
  auto *qn_low_me = dynamic_cast<TH3D *>(combined_file.Get(qn_low_entry->me_object_path.c_str()));
  auto *qn_all_cf = dynamic_cast<TH3D *>(combined_file.Get(qn_all_entry->cf_object_path.c_str()));
  auto *qn_low_cf = dynamic_cast<TH3D *>(combined_file.Get(qn_low_entry->cf_object_path.c_str()));
  Expect(qn_all_me != nullptr && qn_low_me != nullptr && qn_all_cf != nullptr && qn_low_cf != nullptr,
         "combined split/rebin catalog paths should resolve");
  Expect(VisibleIntegral(*qn_low_me) < VisibleIntegral(*qn_all_me),
         "qn-specific ME should be narrower than qn-integrated ME after phi/mT rebin");
  ValidateCfHistogram(*qn_all_cf, "combined qn-all CF closure");
  ValidateCfHistogram(*qn_low_cf, "combined qn-low CF closure");

  // Physical edge, divisibility, and symmetric-map seam failures require a concrete input axis.
  ApplicationConfig nonaligned = config;
  nonaligned.build.mt_rebin.min = 0.25;
  nonaligned.build.mt_rebin.max = 0.55;
  ExpectBuildFailure(nonaligned, logger, "nonaligned.root", "must align");

  ApplicationConfig nondivisible = config;
  nondivisible.build.mt_rebin.factor = 3;
  ExpectBuildFailure(nondivisible, logger, "nondivisible.root", "does not divide");

  ApplicationConfig seam_crossing = config;
  seam_crossing.build.map_pair_phi_to_symmetric_range = true;
  seam_crossing.build.phi_rebin.mode = RebinMode::kRanges;
  seam_crossing.build.phi_rebin.factor.reset();
  seam_crossing.phi_bins = {{0.0, 3.0 * kPi / 4.0, "crosses_seam"}};
  ExpectBuildFailure(seam_crossing, logger, "seam_crossing.root", "crosses the pi/2");

  ApplicationConfig native_seam_crossing = config;
  native_seam_crossing.input.input_root = WriteThreePhiInput(temp_directory / "three_phi_input.root");
  native_seam_crossing.build.map_pair_phi_to_symmetric_range = true;
  native_seam_crossing.build.phi_rebin.configured = true;
  native_seam_crossing.build.phi_rebin.enabled = false;
  native_seam_crossing.build.phi_rebin.mode = RebinMode::kNative;
  native_seam_crossing.build.phi_rebin.factor.reset();
  native_seam_crossing.phi_bins.clear();
  ExpectBuildFailure(native_seam_crossing, logger, "native_seam_crossing.root", "crosses the pi/2");

  std::cout << "build_cf_rebin_test passed\n";
  return 0;
}
