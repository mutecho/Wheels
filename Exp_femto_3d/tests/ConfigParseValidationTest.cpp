#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "exp_femto_3d/Config.h"

namespace {

  void Expect(bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  void ExpectOptionalDouble(const std::optional<double> &value,
                            const double expected,
                            const std::string &message) {
    Expect(value.has_value(), message + " should be configured");
    Expect(std::abs(*value - expected) < 1.0e-12, message + " mismatch");
  }

  void ExpectParameterTriple(const exp_femto_3d::LevyFitParameterOverride &parameter,
                             const double initial,
                             const double min,
                             const double max,
                             const std::string &name) {
    ExpectOptionalDouble(parameter.initial, initial, name + " initial");
    ExpectOptionalDouble(parameter.min, min, name + " min");
    ExpectOptionalDouble(parameter.max, max, name + " max");
  }

  std::string WriteFile(const std::filesystem::path &path, const std::string &contents) {
    std::ofstream output(path);
    output << contents;
    return path.string();
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_3d_config_test";
  std::filesystem::create_directories(temp_dir);

  const std::string valid_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"
cf_root_name = "cf"
fit_root_name = "fit"
fit_summary_name = "summary"
fit_report_directory = "/tmp/report"
fit_report_root_name = "report"
log_level = "debug"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true
split_mixed_event_by_phi = true
progress = false

[fit]
model = "diag"
use_coulomb = false
use_core_halo_lambda = true
use_q2_baseline = false
use_pml = false
fit_q_max = 0.2
map_pair_phi_to_symmetric_range = false
progress = "enabled"

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig config = LoadApplicationConfig(WriteFile(temp_dir / "valid.toml", valid_config));
  Expect(config.input.input_root == "/tmp/input.root", "input_root mismatch");
  Expect(config.fit.model == FitModel::kDiag, "fit model mismatch");
  Expect(config.output.log_level == LogLevel::kDebug, "log level mismatch");
  Expect(config.fit_centrality_bins.size() == 1, "fit centrality fallback failed");
  Expect(config.output.cf_root_name == "cf.root", "root extension normalization failed");
  Expect(config.output.fit_summary_name == "summary.tsv", "summary extension normalization failed");
  Expect(config.output.fit_report_directory == "/tmp/report", "fit report directory should parse");
  Expect(config.output.fit_report_root_name == "report.root", "fit report root extension normalization failed");
  Expect(config.build.split_mixed_event_by_phi, "ME phi split switch should parse");
  Expect(config.build.progress == ProgressMode::kDisabled, "build progress mode mismatch");
  Expect(config.fit.progress == ProgressMode::kEnabled, "fit progress mode mismatch");
  Expect(config.fit.map_pair_phi_to_symmetric_range.has_value(), "fit phi mapping override should parse");
  Expect(!*config.fit.map_pair_phi_to_symmetric_range, "fit phi mapping override should be false");
  Expect(config.fit.options.coulomb_mode == CoulombMode::kNone, "legacy use_coulomb=false should map to none");
  Expect(ToString(CoulombMode::kGamow) == "gamow", "CoulombMode string helper mismatch");
  Expect(!config.build.mt_rebin.configured && !config.build.mt_rebin.enabled,
         "legacy mT configuration should keep the explicit rebin switch absent");
  Expect(config.build.mt_rebin.mode == RebinMode::kLegacyRanges,
         "legacy [[bins.mt]] should retain legacy-range behavior");
  Expect(!config.build.phi_rebin.configured && config.build.phi_rebin.mode == RebinMode::kNative,
         "legacy phi configuration should use native bins");

  const std::string overlapping_bins_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
fit_q_max = 0.15

[[bins.centrality]]
min = 0
max = 10

[[bins.centrality]]
min = 0
max = 30

[[bins.mt]]
min = 0.2
max = 0.4

[[bins.mt]]
min = 0.4
max = 0.6

[[bins.mt]]
min = 0.2
max = 0.6

[[fit_selection.centrality]]
min = 0
max = 30

[[fit_selection.mt]]
min = 0.2
max = 0.4

[[fit_selection.mt]]
min = 0.2
max = 0.6
)toml";

  const ApplicationConfig overlapping_config =
      LoadApplicationConfig(WriteFile(temp_dir / "overlapping_bins.toml", overlapping_bins_config));
  Expect(overlapping_config.centrality_bins.size() == 2, "overlapping centrality bins should be accepted");
  Expect(overlapping_config.mt_bins.size() == 3, "overlapping mt bins should be accepted");
  Expect(overlapping_config.fit_centrality_bins.size() == 1, "fit selection centrality should parse");
  Expect(overlapping_config.fit_mt_bins.size() == 2, "fit selection mt should parse");
  Expect(overlapping_config.output.fit_report_directory == "/tmp/out",
         "fit report directory should default to output directory");
  Expect(overlapping_config.output.fit_report_root_name == "fit_report.root", "fit report root name should default");
  Expect(!overlapping_config.build.split_mixed_event_by_phi, "ME phi split should default to false");
  Expect(!overlapping_config.build.split_same_event_by_qn, "same-event qn split should default to false");
  Expect(!overlapping_config.build.split_mixed_event_by_qn, "ME qn split should default to false");
  Expect(overlapping_config.qn_bins.empty(), "qn bins should default to empty");
  Expect(overlapping_config.build.progress == ProgressMode::kAuto, "build progress should default to auto");
  Expect(overlapping_config.fit.progress == ProgressMode::kAuto, "fit progress should default to auto");
  Expect(!overlapping_config.fit.map_pair_phi_to_symmetric_range.has_value(),
         "fit phi mapping should default to follow-input when omitted");

  const std::string fit_phi_mapping_true_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "diag"
fit_q_max = 0.15
map_pair_phi_to_symmetric_range = true

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig fit_phi_mapping_true =
      LoadApplicationConfig(WriteFile(temp_dir / "fit_phi_mapping_true.toml", fit_phi_mapping_true_config));
  Expect(fit_phi_mapping_true.fit.map_pair_phi_to_symmetric_range.has_value(),
         "fit phi mapping true override should parse");
  Expect(*fit_phi_mapping_true.fit.map_pair_phi_to_symmetric_range, "fit phi mapping true override should be true");

  const std::string explicit_integrated_me_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true
split_mixed_event_by_phi = false

[fit]
model = "diag"
fit_q_max = 0.15

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig explicit_integrated_me =
      LoadApplicationConfig(WriteFile(temp_dir / "explicit_integrated_me.toml", explicit_integrated_me_config));
  Expect(!explicit_integrated_me.build.split_mixed_event_by_phi, "explicit integrated ME mode should parse false");

  const std::filesystem::path project_root = std::filesystem::path(__FILE__).parent_path().parent_path();
  const ApplicationConfig pbpb_example =
      LoadApplicationConfig((project_root / "config/pbpb_build_and_fit.toml").string());
  Expect(pbpb_example.centrality_bins.size() == 5, "pbpb example centrality bins should parse");
  Expect(pbpb_example.mt_bins.size() == 5, "pbpb example merged mt bins should parse");
  Expect(pbpb_example.fit_mt_bins.size() == 3, "pbpb example fit_selection.mt should parse");
  Expect(pbpb_example.build.progress == ProgressMode::kAuto, "pbpb build progress should parse");
  Expect(pbpb_example.fit.progress == ProgressMode::kAuto, "pbpb fit progress should parse");
  Expect(pbpb_example.fit.options.coulomb_mode == CoulombMode::kGamow, "legacy use_coulomb=true should map to gamow");

  // Keep the Wenya production contract explicit across qn splitting, ME policy, rebin, and fit overrides.
  const ApplicationConfig wenya_qn_integrated =
      LoadApplicationConfig((project_root / "config/pbpb_wenya_lhc23_qn_integrated_build_and_fit.toml").string());
  Expect(wenya_qn_integrated.input.input_root
             == "/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb/wenya/3Dfemto_cent_mt_q2_phi_LHC23_merge.root",
         "wenya qn-integrated input path mismatch");
  Expect(wenya_qn_integrated.input.task_name == "femto-dream-pair-task-track-track",
         "wenya qn-integrated task name mismatch");
  Expect(wenya_qn_integrated.input.same_event_subtask == "SameEvent_3Dqn",
         "wenya qn-integrated same-event subtask mismatch");
  Expect(wenya_qn_integrated.input.mixed_event_subtask == "MixedEvent_3Dqn",
         "wenya qn-integrated mixed-event subtask mismatch");
  Expect(wenya_qn_integrated.input.sparse_object_name == "relPair3dRmTMultPercentileQnPairphi",
         "wenya qn-integrated sparse name mismatch");
  Expect(wenya_qn_integrated.output.output_directory == "/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb/wenya",
         "wenya qn-integrated output directory mismatch");
  Expect(wenya_qn_integrated.output.cf_root_name == "EP_dependence_CF_wenya_lhc23_merge_qn_split_plus_integrated.root",
         "wenya qn-split CF output name mismatch");
  Expect(wenya_qn_integrated.output.fit_root_name
             == "EP_dependence_CF_full_fit_wenya_lhc23_merge_qn_split_plus_integrated.root",
         "wenya qn-split fit output name mismatch");
  Expect(wenya_qn_integrated.output.fit_summary_name == "fit_summary_wenya_lhc23_merge_qn_split_plus_integrated.tsv",
         "wenya qn-split fit summary name mismatch");
  Expect(wenya_qn_integrated.output.fit_report_root_name
             == "EP_dependence_CF_full_fit_report_wenya_lhc23_merge_qn_split_plus_integrated.root",
         "wenya qn-split fit report name mismatch");
  Expect(wenya_qn_integrated.centrality_bins.size() == 4,
         "wenya qn-integrated centrality bins should stop at 80 percent");
  Expect(
      wenya_qn_integrated.centrality_bins.back().min == 50.0 && wenya_qn_integrated.centrality_bins.back().max == 80.0,
      "wenya qn-integrated final centrality bin should be 50-80");
  Expect(wenya_qn_integrated.mt_bins.size() == 5, "wenya qn-integrated mT bins should include merged full mT");
  Expect(wenya_qn_integrated.fit_centrality_bins.size() == 3,
         "wenya qn-integrated fit centrality selection should follow PbPb subset");
  Expect(wenya_qn_integrated.fit_mt_bins.size() == 3, "wenya qn-integrated fit mT selection should follow PbPb subset");
  Expect(wenya_qn_integrated.build.split_same_event_by_qn, "wenya config should enable same-event qn splitting");
  Expect(wenya_qn_integrated.qn_bins.size() == 3, "wenya config should define qn1/qn2/qn3 bins");
  Expect(wenya_qn_integrated.qn_bins[0].label == "qn1" && wenya_qn_integrated.qn_bins[0].min == 0.0
             && wenya_qn_integrated.qn_bins[0].max == 3.0,
         "wenya qn1 bin mismatch");
  Expect(wenya_qn_integrated.qn_bins[1].label == "qn2" && wenya_qn_integrated.qn_bins[1].min == 3.0
             && wenya_qn_integrated.qn_bins[1].max == 7.0,
         "wenya qn2 bin mismatch");
  Expect(wenya_qn_integrated.qn_bins[2].label == "qn3" && wenya_qn_integrated.qn_bins[2].min == 7.0
             && wenya_qn_integrated.qn_bins[2].max == 10.0,
         "wenya qn3 bin mismatch");
  Expect(!wenya_qn_integrated.build.map_pair_phi_to_symmetric_range,
         "wenya qn-integrated build should keep raw phi mapping");
  Expect(wenya_qn_integrated.build.split_mixed_event_by_phi,
         "wenya qn-integrated build should keep its configured phi-split ME denominator");
  Expect(!wenya_qn_integrated.build.split_mixed_event_by_qn,
         "wenya qn-integrated build should keep qn-integrated ME denominator");
  Expect(wenya_qn_integrated.build.phi_rebin.configured && !wenya_qn_integrated.build.phi_rebin.enabled
             && wenya_qn_integrated.build.phi_rebin.mode == RebinMode::kNative,
         "wenya qn-integrated build should keep native phi bins");
  Expect(wenya_qn_integrated.build.mt_rebin.configured && wenya_qn_integrated.build.mt_rebin.enabled
             && wenya_qn_integrated.build.mt_rebin.mode == RebinMode::kRanges,
         "wenya qn-integrated build should use explicit mT ranges");
  Expect(wenya_qn_integrated.fit.model == FitModel::kFull, "wenya qn-integrated fit should use full model");
  Expect(wenya_qn_integrated.fit.options.coulomb_mode == CoulombMode::kGamow,
         "wenya qn-integrated fit should keep legacy PbPb Coulomb setting");
  Expect(wenya_qn_integrated.fit.options.use_core_halo_lambda, "wenya qn-integrated fit should keep core-halo lambda");
  Expect(wenya_qn_integrated.fit.options.use_q2_baseline, "wenya qn-integrated fit should keep q2 baseline");
  Expect(wenya_qn_integrated.fit.options.use_pml, "wenya qn-integrated fit should keep PML");
  Expect(!wenya_qn_integrated.fit.map_pair_phi_to_symmetric_range.has_value(),
         "wenya qn-integrated fit should follow input CF phi metadata");
  ExpectParameterTriple(wenya_qn_integrated.fit.options.parameters.alpha, 1.5, 0.5, 2.0,
                        "wenya alpha override");
  ExpectOptionalDouble(wenya_qn_integrated.fit.options.parameters.alpha.fixed_value, 2.0,
                       "wenya alpha fixed value");

  const auto mode_config = [](const std::string &fit_coulomb_lines) {
    return R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "diag"
)toml" + fit_coulomb_lines
           + R"toml(
use_core_halo_lambda = true
use_q2_baseline = false
use_pml = false
fit_q_max = 0.15

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";
  };

  const ApplicationConfig explicit_none =
      LoadApplicationConfig(WriteFile(temp_dir / "coulomb_none.toml", mode_config("coulomb_mode = \"none\"\n")));
  Expect(explicit_none.fit.options.coulomb_mode == CoulombMode::kNone, "explicit none Coulomb mode should parse");

  const ApplicationConfig explicit_gamow =
      LoadApplicationConfig(WriteFile(temp_dir / "coulomb_gamow.toml", mode_config("coulomb_mode = \"gamow\"\n")));
  Expect(explicit_gamow.fit.options.coulomb_mode == CoulombMode::kGamow, "explicit gamow Coulomb mode should parse");

  const ApplicationConfig explicit_finite = LoadApplicationConfig(
      WriteFile(temp_dir / "coulomb_finite.toml",
                mode_config("coulomb_mode = \"finite_source\"\nfinite_source_mode = \"iterative_1d\"\n")));
  Expect(explicit_finite.fit.options.coulomb_mode == CoulombMode::kFiniteSource,
         "explicit finite-source Coulomb mode should parse");
  Expect(explicit_finite.fit.options.finite_source_mode == FiniteSourceMode::kIterative1D,
         "explicit iterative finite-source mode should parse");

  const ApplicationConfig agreeing_legacy = LoadApplicationConfig(
      WriteFile(temp_dir / "coulomb_agree.toml", mode_config("use_coulomb = true\ncoulomb_mode = \"gamow\"\n")));
  Expect(agreeing_legacy.fit.options.coulomb_mode == CoulombMode::kGamow,
         "agreeing legacy and explicit Coulomb fields should parse");

  const auto expect_config_error = [&](const std::string &name, const std::string &contents) {
    bool saw_error = false;
    try {
      (void)LoadApplicationConfig(WriteFile(temp_dir / name, contents));
    } catch (const ConfigError &) {
      saw_error = true;
    }
    Expect(saw_error, name + " should fail config validation");
  };

  // Rebin contract fixtures isolate table-level validation from ROOT-axis boundary checks.
  const auto rebin_config = [](const std::string &rebin_tables, const std::string &bin_tables = std::string()) {
    return R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true
)toml" + rebin_tables
           + R"toml(

[fit]
model = "diag"
fit_q_max = 0.15

[[bins.centrality]]
min = 0
max = 10
)toml" + bin_tables;
  };

  const ApplicationConfig explicit_disabled =
      LoadApplicationConfig(WriteFile(temp_dir / "rebin_explicit_disabled.toml", rebin_config(R"toml(
[build.rebin.mt]
enabled = false

[build.rebin.phi]
enabled = false
)toml")));
  Expect(explicit_disabled.build.mt_rebin.configured && !explicit_disabled.build.mt_rebin.enabled,
         "explicit disabled mT switch should parse");
  Expect(explicit_disabled.build.mt_rebin.mode == RebinMode::kNative, "disabled mT rebin should select native mode");
  Expect(explicit_disabled.build.phi_rebin.configured && !explicit_disabled.build.phi_rebin.enabled,
         "explicit disabled phi switch should parse");
  Expect(explicit_disabled.build.phi_rebin.mode == RebinMode::kNative, "disabled phi rebin should select native mode");
  Expect(explicit_disabled.mt_bins.empty() && explicit_disabled.phi_bins.empty(),
         "native axes should not require configured physical ranges");

  const ApplicationConfig factor_rebin =
      LoadApplicationConfig(WriteFile(temp_dir / "rebin_factor.toml", rebin_config(R"toml(
[build.rebin.mt]
enabled = true
factor = 2
min = 0.2
max = 0.6

[build.rebin.phi]
enabled = true
factor = 4
)toml")));
  Expect(factor_rebin.build.mt_rebin.enabled && factor_rebin.build.mt_rebin.mode == RebinMode::kFactor,
         "mT factor mode should parse");
  Expect(factor_rebin.build.mt_rebin.factor == 2 && factor_rebin.build.mt_rebin.min == 0.2
             && factor_rebin.build.mt_rebin.max == 0.6,
         "mT factor window should parse");
  Expect(factor_rebin.build.phi_rebin.enabled && factor_rebin.build.phi_rebin.mode == RebinMode::kFactor
             && factor_rebin.build.phi_rebin.factor == 4,
         "phi factor mode should parse");

  const ApplicationConfig range_rebin = LoadApplicationConfig(WriteFile(temp_dir / "rebin_ranges.toml",
                                                                        rebin_config(R"toml(
[build.rebin.mt]
enabled = true

[build.rebin.phi]
enabled = true
)toml",
                                                                                     R"toml(

[[bins.mt]]
min = 0.2
max = 0.4
label = "mt_low"

[[bins.mt]]
min = 0.2
max = 0.6
label = "mt_all"

[[bins.phi]]
min = 0.0
max = 1.0
label = "phi_low"

[[bins.phi]]
min = 0.0
max = 2.0
label = "phi_wide"
)toml")));
  Expect(range_rebin.build.mt_rebin.mode == RebinMode::kRanges && range_rebin.mt_bins.size() == 2,
         "overlapping mT ranges should parse in explicit range mode");
  Expect(range_rebin.build.phi_rebin.mode == RebinMode::kRanges && range_rebin.phi_bins.size() == 2,
         "overlapping phi ranges should parse in explicit range mode");

  expect_config_error("rebin_disabled_factor.toml", rebin_config(R"toml(
[build.rebin.mt]
enabled = false
factor = 2

[build.rebin.phi]
enabled = false
)toml"));
  expect_config_error("rebin_disabled_ranges.toml",
                      rebin_config(R"toml(
[build.rebin.mt]
enabled = false

[build.rebin.phi]
enabled = false
)toml",
                                   R"toml(

[[bins.mt]]
min = 0.2
max = 0.4
)toml"));
  expect_config_error("rebin_factor_too_small.toml", rebin_config(R"toml(
[build.rebin.mt]
enabled = true
factor = 1

[build.rebin.phi]
enabled = false
)toml"));
  expect_config_error("rebin_one_sided_window.toml", rebin_config(R"toml(
[build.rebin.mt]
enabled = true
factor = 2
min = 0.2

[build.rebin.phi]
enabled = false
)toml"));
  expect_config_error("rebin_factor_and_ranges.toml",
                      rebin_config(R"toml(
[build.rebin.mt]
enabled = true
factor = 2

[build.rebin.phi]
enabled = false
)toml",
                                   R"toml(

[[bins.mt]]
min = 0.2
max = 0.4
)toml"));
  expect_config_error("rebin_duplicate_ranges.toml",
                      rebin_config(R"toml(
[build.rebin.mt]
enabled = true

[build.rebin.phi]
enabled = false
)toml",
                                   R"toml(

[[bins.mt]]
min = 0.2
max = 0.4

[[bins.mt]]
min = 0.2
max = 0.4
)toml"));

  expect_config_error("coulomb_conflict_false_gamow.toml",
                      mode_config("use_coulomb = false\ncoulomb_mode = \"gamow\"\n"));
  expect_config_error("coulomb_conflict_true_none.toml", mode_config("use_coulomb = true\ncoulomb_mode = \"none\"\n"));
  expect_config_error("coulomb_conflict_true_finite.toml",
                      mode_config("use_coulomb = true\ncoulomb_mode = \"finite_source\"\n"));
  expect_config_error("coulomb_invalid_mode.toml", mode_config("coulomb_mode = \"bogus\"\n"));
  expect_config_error("finite_mode_without_finite_source.toml",
                      mode_config("coulomb_mode = \"gamow\"\nfinite_source_mode = \"fixed_1d\"\n"));
  expect_config_error("finite_mode_invalid.toml",
                      mode_config("coulomb_mode = \"finite_source\"\nfinite_source_mode = \"bogus\"\n"));

  const std::string parameter_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
use_core_halo_lambda = true
use_q2_baseline = true
use_pml = false
fit_q_max = 0.15

[fit.parameters.norm]
initial = 1.02
min = 0.8
max = 1.2

[fit.parameters.lambda]
initial = 0.55
min = 0.1
max = 0.9
fixed_value = 0.65

[fit.parameters.rout2]
initial = 30.0
min = 0.5
max = 500.0

[fit.parameters.rside2]
initial = 31.0
min = 0.5
max = 500.0

[fit.parameters.rlong2]
initial = 32.0
min = 0.5
max = 500.0

[fit.parameters.routside2]
initial = 1.0
min = -100.0
max = 100.0

[fit.parameters.routlong2]
initial = -2.0
min = -100.0
max = 100.0

[fit.parameters.rsidelong2]
initial = 3.0
min = -100.0
max = 100.0

[fit.parameters.alpha]
initial = 1.4
min = 0.8
max = 1.8
fixed_value = 1.2

[fit.parameters.baseline_q2]
initial = 0.5
min = -10.0
max = 20.0

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig parameter_overrides =
      LoadApplicationConfig(WriteFile(temp_dir / "parameter_overrides.toml", parameter_config));
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.norm, 1.02, 0.8, 1.2, "norm");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.lambda, 0.55, 0.1, 0.9, "lambda");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.rout2, 30.0, 0.5, 500.0, "rout2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.rside2, 31.0, 0.5, 500.0, "rside2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.rlong2, 32.0, 0.5, 500.0, "rlong2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.routside2, 1.0, -100.0, 100.0, "routside2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.routlong2, -2.0, -100.0, 100.0, "routlong2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.rsidelong2, 3.0, -100.0, 100.0, "rsidelong2");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.alpha, 1.4, 0.8, 1.8, "alpha");
  ExpectParameterTriple(parameter_overrides.fit.options.parameters.baseline_q2, 0.5, -10.0, 20.0, "baseline_q2");
  ExpectOptionalDouble(parameter_overrides.fit.options.parameters.lambda.fixed_value, 0.65, "lambda fixed value");
  ExpectOptionalDouble(parameter_overrides.fit.options.parameters.alpha.fixed_value, 1.2, "alpha fixed value");

  const auto parameter_error_config = [](const std::string &parameter_lines) {
    return R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
use_core_halo_lambda = true
use_q2_baseline = true
use_pml = false
fit_q_max = 0.15

)toml" + parameter_lines + R"toml(
[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";
  };

  expect_config_error("parameter_unknown_name.toml",
                      parameter_error_config("[fit.parameters.bogus]\ninitial = 1.0\n"));
  expect_config_error("parameter_unknown_field.toml",
                      parameter_error_config("[fit.parameters.lambda]\nstart = 0.5\n"));
  expect_config_error("parameter_non_finite.toml",
                      parameter_error_config("[fit.parameters.lambda]\ninitial = inf\n"));
  expect_config_error("parameter_single_limit.toml",
                      parameter_error_config("[fit.parameters.lambda]\nmin = 0.0\n"));
  expect_config_error("parameter_invalid_limit.toml",
                      parameter_error_config("[fit.parameters.lambda]\nmin = 1.0\nmax = 0.0\n"));
  expect_config_error("parameter_fixed_non_fixable.toml",
                      parameter_error_config("[fit.parameters.norm]\nfixed_value = 1.0\n"));
  expect_config_error("parameter_fixed_out_of_default_bounds.toml",
                      parameter_error_config("[fit.parameters.lambda]\nfixed_value = 1.2\n"));

  const auto parameter_error_config_with_switches =
      [](const bool use_core_halo_lambda, const bool use_q2_baseline, const std::string &parameter_lines) {
        return R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = false
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
use_core_halo_lambda = )toml"
               + std::string(use_core_halo_lambda ? "true\n" : "false\n") + "use_q2_baseline = "
               + std::string(use_q2_baseline ? "true\n" : "false\n") + R"toml(
use_pml = false
fit_q_max = 0.15

)toml" + parameter_lines + R"toml(
[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";
      };

  expect_config_error("parameter_lambda_disabled.toml",
                      parameter_error_config_with_switches(false, true, "[fit.parameters.lambda]\ninitial = 0.6\n"));
  expect_config_error(
      "parameter_baseline_disabled.toml",
      parameter_error_config_with_switches(true, false, "[fit.parameters.baseline_q2]\ninitial = 1.0\n"));

  const std::string invalid_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
fit_q_max = -1.0

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_invalid = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid.toml", invalid_config));
  } catch (const ConfigError &) {
    saw_invalid = true;
  }
  Expect(saw_invalid, "invalid config should fail");

  const std::string invalid_duplicate_bin_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
fit_q_max = 0.2

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_duplicate = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_duplicate_bin.toml", invalid_duplicate_bin_config));
  } catch (const ConfigError &) {
    saw_duplicate = true;
  }
  Expect(saw_duplicate, "duplicate bins should fail");

  const std::string invalid_fit_selection_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
fit_q_max = 0.2

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4

[[fit_selection.mt]]
min = 0.2
max = 0.5
)toml";

  bool saw_invalid_fit_selection = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_fit_selection.toml", invalid_fit_selection_config));
  } catch (const ConfigError &) {
    saw_invalid_fit_selection = true;
  }
  Expect(saw_invalid_fit_selection, "fit_selection bins must exactly match build bins");

  const std::string invalid_qn_split_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true
split_same_event_by_qn = true

[fit]
model = "full"
fit_q_max = 0.2

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_invalid_qn_split = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_qn_split.toml", invalid_qn_split_config));
  } catch (const ConfigError &) {
    saw_invalid_qn_split = true;
  }
  Expect(saw_invalid_qn_split, "qn split without qn bins should fail");

  const std::string invalid_mixed_qn_split_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true
split_mixed_event_by_qn = true

[fit]
model = "full"
fit_q_max = 0.2

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4

[[bins.qn]]
label = "qn1"
min = 0
max = 3
)toml";

  bool saw_invalid_mixed_qn_split = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_mixed_qn_split.toml", invalid_mixed_qn_split_config));
  } catch (const ConfigError &) {
    saw_invalid_mixed_qn_split = true;
  }
  Expect(saw_invalid_mixed_qn_split, "ME qn split without same-event qn split should fail");

  const std::string invalid_range_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "full"
fit_q_max = 0.2

[[bins.centrality]]
min = 10
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_invalid_range = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_range.toml", invalid_range_config));
  } catch (const ConfigError &) {
    saw_invalid_range = true;
  }
  Expect(saw_invalid_range, "invalid range (min >= max) should fail");

  std::cout << "config_parse_validation_test passed\n";
  return 0;
}
