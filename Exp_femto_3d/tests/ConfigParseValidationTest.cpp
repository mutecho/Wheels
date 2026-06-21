#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "exp_femto_3d/Config.h"

namespace {

  void Expect(bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
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
  Expect(overlapping_config.output.fit_report_root_name == "fit_report.root",
         "fit report root name should default");
  Expect(!overlapping_config.build.split_mixed_event_by_phi, "ME phi split should default to false");
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
  Expect(*fit_phi_mapping_true.fit.map_pair_phi_to_symmetric_range,
         "fit phi mapping true override should be true");

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
  const ApplicationConfig pbpb_example = LoadApplicationConfig((project_root / "config/pbpb_build_and_fit.toml").string());
  Expect(pbpb_example.centrality_bins.size() == 5, "pbpb example centrality bins should parse");
  Expect(pbpb_example.mt_bins.size() == 5, "pbpb example merged mt bins should parse");
  Expect(pbpb_example.fit_mt_bins.size() == 3, "pbpb example fit_selection.mt should parse");
  Expect(pbpb_example.build.progress == ProgressMode::kAuto, "pbpb build progress should parse");
  Expect(pbpb_example.fit.progress == ProgressMode::kAuto, "pbpb fit progress should parse");
  Expect(pbpb_example.fit.options.coulomb_mode == CoulombMode::kGamow,
         "legacy use_coulomb=true should map to gamow");

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
)toml" + fit_coulomb_lines + R"toml(
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

  const ApplicationConfig explicit_finite = LoadApplicationConfig(WriteFile(
      temp_dir / "coulomb_finite.toml",
      mode_config("coulomb_mode = \"finite_source\"\nfinite_source_mode = \"iterative_1d\"\n")));
  Expect(explicit_finite.fit.options.coulomb_mode == CoulombMode::kFiniteSource,
         "explicit finite-source Coulomb mode should parse");
  Expect(explicit_finite.fit.options.finite_source_mode == FiniteSourceMode::kIterative1D,
         "explicit iterative finite-source mode should parse");

  const ApplicationConfig agreeing_legacy = LoadApplicationConfig(WriteFile(
      temp_dir / "coulomb_agree.toml", mode_config("use_coulomb = true\ncoulomb_mode = \"gamow\"\n")));
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

  expect_config_error("coulomb_conflict_false_gamow.toml",
                      mode_config("use_coulomb = false\ncoulomb_mode = \"gamow\"\n"));
  expect_config_error("coulomb_conflict_true_none.toml",
                      mode_config("use_coulomb = true\ncoulomb_mode = \"none\"\n"));
  expect_config_error("coulomb_conflict_true_finite.toml",
                      mode_config("use_coulomb = true\ncoulomb_mode = \"finite_source\"\n"));
  expect_config_error("coulomb_invalid_mode.toml", mode_config("coulomb_mode = \"bogus\"\n"));
  expect_config_error("finite_mode_without_finite_source.toml",
                      mode_config("coulomb_mode = \"gamow\"\nfinite_source_mode = \"fixed_1d\"\n"));
  expect_config_error("finite_mode_invalid.toml",
                      mode_config("coulomb_mode = \"finite_source\"\nfinite_source_mode = \"bogus\"\n"));

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
