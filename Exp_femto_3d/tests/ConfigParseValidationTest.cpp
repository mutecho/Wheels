#include "exp_femto_3d/Config.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void Expect(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string WriteFile(const std::filesystem::path& path, const std::string& contents) {
  std::ofstream output(path);
  output << contents;
  return path.string();
}

}  // namespace

int main() {
  using namespace exp_femto_3d;

  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "exp_femto_3d_config_test";
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
log_level = "debug"

[build]
map_pair_phi_to_symmetric_range = true
write_normalized_se_me_1d_projections = false
reopen_output_file_per_slice = true

[fit]
model = "diag"
use_coulomb = false
use_core_halo_lambda = true
use_q2_baseline = false
use_pml = false
fit_q_max = 0.2

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig config =
      LoadApplicationConfig(WriteFile(temp_dir / "valid.toml", valid_config));
  Expect(config.input.input_root == "/tmp/input.root", "input_root mismatch");
  Expect(config.fit.model == FitModel::kDiag, "fit model mismatch");
  Expect(config.output.log_level == LogLevel::kDebug, "log level mismatch");
  Expect(config.fit_centrality_bins.size() == 1, "fit centrality fallback failed");
  Expect(config.output.cf_root_name == "cf.root", "root extension normalization failed");
  Expect(config.output.fit_summary_name == "summary.tsv",
         "summary extension normalization failed");

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

  const ApplicationConfig overlapping_config = LoadApplicationConfig(
      WriteFile(temp_dir / "overlapping_bins.toml", overlapping_bins_config));
  Expect(overlapping_config.centrality_bins.size() == 2,
         "overlapping centrality bins should be accepted");
  Expect(overlapping_config.mt_bins.size() == 3,
         "overlapping mt bins should be accepted");
  Expect(overlapping_config.fit_centrality_bins.size() == 1,
         "fit selection centrality should parse");
  Expect(overlapping_config.fit_mt_bins.size() == 2,
         "fit selection mt should parse");

  const std::filesystem::path project_root =
      std::filesystem::path(__FILE__).parent_path().parent_path();
  const ApplicationConfig pbpb_example = LoadApplicationConfig(
      (project_root / "config/examples/pbpb_build_and_fit.toml").string());
  Expect(pbpb_example.centrality_bins.size() == 5,
         "pbpb example centrality bins should parse");
  Expect(pbpb_example.mt_bins.size() == 5, "pbpb example merged mt bins should parse");
  Expect(pbpb_example.fit_mt_bins.size() == 3,
         "pbpb example fit_selection.mt should parse");

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
  } catch (const ConfigError&) {
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
    (void)LoadApplicationConfig(
        WriteFile(temp_dir / "invalid_duplicate_bin.toml", invalid_duplicate_bin_config));
  } catch (const ConfigError&) {
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
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_fit_selection.toml",
                                          invalid_fit_selection_config));
  } catch (const ConfigError&) {
    saw_invalid_fit_selection = true;
  }
  Expect(saw_invalid_fit_selection,
         "fit_selection bins must exactly match build bins");

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
    (void)LoadApplicationConfig(
        WriteFile(temp_dir / "invalid_range.toml", invalid_range_config));
  } catch (const ConfigError&) {
    saw_invalid_range = true;
  }
  Expect(saw_invalid_range, "invalid range (min >= max) should fail");

  std::cout << "config_parse_validation_test passed\n";
  return 0;
}
