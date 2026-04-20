#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "exp_femto_1d/Config.h"

namespace {

  void Expect(const bool condition, const std::string &message) {
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
  using namespace exp_femto_1d;

  const std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "exp_femto_1d_config_test";
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
norm_low = 0.5
norm_high = 0.8
kstar_min = 0.0
kstar_max = 0.8
reopen_output_file_per_slice = true
progress = false

[fit]
fit_kstar_max = 0.2
use_coulomb = false
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
  Expect(config.output.log_level == LogLevel::kDebug, "log level mismatch");
  Expect(config.output.cf_root_name == "cf.root", "root extension normalization failed");
  Expect(config.output.fit_summary_name == "summary.tsv", "summary extension normalization failed");
  Expect(config.build.progress == ProgressMode::kDisabled, "build progress mode mismatch");
  Expect(config.fit.progress == ProgressMode::kEnabled, "fit progress mode mismatch");
  Expect(config.fit_centrality_bins.size() == 1, "fit centrality fallback failed");
  Expect(config.fit_mt_bins.size() == 1, "fit mT fallback failed");

  const std::string progress_alias_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]
progress = "enabled"

[fit]
progress = "disabled"

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  const ApplicationConfig progress_alias =
      LoadApplicationConfig(WriteFile(temp_dir / "progress_alias.toml", progress_alias_config));
  Expect(progress_alias.build.progress == ProgressMode::kEnabled, "enabled alias should parse");
  Expect(progress_alias.fit.progress == ProgressMode::kDisabled, "disabled alias should parse");

  const std::string invalid_fit_limit_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]

[fit]
fit_kstar_max = -1

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_invalid_fit_limit = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_fit_limit.toml", invalid_fit_limit_config));
  } catch (const ConfigError &) {
    saw_invalid_fit_limit = true;
  }
  Expect(saw_invalid_fit_limit, "negative fit_kstar_max should fail");

  const std::string invalid_parameter_limits = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]

[fit]
baseline_p0_min = 1.1
baseline_p0_max = 0.9

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4
)toml";

  bool saw_invalid_limits = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_limits.toml", invalid_parameter_limits));
  } catch (const ConfigError &) {
    saw_invalid_limits = true;
  }
  Expect(saw_invalid_limits, "invalid fit parameter limits should fail");

  const std::string duplicate_bins_config = R"toml(
[input]
input_root = "/tmp/input.root"
task_name = "task"
same_event_subtask = "Same"
mixed_event_subtask = "Mixed"
sparse_object_name = "sparse"

[output]
output_directory = "/tmp/out"

[build]

[fit]

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

  bool saw_duplicate_bins = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "duplicate_bins.toml", duplicate_bins_config));
  } catch (const ConfigError &) {
    saw_duplicate_bins = true;
  }
  Expect(saw_duplicate_bins, "duplicate bins should fail");

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

[fit]

[[bins.centrality]]
min = 0
max = 10

[[bins.mt]]
min = 0.2
max = 0.4

[[fit_selection.mt]]
min = 0.3
max = 0.4
)toml";

  bool saw_invalid_fit_selection = false;
  try {
    (void)LoadApplicationConfig(WriteFile(temp_dir / "invalid_fit_selection.toml", invalid_fit_selection_config));
  } catch (const ConfigError &) {
    saw_invalid_fit_selection = true;
  }
  Expect(saw_invalid_fit_selection, "fit selection must match build bins exactly");

  const std::filesystem::path project_root = std::filesystem::path(__FILE__).parent_path().parent_path();
  const ApplicationConfig example_config =
      LoadApplicationConfig((project_root / "config/examples/exp_femto_1d.example.toml").string());
  Expect(example_config.centrality_bins.size() == 2, "example config centrality bins should parse");
  Expect(example_config.mt_bins.size() == 2, "example config mt bins should parse");

  const ApplicationConfig pbpb_config = LoadApplicationConfig((project_root / "config/pbpb_build_and_fit.toml").string());
  Expect(pbpb_config.centrality_bins.size() == 5, "pbpb config centrality bins should parse");
  Expect(pbpb_config.mt_bins.size() == 5, "pbpb config mt bins should parse");
  Expect(pbpb_config.fit_mt_bins.size() == 3, "pbpb fit_selection.mt should parse");

  return 0;
}
