#include "femto3d/Config.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string WriteFile(const std::filesystem::path& path, const std::string& content) {
  std::ofstream stream(path);
  if (!stream.is_open()) {
    throw std::runtime_error("Failed to open file for writing: " + path.string());
  }
  stream << content;
  stream.close();
  return path.string();
}

bool ExpectConfigError(const std::string& config_path) {
  try {
    (void)femto3d::LoadApplicationConfig(config_path);
  } catch (const femto3d::ConfigError&) {
    return true;
  }
  return false;
}

bool ExpectCliRuntimeError(const std::vector<std::string>& args) {
  std::vector<std::vector<char>> storage;
  storage.reserve(args.size());
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (const std::string& token : args) {
    storage.emplace_back(token.begin(), token.end());
    storage.back().push_back('\0');
    argv.push_back(storage.back().data());
  }

  try {
    (void)femto3d::ParseCliArgs(static_cast<int>(argv.size()), argv.data());
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

femto3d::CliOptions ParseCli(const std::vector<std::string>& args) {
  std::vector<std::vector<char>> storage;
  storage.reserve(args.size());
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (const std::string& token : args) {
    storage.emplace_back(token.begin(), token.end());
    storage.back().push_back('\0');
    argv.push_back(storage.back().data());
  }

  return femto3d::ParseCliArgs(static_cast<int>(argv.size()), argv.data());
}

}  // namespace

int main() {
  if (!ExpectCliRuntimeError({"eventgen_femto_3d"})) {
    std::cerr << "Expected missing --config to fail CLI parse.\n";
    return 1;
  }
  if (!ExpectCliRuntimeError({"eventgen_femto_3d", "--config"})) {
    std::cerr << "Expected missing --config value to fail CLI parse.\n";
    return 2;
  }
  if (!ExpectCliRuntimeError(
          {"eventgen_femto_3d", "--config", "cfg.toml", "--bogus"})) {
    std::cerr << "Expected unknown CLI flag to fail CLI parse.\n";
    return 3;
  }

  const std::filesystem::path temp_dir =
      std::filesystem::temp_directory_path() / "eventgen_femto_3d_config_parse_test";
  std::filesystem::create_directories(temp_dir);

  const std::string valid_config = R"toml(
[input]
schema = "blastwave_flat_trees"
input_root = "input.root"

[output]
output_root = "output.root"
)toml";

  const femto3d::ApplicationConfig valid =
      femto3d::LoadApplicationConfig(WriteFile(temp_dir / "valid.toml", valid_config));
  if (valid.input_schema != femto3d::InputSchema::kBlastwaveFlatTrees) {
    std::cerr << "Expected blastwave schema in valid config.\n";
    return 4;
  }
  if (valid.analysis.event_plane.use_internal_reconstruction ||
      !valid.analysis.event_plane.fallback_to_input_branch) {
    std::cerr << "Expected blastwave schema defaults to prefer input psi2.\n";
    return 5;
  }
  if (!valid.analysis.projection_fit
           .accept_hbt_central_value_only_for_summary) {
    std::cerr << "Expected HBT summary central-value-only fallback to default "
                 "to enabled.\n";
    return 16;
  }

  const std::string legacy_config = R"toml(
[input]
schema = "legacy_vector_tree"
input_root = "legacy_input.root"

[output]
output_root = "legacy_output.root"
)toml";

  const femto3d::ApplicationConfig legacy =
      femto3d::LoadApplicationConfig(WriteFile(temp_dir / "legacy.toml", legacy_config));
  if (legacy.input_schema != femto3d::InputSchema::kLegacyVectorTree) {
    std::cerr << "Expected legacy schema in legacy config.\n";
    return 6;
  }
  if (!legacy.analysis.event_plane.use_internal_reconstruction ||
      legacy.analysis.event_plane.fallback_to_input_branch) {
    std::cerr << "Expected legacy schema defaults to keep internal reconstruction.\n";
    return 7;
  }

  const std::string missing_input_key = R"toml(
[input]
schema = "blastwave_flat_trees"

[output]
output_root = "output.root"
)toml";

  if (!ExpectConfigError(
          WriteFile(temp_dir / "missing_input_key.toml", missing_input_key))) {
    std::cerr << "Expected missing input_root key to fail config parse.\n";
    return 8;
  }

  const std::string invalid_schema = R"toml(
[input]
schema = "unknown_schema"
input_root = "input.root"

[output]
output_root = "output.root"
)toml";

  if (!ExpectConfigError(WriteFile(temp_dir / "invalid_schema.toml", invalid_schema))) {
    std::cerr << "Expected invalid input schema to fail config parse.\n";
    return 9;
  }

  const std::string empty_bins = R"toml(
[input]
schema = "legacy_vector_tree"
input_root = "input.root"

[output]
output_root = "output.root"

[bins]
centrality = []
mt = [{ min = 0.2, max = 0.4 }]
phi = [{ min = -1.0, max = 1.0 }]
)toml";

  if (!ExpectConfigError(WriteFile(temp_dir / "empty_bins.toml", empty_bins))) {
    std::cerr << "Expected empty centrality bins to fail config parse.\n";
    return 10;
  }

  const std::string overlap_bins = R"toml(
[input]
schema = "legacy_vector_tree"
input_root = "input.root"

[output]
output_root = "output.root"

[bins]
centrality = [{ min = 0.0, max = 10.0 }, { min = 5.0, max = 20.0 }]
mt = [{ min = 0.2, max = 0.4 }]
phi = [{ min = -1.0, max = 1.0 }]
)toml";

  if (!ExpectConfigError(WriteFile(temp_dir / "overlap_bins.toml", overlap_bins))) {
    std::cerr << "Expected overlapping bins to fail config parse.\n";
    return 11;
  }

  const std::string invalid_range = R"toml(
[input]
schema = "legacy_vector_tree"
input_root = "input.root"

[output]
output_root = "output.root"

[bins]
centrality = [{ min = 0.0, max = 10.0 }]
mt = [{ min = 0.4, max = 0.2 }]
phi = [{ min = -1.0, max = 1.0 }]
)toml";

  if (!ExpectConfigError(WriteFile(temp_dir / "invalid_range.toml", invalid_range))) {
    std::cerr << "Expected invalid mt range to fail config parse.\n";
    return 12;
  }

  const femto3d::CliOptions overrides = ParseCli(
      {"eventgen_femto_3d",
       "--config",
       "config.toml",
       "--input-root",
       "/tmp/override_input.root",
       "--output-root",
       "/tmp/override_output.root",
       "--input-schema",
       "legacy_vector_tree"});
  femto3d::ApplicationConfig overridden = valid;
  femto3d::ApplyCliOverrides(overrides, overridden);
  if (overridden.input_root_path != "/tmp/override_input.root" ||
      overridden.output_root_path != "/tmp/override_output.root" ||
      overridden.input_schema != femto3d::InputSchema::kLegacyVectorTree) {
    std::cerr << "Expected CLI overrides to replace configured run inputs.\n";
    return 13;
  }
  if (!overridden.analysis.event_plane.use_internal_reconstruction ||
      overridden.analysis.event_plane.fallback_to_input_branch) {
    std::cerr << "Expected schema CLI override to update default event-plane source.\n";
    return 14;
  }
  if (overridden.analysis.top_directory_name != valid.analysis.top_directory_name ||
      overridden.analysis.r2_summary_directory_name !=
          valid.analysis.r2_summary_directory_name) {
    std::cerr << "Expected CLI overrides to preserve unrelated config fields.\n";
    return 15;
  }

  const std::string strict_summary_policy = R"toml(
[input]
schema = "blastwave_flat_trees"
input_root = "input.root"

[output]
output_root = "output.root"

[projection_fit]
accept_hbt_central_value_only_for_summary = false
)toml";

  const femto3d::ApplicationConfig strict_summary_config =
      femto3d::LoadApplicationConfig(
          WriteFile(temp_dir / "strict_summary_policy.toml",
                    strict_summary_policy));
  if (strict_summary_config.analysis.projection_fit
          .accept_hbt_central_value_only_for_summary) {
    std::cerr << "Expected explicit HBT summary policy override to disable "
                 "central-value-only fallback.\n";
    return 17;
  }

  const std::filesystem::path project_root =
      std::filesystem::path(__FILE__).parent_path().parent_path();
  (void)femto3d::LoadApplicationConfig(
      (project_root / "config/examples/legacy_vector_tree.example.toml").string());
  (void)femto3d::LoadApplicationConfig(
      (project_root / "config/examples/blastwave_flat_trees.example.toml").string());

  return 0;
}
