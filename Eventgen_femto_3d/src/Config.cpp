#include "femto3d/Config.h"

#include <toml++/toml.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace femto3d {

namespace {

std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

const toml::table& GetRequiredTable(const toml::table& table,
                                    const std::string& key,
                                    const std::string& context) {
  const auto* child = table[key].as_table();
  if (child == nullptr) {
    throw ConfigError("Missing required table '" + key + "' in " + context);
  }
  return *child;
}

const toml::table* GetOptionalTable(const toml::table& table,
                                    const std::string& key) {
  return table[key].as_table();
}

const toml::array* GetOptionalArray(const toml::table& table,
                                    const std::string& key) {
  return table[key].as_array();
}

std::string ReadRequiredString(const toml::table& table,
                               const std::string& key,
                               const std::string& context) {
  if (const auto value = table[key].value<std::string>(); value.has_value()) {
    return *value;
  }
  throw ConfigError("Missing required string '" + key + "' in " + context);
}

std::string ReadOptionalString(const toml::table& table,
                               const std::string& key,
                               const std::string& fallback) {
  if (const auto value = table[key].value<std::string>(); value.has_value()) {
    return *value;
  }
  return fallback;
}

bool ReadOptionalBool(const toml::table& table,
                      const std::string& key,
                      const bool fallback) {
  if (const auto value = table[key].value<bool>(); value.has_value()) {
    return *value;
  }
  return fallback;
}

double ReadOptionalDouble(const toml::table& table,
                          const std::string& key,
                          const double fallback) {
  if (const auto value = table[key].value<double>(); value.has_value()) {
    return *value;
  }
  return fallback;
}

int ReadOptionalInt(const toml::table& table,
                    const std::string& key,
                    const int fallback) {
  if (const auto value = table[key].value<int64_t>(); value.has_value()) {
    return static_cast<int>(*value);
  }
  return fallback;
}

std::size_t ReadOptionalSize(const toml::table& table,
                             const std::string& key,
                             const std::size_t fallback) {
  if (const auto value = table[key].value<int64_t>(); value.has_value()) {
    if (*value < 0) {
      throw ConfigError("Expected non-negative integer for key '" + key + "'.");
    }
    return static_cast<std::size_t>(*value);
  }
  return fallback;
}

std::vector<int> ReadOptionalIntArray(const toml::table& table,
                                      const std::string& key,
                                      const std::vector<int>& fallback) {
  const toml::array* array = table[key].as_array();
  if (array == nullptr) {
    return fallback;
  }

  std::vector<int> values;
  values.reserve(array->size());
  for (std::size_t index = 0; index < array->size(); ++index) {
    const toml::node* node = array->get(index);
    if (node == nullptr || !node->is_integer()) {
      throw ConfigError("Expected integer entry in array '" + key + "'.");
    }
    values.push_back(static_cast<int>(node->value_or<int64_t>(0)));
  }
  return values;
}

RangeBin ParseRangeBin(const toml::table& table, const std::string& context) {
  RangeBin bin;
  if (const auto value = table["min"].value<double>(); value.has_value()) {
    bin.min = *value;
  } else if (const auto value = table["low"].value<double>(); value.has_value()) {
    bin.min = *value;
  } else {
    throw ConfigError("Missing min/low in " + context);
  }

  if (const auto value = table["max"].value<double>(); value.has_value()) {
    bin.max = *value;
  } else if (const auto value = table["high"].value<double>(); value.has_value()) {
    bin.max = *value;
  } else {
    throw ConfigError("Missing max/high in " + context);
  }

  if (const auto label = table["label"].value<std::string>(); label.has_value()) {
    bin.label = *label;
  }
  return bin;
}

std::vector<RangeBin> ParseRangeBinArray(const toml::array* array,
                                         const std::string& context) {
  std::vector<RangeBin> bins;
  if (array == nullptr) {
    return bins;
  }

  for (std::size_t index = 0; index < array->size(); ++index) {
    const toml::node* node = array->get(index);
    const auto* table = node != nullptr ? node->as_table() : nullptr;
    if (table == nullptr) {
      throw ConfigError("Expected table entry in " + context);
    }
    bins.push_back(ParseRangeBin(*table, context));
  }

  return bins;
}

AxisSpec ParseAxisSpec(const toml::table& table,
                       const AxisSpec& fallback,
                       const std::string& context) {
  AxisSpec axis = fallback;
  axis.title = ReadOptionalString(table, "title", axis.title);
  axis.min = ReadOptionalDouble(table, "min", axis.min);
  axis.max = ReadOptionalDouble(table, "max", axis.max);
  axis.bin_width = ReadOptionalDouble(table, "bin_width", axis.bin_width);
  if (axis.title.empty()) {
    throw ConfigError(context + ".title must not be empty.");
  }
  return axis;
}

bool MatchesRangeBin(const RangeBin& left,
                     const RangeBin& right,
                     const double tolerance = 1.0e-12) {
  return std::abs(left.min - right.min) < tolerance &&
         std::abs(left.max - right.max) < tolerance;
}

bool IsValidRangeBin(const RangeBin& bin) {
  return std::isfinite(bin.min) && std::isfinite(bin.max) && bin.max > bin.min;
}

void ValidateRangeCollection(const std::string& label,
                             std::vector<RangeBin>& bins) {
  if (bins.empty()) {
    throw ConfigError("Missing required bin collection: " + label);
  }

  for (std::size_t index = 0; index < bins.size(); ++index) {
    RangeBin& bin = bins[index];
    if (!IsValidRangeBin(bin)) {
      throw ConfigError("Invalid range bin in " + label + " at index " +
                        std::to_string(index));
    }
    if (bin.label.empty()) {
      bin.label = FormatBinLabel(label, bin.min, bin.max);
    }
  }

  for (std::size_t i = 0; i < bins.size(); ++i) {
    for (std::size_t j = i + 1U; j < bins.size(); ++j) {
      if (MatchesRangeBin(bins[i], bins[j])) {
        throw ConfigError("Duplicate range bin found in " + label);
      }
      const bool overlap = bins[i].min < bins[j].max && bins[j].min < bins[i].max;
      if (overlap) {
        throw ConfigError("Overlapping range bins found in " + label);
      }
    }
  }
}

void ValidateAxis(const std::string& label, const AxisSpec& axis) {
  if (axis.title.empty()) {
    throw ConfigError(label + ".title must not be empty.");
  }
  if (!std::isfinite(axis.min) || !std::isfinite(axis.max) ||
      axis.max <= axis.min) {
    throw ConfigError(label + " requires finite min/max with max > min.");
  }
  if (!std::isfinite(axis.bin_width) || axis.bin_width <= 0.0) {
    throw ConfigError(label + ".bin_width must be finite and positive.");
  }
}

void ValidateNonEmpty(const std::string& label, const std::string& value) {
  if (value.empty()) {
    throw ConfigError(label + " must not be empty.");
  }
}

void ValidateLegacyInputConfig(const LegacyInputTreeConfig& config) {
  ValidateNonEmpty("input.legacy.tree_name", config.tree_name);
  ValidateNonEmpty("input.legacy.centrality_branch", config.centrality_branch);
  ValidateNonEmpty("input.legacy.event_plane_branch", config.event_plane_branch);
  ValidateNonEmpty("input.legacy.pdg_branch", config.pdg_branch);
  ValidateNonEmpty("input.legacy.px_branch", config.px_branch);
  ValidateNonEmpty("input.legacy.py_branch", config.py_branch);
  ValidateNonEmpty("input.legacy.pz_branch", config.pz_branch);
  ValidateNonEmpty("input.legacy.mass_branch", config.mass_branch);
  ValidateNonEmpty("input.legacy.x_branch", config.x_branch);
  ValidateNonEmpty("input.legacy.y_branch", config.y_branch);
  ValidateNonEmpty("input.legacy.z_branch", config.z_branch);
  ValidateNonEmpty("input.legacy.t_branch", config.t_branch);
}

void ValidateBlastwaveInputConfig(const BlastwaveInputTreeConfig& config) {
  ValidateNonEmpty("input.blastwave.events_tree", config.events_tree);
  ValidateNonEmpty("input.blastwave.particles_tree", config.particles_tree);
  ValidateNonEmpty("input.blastwave.event_id_branch", config.event_id_branch);
  ValidateNonEmpty("input.blastwave.centrality_branch", config.centrality_branch);
  ValidateNonEmpty("input.blastwave.event_plane_branch", config.event_plane_branch);
  ValidateNonEmpty("input.blastwave.pid_branch", config.pid_branch);
  ValidateNonEmpty("input.blastwave.px_branch", config.px_branch);
  ValidateNonEmpty("input.blastwave.py_branch", config.py_branch);
  ValidateNonEmpty("input.blastwave.pz_branch", config.pz_branch);
  ValidateNonEmpty("input.blastwave.mass_branch", config.mass_branch);
  ValidateNonEmpty("input.blastwave.x_branch", config.x_branch);
  ValidateNonEmpty("input.blastwave.y_branch", config.y_branch);
  ValidateNonEmpty("input.blastwave.z_branch", config.z_branch);
  ValidateNonEmpty("input.blastwave.t_branch", config.t_branch);
}

void ParseLegacyInputConfig(const toml::table* table,
                            LegacyInputTreeConfig& config) {
  if (table == nullptr) {
    return;
  }
  config.tree_name = ReadOptionalString(*table, "tree_name", config.tree_name);
  config.centrality_branch =
      ReadOptionalString(*table, "centrality_branch", config.centrality_branch);
  config.event_plane_branch =
      ReadOptionalString(*table, "event_plane_branch", config.event_plane_branch);
  config.pdg_branch = ReadOptionalString(*table, "pdg_branch", config.pdg_branch);
  config.px_branch = ReadOptionalString(*table, "px_branch", config.px_branch);
  config.py_branch = ReadOptionalString(*table, "py_branch", config.py_branch);
  config.pz_branch = ReadOptionalString(*table, "pz_branch", config.pz_branch);
  config.mass_branch =
      ReadOptionalString(*table, "mass_branch", config.mass_branch);
  config.x_branch = ReadOptionalString(*table, "x_branch", config.x_branch);
  config.y_branch = ReadOptionalString(*table, "y_branch", config.y_branch);
  config.z_branch = ReadOptionalString(*table, "z_branch", config.z_branch);
  config.t_branch = ReadOptionalString(*table, "t_branch", config.t_branch);
}

void ParseBlastwaveInputConfig(const toml::table* table,
                               BlastwaveInputTreeConfig& config) {
  if (table == nullptr) {
    return;
  }
  config.events_tree = ReadOptionalString(*table, "events_tree", config.events_tree);
  config.particles_tree =
      ReadOptionalString(*table, "particles_tree", config.particles_tree);
  config.event_id_branch =
      ReadOptionalString(*table, "event_id_branch", config.event_id_branch);
  config.centrality_branch =
      ReadOptionalString(*table, "centrality_branch", config.centrality_branch);
  config.event_plane_branch =
      ReadOptionalString(*table, "event_plane_branch", config.event_plane_branch);
  config.pid_branch = ReadOptionalString(
      *table,
      "pdg_branch",
      ReadOptionalString(*table, "pid_branch", config.pid_branch));
  config.px_branch = ReadOptionalString(*table, "px_branch", config.px_branch);
  config.py_branch = ReadOptionalString(*table, "py_branch", config.py_branch);
  config.pz_branch = ReadOptionalString(*table, "pz_branch", config.pz_branch);
  config.mass_branch =
      ReadOptionalString(*table, "mass_branch", config.mass_branch);
  config.x_branch = ReadOptionalString(*table, "x_branch", config.x_branch);
  config.y_branch = ReadOptionalString(*table, "y_branch", config.y_branch);
  config.z_branch = ReadOptionalString(*table, "z_branch", config.z_branch);
  config.t_branch = ReadOptionalString(*table, "t_branch", config.t_branch);
}

}  // namespace

CliOptions ParseCliArgs(const int argc, char** argv) {
  CliOptions options;

  for (int index = 1; index < argc; ++index) {
    const std::string token = argv[index];
    if (token == "--config") {
      if (index + 1 >= argc) {
        throw std::runtime_error("Missing value after --config.");
      }
      options.config_path = argv[++index];
      continue;
    }
    if (token == "--input-root") {
      if (index + 1 >= argc) {
        throw std::runtime_error("Missing value after --input-root.");
      }
      options.input_root_override = std::string(argv[++index]);
      continue;
    }
    if (token == "--output-root") {
      if (index + 1 >= argc) {
        throw std::runtime_error("Missing value after --output-root.");
      }
      options.output_root_override = std::string(argv[++index]);
      continue;
    }
    if (token == "--input-schema") {
      if (index + 1 >= argc) {
        throw std::runtime_error("Missing value after --input-schema.");
      }
      options.input_schema_override = ParseInputSchema(argv[++index]);
      continue;
    }
    throw std::runtime_error("Unknown argument: " + token);
  }

  if (options.config_path.empty()) {
    throw std::runtime_error("--config is required.");
  }

  return options;
}

void PrintUsage(std::ostream& stream) {
  stream
      << "Usage:\n"
      << "  eventgen_femto_3d --config <file.toml> "
         "[--input-root <path>] [--output-root <path>] "
         "[--input-schema legacy_vector_tree|blastwave_flat_trees]\n";
}

ApplicationConfig LoadApplicationConfig(const std::string& path) {
  ApplicationConfig config = MakeDefaultApplicationConfig();

  toml::table root;
  try {
    root = toml::parse_file(path);
  } catch (const toml::parse_error& error) {
    throw ConfigError(std::string("TOML parse error: ") +
                      std::string(error.description()));
  }

  const toml::table& input = GetRequiredTable(root, "input", "root");
  const toml::table& output = GetRequiredTable(root, "output", "root");

  config.input_schema = ParseInputSchema(
      ReadRequiredString(input, "schema", "input"));
  config.analysis.event_plane =
      MakeDefaultEventPlaneConfigForSchema(config.input_schema);
  config.input_root_path = ReadRequiredString(input, "input_root", "input");
  ParseLegacyInputConfig(GetOptionalTable(input, "legacy"), config.analysis.input.legacy);
  ParseBlastwaveInputConfig(GetOptionalTable(input, "blastwave"),
                            config.analysis.input.blastwave);

  config.output_root_path =
      ReadRequiredString(output, "output_root", "output");
  config.analysis.top_directory_name = ReadOptionalString(
      output, "top_directory_name", config.analysis.top_directory_name);
  config.analysis.r2_summary_directory_name = ReadOptionalString(
      output,
      "r2_summary_directory_name",
      config.analysis.r2_summary_directory_name);

  if (const toml::table* event_plane = root["event_plane"].as_table();
      event_plane != nullptr) {
    config.analysis.event_plane.enabled =
        ReadOptionalBool(*event_plane, "enabled", config.analysis.event_plane.enabled);
    config.analysis.event_plane.harmonic_order = ReadOptionalInt(
        *event_plane,
        "harmonic_order",
        config.analysis.event_plane.harmonic_order);
    config.analysis.event_plane.use_internal_reconstruction =
        ReadOptionalBool(*event_plane,
                         "use_internal_reconstruction",
                         config.analysis.event_plane.use_internal_reconstruction);
    config.analysis.event_plane.fallback_to_input_branch = ReadOptionalBool(
        *event_plane,
        "fallback_to_input_branch",
        config.analysis.event_plane.fallback_to_input_branch);
    config.analysis.event_plane.eta_min =
        ReadOptionalDouble(*event_plane, "eta_min", config.analysis.event_plane.eta_min);
    config.analysis.event_plane.eta_max =
        ReadOptionalDouble(*event_plane, "eta_max", config.analysis.event_plane.eta_max);
    config.analysis.event_plane.weight_mode = ParseEventPlaneWeightMode(
        ReadOptionalString(*event_plane,
                           "weight_mode",
                           ToString(config.analysis.event_plane.weight_mode)));
    config.analysis.event_plane.allowed_abs_pdg = ReadOptionalIntArray(
        *event_plane,
        "allowed_abs_pdg",
        config.analysis.event_plane.allowed_abs_pdg);
    config.analysis.event_plane.min_candidates = ReadOptionalSize(
        *event_plane,
        "min_candidates",
        config.analysis.event_plane.min_candidates);
    config.analysis.event_plane.min_q_magnitude = ReadOptionalDouble(
        *event_plane,
        "min_q_magnitude",
        config.analysis.event_plane.min_q_magnitude);
  }

  if (const toml::table* selection = root["selection"].as_table();
      selection != nullptr) {
    config.analysis.selection.target_pdg =
        ReadOptionalInt(*selection, "target_pdg", config.analysis.selection.target_pdg);
    config.analysis.selection.femto_eta_min = ReadOptionalDouble(
        *selection,
        "femto_eta_min",
        config.analysis.selection.femto_eta_min);
    config.analysis.selection.femto_eta_max = ReadOptionalDouble(
        *selection,
        "femto_eta_max",
        config.analysis.selection.femto_eta_max);
    config.analysis.selection.femto_pt_min = ReadOptionalDouble(
        *selection,
        "femto_pt_min",
        config.analysis.selection.femto_pt_min);
    config.analysis.selection.femto_pt_max = ReadOptionalDouble(
        *selection,
        "femto_pt_max",
        config.analysis.selection.femto_pt_max);
    config.analysis.selection.femto_mt_reference_mass = ReadOptionalDouble(
        *selection,
        "femto_mt_reference_mass",
        config.analysis.selection.femto_mt_reference_mass);
    config.analysis.selection.femto_mt_mass_tolerance = ReadOptionalDouble(
        *selection,
        "femto_mt_mass_tolerance",
        config.analysis.selection.femto_mt_mass_tolerance);
  }

  if (const toml::table* histograms = root["histograms"].as_table();
      histograms != nullptr) {
    config.analysis.histograms.warn_on_overflow =
        ReadOptionalBool(*histograms,
                         "warn_on_overflow",
                         config.analysis.histograms.warn_on_overflow);
    if (const toml::table* axis = GetOptionalTable(*histograms, "rho_out_axis");
        axis != nullptr) {
      config.analysis.histograms.rho_out_axis = ParseAxisSpec(
          *axis, config.analysis.histograms.rho_out_axis, "histograms.rho_out_axis");
    }
    if (const toml::table* axis = GetOptionalTable(*histograms, "rho_side_axis");
        axis != nullptr) {
      config.analysis.histograms.rho_side_axis = ParseAxisSpec(
          *axis, config.analysis.histograms.rho_side_axis, "histograms.rho_side_axis");
    }
    if (const toml::table* axis = GetOptionalTable(*histograms, "rho_long_axis");
        axis != nullptr) {
      config.analysis.histograms.rho_long_axis = ParseAxisSpec(
          *axis, config.analysis.histograms.rho_long_axis, "histograms.rho_long_axis");
    }
    if (const toml::table* axis = GetOptionalTable(*histograms, "projection_axis");
        axis != nullptr) {
      config.analysis.histograms.projection_axis = ParseAxisSpec(
          *axis, config.analysis.histograms.projection_axis, "histograms.projection_axis");
    }
  }

  if (const toml::table* projection_fit = root["projection_fit"].as_table();
      projection_fit != nullptr) {
    config.analysis.projection_fit.use_adaptive_integration =
        ReadOptionalBool(*projection_fit,
                         "use_adaptive_integration",
                         config.analysis.projection_fit.use_adaptive_integration);
    config.analysis.projection_fit.accept_forced_posdef_covariance_as_valid =
        ReadOptionalBool(
            *projection_fit,
            "accept_forced_posdef_covariance_as_valid",
            config.analysis.projection_fit
                .accept_forced_posdef_covariance_as_valid);
    config.analysis.projection_fit.fail_alpha_result_when_error_invalid =
        ReadOptionalBool(
            *projection_fit,
            "fail_alpha_result_when_error_invalid",
            config.analysis.projection_fit.fail_alpha_result_when_error_invalid);
    config.analysis.projection_fit.fail_hbt_results_when_error_invalid =
        ReadOptionalBool(
            *projection_fit,
            "fail_hbt_results_when_error_invalid",
            config.analysis.projection_fit.fail_hbt_results_when_error_invalid);
    config.analysis.projection_fit.fail_directional_results_when_error_invalid =
        ReadOptionalBool(
            *projection_fit,
            "fail_directional_results_when_error_invalid",
            config.analysis.projection_fit
                .fail_directional_results_when_error_invalid);
    config.analysis.projection_fit.accept_hbt_central_value_only_for_summary =
        ReadOptionalBool(
            *projection_fit,
            "accept_hbt_central_value_only_for_summary",
            config.analysis.projection_fit
                .accept_hbt_central_value_only_for_summary);
  }

  if (const toml::table* bins = root["bins"].as_table(); bins != nullptr) {
    if (const toml::array* array = GetOptionalArray(*bins, "centrality");
        array != nullptr) {
      config.analysis.centrality_bins =
          ParseRangeBinArray(array, "bins.centrality");
    }

    if (const toml::array* array = GetOptionalArray(*bins, "mt");
        array != nullptr) {
      config.analysis.mt_bins = ParseRangeBinArray(array, "bins.mt");
    }

    if (const toml::array* array = GetOptionalArray(*bins, "phi");
        array != nullptr) {
      config.analysis.phi_bins = ParseRangeBinArray(array, "bins.phi");
    }
  }

  ValidateApplicationConfig(config);
  return config;
}

void ValidateApplicationConfig(ApplicationConfig& config) {
  ValidateNonEmpty("input.input_root", config.input_root_path);
  ValidateNonEmpty("output.output_root", config.output_root_path);
  ValidateNonEmpty("output.top_directory_name", config.analysis.top_directory_name);
  ValidateNonEmpty("output.r2_summary_directory_name",
                   config.analysis.r2_summary_directory_name);

  ValidateLegacyInputConfig(config.analysis.input.legacy);
  ValidateBlastwaveInputConfig(config.analysis.input.blastwave);

  if (config.analysis.selection.femto_eta_max <=
      config.analysis.selection.femto_eta_min) {
    throw ConfigError("selection.femto_eta_max must be greater than femto_eta_min.");
  }
  if (config.analysis.selection.femto_pt_max <=
      config.analysis.selection.femto_pt_min) {
    throw ConfigError("selection.femto_pt_max must be greater than femto_pt_min.");
  }
  if (!std::isfinite(config.analysis.selection.femto_mt_reference_mass) ||
      config.analysis.selection.femto_mt_reference_mass < 0.0) {
    throw ConfigError("selection.femto_mt_reference_mass must be finite and non-negative.");
  }
  if (!std::isfinite(config.analysis.selection.femto_mt_mass_tolerance) ||
      config.analysis.selection.femto_mt_mass_tolerance < 0.0) {
    throw ConfigError("selection.femto_mt_mass_tolerance must be finite and non-negative.");
  }
  if (config.analysis.event_plane.harmonic_order <= 0) {
    throw ConfigError("event_plane.harmonic_order must be positive.");
  }
  if (config.analysis.event_plane.eta_max <= config.analysis.event_plane.eta_min) {
    throw ConfigError("event_plane.eta_max must be greater than eta_min.");
  }
  if (!std::isfinite(config.analysis.event_plane.min_q_magnitude) ||
      config.analysis.event_plane.min_q_magnitude < 0.0) {
    throw ConfigError("event_plane.min_q_magnitude must be finite and non-negative.");
  }

  ValidateAxis("histograms.rho_out_axis", config.analysis.histograms.rho_out_axis);
  ValidateAxis("histograms.rho_side_axis", config.analysis.histograms.rho_side_axis);
  ValidateAxis("histograms.rho_long_axis", config.analysis.histograms.rho_long_axis);
  ValidateAxis("histograms.projection_axis",
               config.analysis.histograms.projection_axis);
  ValidateRangeCollection("bins.centrality", config.analysis.centrality_bins);
  ValidateRangeCollection("bins.mt", config.analysis.mt_bins);
  ValidateRangeCollection("bins.phi", config.analysis.phi_bins);
}

void ApplyCliOverrides(const CliOptions& cli_options, ApplicationConfig& config) {
  if (cli_options.input_root_override.has_value()) {
    config.input_root_path = *cli_options.input_root_override;
  }
  if (cli_options.output_root_override.has_value()) {
    config.output_root_path = *cli_options.output_root_override;
  }
  if (cli_options.input_schema_override.has_value()) {
    const EventPlaneConfig current_default =
        MakeDefaultEventPlaneConfigForSchema(config.input_schema);
    if (config.analysis.event_plane.use_internal_reconstruction ==
            current_default.use_internal_reconstruction &&
        config.analysis.event_plane.fallback_to_input_branch ==
            current_default.fallback_to_input_branch) {
      config.analysis.event_plane =
          MakeDefaultEventPlaneConfigForSchema(*cli_options.input_schema_override);
    }
    config.input_schema = *cli_options.input_schema_override;
  }
  ValidateApplicationConfig(config);
}

std::string ToString(const InputSchema schema) {
  switch (schema) {
    case InputSchema::kLegacyVectorTree:
      return "legacy_vector_tree";
    case InputSchema::kBlastwaveFlatTrees:
      return "blastwave_flat_trees";
  }
  return "blastwave_flat_trees";
}

InputSchema ParseInputSchema(const std::string& token) {
  const std::string lowered = ToLower(token);
  if (lowered == "legacy_vector_tree") {
    return InputSchema::kLegacyVectorTree;
  }
  if (lowered == "blastwave_flat_trees") {
    return InputSchema::kBlastwaveFlatTrees;
  }
  throw ConfigError("Unsupported input schema: " + token);
}

std::string ToString(const EventPlaneWeightMode mode) {
  switch (mode) {
    case EventPlaneWeightMode::kUnit:
      return "unit";
    case EventPlaneWeightMode::kPt:
      return "pt";
  }
  return "pt";
}

EventPlaneWeightMode ParseEventPlaneWeightMode(const std::string& token) {
  const std::string lowered = ToLower(token);
  if (lowered == "unit") {
    return EventPlaneWeightMode::kUnit;
  }
  if (lowered == "pt") {
    return EventPlaneWeightMode::kPt;
  }
  throw ConfigError("Unsupported event-plane weight mode: " + token);
}

}  // namespace femto3d
