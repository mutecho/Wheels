#include "exp_femto_1d/Config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

namespace exp_femto_1d {

  namespace {

    std::string ToLower(std::string value) {
      std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
      });
      return value;
    }

    bool ContainsExactRange(const std::vector<RangeBin> &bins, const RangeBin &bin) {
      return std::any_of(bins.begin(), bins.end(), [&](const RangeBin &candidate) {
        return MatchesRangeBin(candidate, bin);
      });
    }

    // Validate bin collections while preserving the 3D project's overlap policy.
    void ValidateRangeCollection(const std::string &label,
                                 std::vector<RangeBin> &bins,
                                 const bool allow_overlap = false) {
      if (bins.empty()) {
        throw ConfigError("Missing required bin collection: " + label);
      }

      for (std::size_t index = 0; index < bins.size(); ++index) {
        RangeBin &bin = bins[index];
        if (!IsValidRangeBin(bin)) {
          throw ConfigError("Invalid range bin in " + label + " at index " + std::to_string(index));
        }
        if (bin.label.empty()) {
          std::ostringstream label_stream;
          label_stream << label << "_" << bin.min << "_" << bin.max;
          bin.label = label_stream.str();
        }
      }

      for (std::size_t i = 0; i < bins.size(); ++i) {
        for (std::size_t j = i + 1U; j < bins.size(); ++j) {
          if (MatchesRangeBin(bins[i], bins[j])) {
            throw ConfigError("Duplicate range bin found in " + label);
          }
          const bool overlap = bins[i].min < bins[j].max && bins[j].min < bins[i].max;
          if (overlap && !allow_overlap) {
            throw ConfigError("Overlapping range bins found in " + label);
          }
        }
      }
    }

    // Normalize user-facing output names to the fixed extensions used by the contract.
    void EnsureExtension(std::string &value, const std::string &extension) {
      if (value.empty()) {
        return;
      }
      if (value.size() >= extension.size()
          && value.compare(value.size() - extension.size(), extension.size(), extension) == 0) {
        return;
      }
      value += extension;
    }

    const toml::table &GetRequiredTable(const toml::table &table, const std::string &key, const std::string &context) {
      const auto *child = table[key].as_table();
      if (child == nullptr) {
        throw ConfigError("Missing required table '" + key + "' in " + context);
      }
      return *child;
    }

    const toml::array *GetOptionalArray(const toml::table &table, const std::string &key) {
      return table[key].as_array();
    }

    std::string ReadRequiredString(const toml::table &table, const std::string &key, const std::string &context) {
      if (const auto value = table[key].value<std::string>(); value.has_value()) {
        return *value;
      }
      throw ConfigError("Missing required string '" + key + "' in " + context);
    }

    std::string ReadOptionalString(const toml::table &table, const std::string &key, const std::string &fallback) {
      if (const auto value = table[key].value<std::string>(); value.has_value()) {
        return *value;
      }
      return fallback;
    }

    bool ReadOptionalBool(const toml::table &table, const std::string &key, const bool fallback) {
      if (const auto value = table[key].value<bool>(); value.has_value()) {
        return *value;
      }
      return fallback;
    }

    double ReadOptionalDouble(const toml::table &table, const std::string &key, const double fallback) {
      if (const auto value = table[key].value<double>(); value.has_value()) {
        return *value;
      }
      return fallback;
    }

    unsigned ReadOptionalUnsigned(const toml::table &table, const std::string &key, const unsigned fallback) {
      if (const auto value = table[key].value<int64_t>(); value.has_value()) {
        if (*value < 0) {
          throw ConfigError("Expected non-negative integer for key '" + key + "'.");
        }
        return static_cast<unsigned>(*value);
      }
      return fallback;
    }

    ProgressMode ReadOptionalProgressMode(const toml::table &table,
                                          const std::string &key,
                                          const ProgressMode fallback) {
      if (const auto value = table[key].value<bool>(); value.has_value()) {
        return *value ? ProgressMode::kEnabled : ProgressMode::kDisabled;
      }
      if (const auto value = table[key].value<std::string>(); value.has_value()) {
        return ParseProgressMode(*value);
      }
      return fallback;
    }

    RangeBin ParseRangeBin(const toml::table &table, const std::string &context) {
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

    std::vector<RangeBin> ParseRangeBinArray(const toml::array *array, const std::string &context) {
      std::vector<RangeBin> bins;
      if (array == nullptr) {
        return bins;
      }
      for (std::size_t index = 0; index < array->size(); ++index) {
        const toml::node *node = array->get(index);
        const auto *table = node != nullptr ? node->as_table() : nullptr;
        if (table == nullptr) {
          throw ConfigError("Expected table entry in " + context);
        }
        bins.push_back(ParseRangeBin(*table, context));
      }
      return bins;
    }

    void ValidateFinitePositive(const double value, const std::string &name) {
      if (!std::isfinite(value) || value <= 0.0) {
        throw ConfigError(name + " must be finite and positive.");
      }
    }

    void ValidateFiniteIncreasing(const double low, const double high, const std::string &name) {
      if (!std::isfinite(low) || !std::isfinite(high) || !(low < high)) {
        throw ConfigError(name + " must satisfy low < high.");
      }
    }

    void ValidateInitInRange(const double value, const double low, const double high, const std::string &name) {
      if (!std::isfinite(value) || value < low || value > high) {
        throw ConfigError(name + " init value must lie within its configured range.");
      }
    }

  }  // namespace

  ApplicationConfig LoadApplicationConfig(const std::string &path) {
    ApplicationConfig config;
    toml::table root;
    try {
      root = toml::parse_file(path);
    } catch (const toml::parse_error &error) {
      throw ConfigError(std::string("TOML parse error: ") + std::string(error.description()));
    }

    const toml::table &input = GetRequiredTable(root, "input", "root");
    const toml::table &output = GetRequiredTable(root, "output", "root");
    const toml::table &build = GetRequiredTable(root, "build", "root");
    const toml::table &fit = GetRequiredTable(root, "fit", "root");

    config.input.input_root = ReadRequiredString(input, "input_root", "input");
    config.input.task_name = ReadRequiredString(input, "task_name", "input");
    config.input.same_event_subtask = ReadRequiredString(input, "same_event_subtask", "input");
    config.input.mixed_event_subtask = ReadRequiredString(input, "mixed_event_subtask", "input");
    config.input.sparse_object_name = ReadRequiredString(input, "sparse_object_name", "input");

    config.output.output_directory = ReadRequiredString(output, "output_directory", "output");
    config.output.cf_root_name = ReadOptionalString(output, "cf_root_name", config.output.cf_root_name);
    config.output.fit_root_name = ReadOptionalString(output, "fit_root_name", config.output.fit_root_name);
    config.output.fit_summary_name = ReadOptionalString(output, "fit_summary_name", config.output.fit_summary_name);
    config.output.log_level = ParseLogLevel(ReadOptionalString(output, "log_level", ToString(config.output.log_level)));

    config.build.norm_low = ReadOptionalDouble(build, "norm_low", config.build.norm_low);
    config.build.norm_high = ReadOptionalDouble(build, "norm_high", config.build.norm_high);
    config.build.kstar_min = ReadOptionalDouble(build, "kstar_min", config.build.kstar_min);
    config.build.kstar_max = ReadOptionalDouble(build, "kstar_max", config.build.kstar_max);
    config.build.reopen_output_file_per_slice =
        ReadOptionalBool(build, "reopen_output_file_per_slice", config.build.reopen_output_file_per_slice);
    config.build.progress = ReadOptionalProgressMode(build, "progress", config.build.progress);

    config.fit.fit_kstar_max = ReadOptionalDouble(fit, "fit_kstar_max", config.fit.fit_kstar_max);
    config.fit.use_coulomb = ReadOptionalBool(fit, "use_coulomb", config.fit.use_coulomb);
    config.fit.reopen_output_file_per_slice =
        ReadOptionalBool(fit, "reopen_output_file_per_slice", config.fit.reopen_output_file_per_slice);
    config.fit.progress = ReadOptionalProgressMode(fit, "progress", config.fit.progress);

    config.fit.baseline_p0_init = ReadOptionalDouble(fit, "baseline_p0_init", config.fit.baseline_p0_init);
    config.fit.baseline_p0_min = ReadOptionalDouble(fit, "baseline_p0_min", config.fit.baseline_p0_min);
    config.fit.baseline_p0_max = ReadOptionalDouble(fit, "baseline_p0_max", config.fit.baseline_p0_max);
    config.fit.baseline_p1_init = ReadOptionalDouble(fit, "baseline_p1_init", config.fit.baseline_p1_init);
    config.fit.baseline_p1_min = ReadOptionalDouble(fit, "baseline_p1_min", config.fit.baseline_p1_min);
    config.fit.baseline_p1_max = ReadOptionalDouble(fit, "baseline_p1_max", config.fit.baseline_p1_max);
    config.fit.baseline_p2_init = ReadOptionalDouble(fit, "baseline_p2_init", config.fit.baseline_p2_init);
    config.fit.baseline_p2_min = ReadOptionalDouble(fit, "baseline_p2_min", config.fit.baseline_p2_min);
    config.fit.baseline_p2_max = ReadOptionalDouble(fit, "baseline_p2_max", config.fit.baseline_p2_max);
    config.fit.baseline_p3_fixed = ReadOptionalBool(fit, "baseline_p3_fixed", config.fit.baseline_p3_fixed);
    config.fit.baseline_p3_value = ReadOptionalDouble(fit, "baseline_p3_value", config.fit.baseline_p3_value);
    config.fit.baseline_p4_fixed = ReadOptionalBool(fit, "baseline_p4_fixed", config.fit.baseline_p4_fixed);
    config.fit.baseline_p4_value = ReadOptionalDouble(fit, "baseline_p4_value", config.fit.baseline_p4_value);
    config.fit.source_size_init = ReadOptionalDouble(fit, "source_size_init", config.fit.source_size_init);
    config.fit.source_size_min = ReadOptionalDouble(fit, "source_size_min", config.fit.source_size_min);
    config.fit.source_size_max = ReadOptionalDouble(fit, "source_size_max", config.fit.source_size_max);
    config.fit.cats_num_mom_bins = ReadOptionalUnsigned(fit, "cats_num_mom_bins", config.fit.cats_num_mom_bins);
    config.fit.cats_kmin_mev = ReadOptionalDouble(fit, "cats_kmin_mev", config.fit.cats_kmin_mev);
    config.fit.cats_kmax_mev = ReadOptionalDouble(fit, "cats_kmax_mev", config.fit.cats_kmax_mev);

    if (const auto *bins = root["bins"].as_table(); bins != nullptr) {
      config.centrality_bins = ParseRangeBinArray(GetOptionalArray(*bins, "centrality"), "bins.centrality");
      config.mt_bins = ParseRangeBinArray(GetOptionalArray(*bins, "mt"), "bins.mt");
    }

    if (const auto *fit_selection = root["fit_selection"].as_table(); fit_selection != nullptr) {
      config.fit_centrality_bins =
          ParseRangeBinArray(GetOptionalArray(*fit_selection, "centrality"), "fit_selection.centrality");
      config.fit_mt_bins = ParseRangeBinArray(GetOptionalArray(*fit_selection, "mt"), "fit_selection.mt");
    }

    ValidateApplicationConfig(config);
    return config;
  }

  void ValidateApplicationConfig(ApplicationConfig &config) {
    if (config.input.input_root.empty()) {
      throw ConfigError("input.input_root is required.");
    }
    if (config.input.task_name.empty()) {
      throw ConfigError("input.task_name is required.");
    }
    if (config.input.same_event_subtask.empty()) {
      throw ConfigError("input.same_event_subtask is required.");
    }
    if (config.input.mixed_event_subtask.empty()) {
      throw ConfigError("input.mixed_event_subtask is required.");
    }
    if (config.input.sparse_object_name.empty()) {
      throw ConfigError("input.sparse_object_name is required.");
    }
    if (config.output.output_directory.empty()) {
      throw ConfigError("output.output_directory is required.");
    }

    ValidateFiniteIncreasing(config.build.norm_low, config.build.norm_high, "build normalization range");
    ValidateFiniteIncreasing(config.build.kstar_min, config.build.kstar_max, "build k* range");
    ValidateFinitePositive(config.fit.fit_kstar_max, "fit.fit_kstar_max");
    ValidateFiniteIncreasing(config.fit.baseline_p0_min, config.fit.baseline_p0_max, "fit baseline_p0 range");
    ValidateFiniteIncreasing(config.fit.baseline_p1_min, config.fit.baseline_p1_max, "fit baseline_p1 range");
    ValidateFiniteIncreasing(config.fit.baseline_p2_min, config.fit.baseline_p2_max, "fit baseline_p2 range");
    ValidateFiniteIncreasing(config.fit.source_size_min, config.fit.source_size_max, "fit source_size range");
    ValidateFinitePositive(config.fit.cats_kmax_mev - config.fit.cats_kmin_mev, "fit CATS momentum span");
    ValidateInitInRange(config.fit.baseline_p0_init,
                        config.fit.baseline_p0_min,
                        config.fit.baseline_p0_max,
                        "fit baseline_p0");
    ValidateInitInRange(config.fit.baseline_p1_init,
                        config.fit.baseline_p1_min,
                        config.fit.baseline_p1_max,
                        "fit baseline_p1");
    ValidateInitInRange(config.fit.baseline_p2_init,
                        config.fit.baseline_p2_min,
                        config.fit.baseline_p2_max,
                        "fit baseline_p2");
    ValidateInitInRange(config.fit.source_size_init,
                        config.fit.source_size_min,
                        config.fit.source_size_max,
                        "fit source_size");

    ValidateRangeCollection("bins.centrality", config.centrality_bins, true);
    ValidateRangeCollection("bins.mt", config.mt_bins, true);

    if (config.fit_centrality_bins.empty()) {
      config.fit_centrality_bins = config.centrality_bins;
    } else {
      ValidateRangeCollection("fit_selection.centrality", config.fit_centrality_bins, true);
      for (const RangeBin &bin : config.fit_centrality_bins) {
        if (!ContainsExactRange(config.centrality_bins, bin)) {
          throw ConfigError("fit_selection.centrality must match an existing build centrality bin exactly.");
        }
      }
    }

    if (config.fit_mt_bins.empty()) {
      config.fit_mt_bins = config.mt_bins;
    } else {
      ValidateRangeCollection("fit_selection.mt", config.fit_mt_bins, true);
      for (const RangeBin &bin : config.fit_mt_bins) {
        if (!ContainsExactRange(config.mt_bins, bin)) {
          throw ConfigError("fit_selection.mt must match an existing build mT bin exactly.");
        }
      }
    }

    EnsureExtension(config.output.cf_root_name, ".root");
    EnsureExtension(config.output.fit_root_name, ".root");
    EnsureExtension(config.output.fit_summary_name, ".tsv");
  }

  std::string ToString(const LogLevel level) {
    switch (level) {
      case LogLevel::kDebug:
        return "debug";
      case LogLevel::kInfo:
        return "info";
      case LogLevel::kWarn:
        return "warn";
      case LogLevel::kError:
        return "error";
    }
    return "info";
  }

  LogLevel ParseLogLevel(const std::string &token) {
    const std::string lowered = ToLower(token);
    if (lowered == "debug") {
      return LogLevel::kDebug;
    }
    if (lowered == "info") {
      return LogLevel::kInfo;
    }
    if (lowered == "warn" || lowered == "warning") {
      return LogLevel::kWarn;
    }
    if (lowered == "error") {
      return LogLevel::kError;
    }
    throw ConfigError("Unknown log level: " + token);
  }

  std::string ToString(const ProgressMode mode) {
    switch (mode) {
      case ProgressMode::kAuto:
        return "auto";
      case ProgressMode::kEnabled:
        return "enabled";
      case ProgressMode::kDisabled:
        return "disabled";
    }
    return "auto";
  }

  ProgressMode ParseProgressMode(const std::string &token) {
    const std::string lowered = ToLower(token);
    if (lowered == "auto") {
      return ProgressMode::kAuto;
    }
    if (lowered == "enabled" || lowered == "true") {
      return ProgressMode::kEnabled;
    }
    if (lowered == "disabled" || lowered == "false") {
      return ProgressMode::kDisabled;
    }
    throw ConfigError("Unknown progress mode: " + token);
  }

  std::string ToString(const RegionKind kind) {
    switch (kind) {
      case RegionKind::kMinBias:
        return "MinBias";
      case RegionKind::kInPlane:
        return "InPlane";
      case RegionKind::kOutOfPlane:
        return "OutOfPlane";
    }
    return "MinBias";
  }

}  // namespace exp_femto_1d
