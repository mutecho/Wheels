#include "exp_femto_3d/Config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

namespace exp_femto_3d {

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
          // Build and fit-selection lists may intentionally include merged ranges
          // that overlap smaller base bins, so only exact duplicates are rejected.
          const bool overlap = bins[i].min < bins[j].max && bins[j].min < bins[i].max;
          if (overlap && !allow_overlap) {
            throw ConfigError("Overlapping range bins found in " + label);
          }
        }
      }
    }

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

    std::optional<bool> ReadOptionalNullableBool(const toml::table &table, const std::string &key) {
      if (const auto value = table[key].value<bool>(); value.has_value()) {
        return *value;
      }
      return std::nullopt;
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

    double ReadOptionalDouble(const toml::table &table, const std::string &key, const double fallback) {
      if (const auto value = table[key].value<double>(); value.has_value()) {
        return *value;
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

    config.build.map_pair_phi_to_symmetric_range =
        ReadOptionalBool(build, "map_pair_phi_to_symmetric_range", config.build.map_pair_phi_to_symmetric_range);
    config.build.write_normalized_se_me_1d_projections = ReadOptionalBool(
        build, "write_normalized_se_me_1d_projections", config.build.write_normalized_se_me_1d_projections);
    config.build.reopen_output_file_per_slice =
        ReadOptionalBool(build, "reopen_output_file_per_slice", config.build.reopen_output_file_per_slice);
    config.build.progress = ReadOptionalProgressMode(build, "progress", config.build.progress);

    config.fit.model = ParseFitModel(ReadOptionalString(fit, "model", ToString(config.fit.model)));
    config.fit.options.use_coulomb = ReadOptionalBool(fit, "use_coulomb", config.fit.options.use_coulomb);
    config.fit.options.use_core_halo_lambda =
        ReadOptionalBool(fit, "use_core_halo_lambda", config.fit.options.use_core_halo_lambda);
    config.fit.options.use_q2_baseline = ReadOptionalBool(fit, "use_q2_baseline", config.fit.options.use_q2_baseline);
    config.fit.options.use_pml = ReadOptionalBool(fit, "use_pml", config.fit.options.use_pml);
    config.fit.options.fit_q_max = ReadOptionalDouble(fit, "fit_q_max", config.fit.options.fit_q_max);
    config.fit.map_pair_phi_to_symmetric_range = ReadOptionalNullableBool(fit, "map_pair_phi_to_symmetric_range");
    config.fit.reopen_output_file_per_slice =
        ReadOptionalBool(fit, "reopen_output_file_per_slice", config.fit.reopen_output_file_per_slice);
    config.fit.progress = ReadOptionalProgressMode(fit, "progress", config.fit.progress);

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
    if (config.fit.options.fit_q_max <= 0.0 || !std::isfinite(config.fit.options.fit_q_max)) {
      throw ConfigError("fit.fit_q_max must be finite and positive.");
    }

    ValidateRangeCollection("centrality", config.centrality_bins, true);
    ValidateRangeCollection("mt", config.mt_bins, true);

    if (config.fit_centrality_bins.empty()) {
      config.fit_centrality_bins = config.centrality_bins;
    }
    if (config.fit_mt_bins.empty()) {
      config.fit_mt_bins = config.mt_bins;
    }

    ValidateRangeCollection("fit_selection.centrality", config.fit_centrality_bins, true);
    ValidateRangeCollection("fit_selection.mt", config.fit_mt_bins, true);

    for (const RangeBin &bin : config.fit_centrality_bins) {
      if (!ContainsExactRange(config.centrality_bins, bin)) {
        throw ConfigError("Each fit_selection.centrality bin must exactly match a build bin.");
      }
    }
    for (const RangeBin &bin : config.fit_mt_bins) {
      if (!ContainsExactRange(config.mt_bins, bin)) {
        throw ConfigError("Each fit_selection.mt bin must exactly match a build bin.");
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
    throw ConfigError("Unsupported log level: " + token);
  }

  std::string ToString(const FitModel model) {
    switch (model) {
      case FitModel::kDiag:
        return "diag";
      case FitModel::kFull:
        return "full";
    }
    return "full";
  }

  FitModel ParseFitModel(const std::string &token) {
    const std::string lowered = ToLower(token);
    if (lowered == "diag") {
      return FitModel::kDiag;
    }
    if (lowered == "full") {
      return FitModel::kFull;
    }
    throw ConfigError("Unsupported fit model: " + token);
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
    if (lowered == "enabled" || lowered == "enable" || lowered == "on" || lowered == "true") {
      return ProgressMode::kEnabled;
    }
    if (lowered == "disabled" || lowered == "disable" || lowered == "off" || lowered == "false") {
      return ProgressMode::kDisabled;
    }
    throw ConfigError("Unsupported progress mode: " + token);
  }

}  // namespace exp_femto_3d
