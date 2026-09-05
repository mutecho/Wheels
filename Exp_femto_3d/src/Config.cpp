#include "exp_femto_3d/Config.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
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

    std::optional<double> ReadOptionalFiniteDouble(const toml::table &table,
                                                   const std::string &key,
                                                   const std::string &context) {
      if (!table.contains(key)) {
        return std::nullopt;
      }
      if (const auto value = table[key].value<double>(); value.has_value() && std::isfinite(*value)) {
        return *value;
      }
      throw ConfigError("Expected finite numeric field '" + key + "' in " + context + ".");
    }

    std::vector<std::string> ReadStringArray(const toml::table &table,
                                             const std::string &key,
                                             const std::string &context) {
      const auto *array = table[key].as_array();
      if (array == nullptr) {
        throw ConfigError("Expected string array '" + key + "' in " + context + ".");
      }
      std::vector<std::string> result;
      result.reserve(array->size());
      for (const toml::node &node : *array) {
        const auto value = node.value<std::string>();
        if (!value.has_value() || value->empty()) {
          throw ConfigError("Expected non-empty string values in " + context + "." + key + ".");
        }
        result.push_back(*value);
      }
      return result;
    }

    std::vector<int> ReadPointArray(const toml::table &table, const std::string &key, const std::string &context) {
      const auto *array = table[key].as_array();
      if (array == nullptr) {
        throw ConfigError("Expected integer array '" + key + "' in " + context + ".");
      }
      std::vector<int> result;
      result.reserve(array->size());
      for (const toml::node &node : *array) {
        const auto value = node.value<std::int64_t>();
        if (!value.has_value() || *value < 3 || *value > std::numeric_limits<int>::max()) {
          throw ConfigError(context + "." + key + " values must be integers >= 3.");
        }
        result.push_back(static_cast<int>(*value));
      }
      return result;
    }

    std::vector<double> ReadFiniteDoubleArray(const toml::table &table,
                                              const std::string &key,
                                              const std::string &context) {
      const auto *array = table[key].as_array();
      if (array == nullptr) {
        throw ConfigError("Expected numeric array '" + key + "' in " + context + ".");
      }
      std::vector<double> result;
      result.reserve(array->size());
      for (const toml::node &node : *array) {
        const auto value = node.value<double>();
        if (!value.has_value() || !std::isfinite(*value)) {
          throw ConfigError(context + "." + key + " must contain only finite values.");
        }
        result.push_back(*value);
      }
      return result;
    }

    bool IsSafeProfileScanId(const std::string &id) {
      if (id.empty()) {
        return false;
      }
      return std::all_of(id.begin(), id.end(), [](const unsigned char character) {
        return std::isalnum(character) || character == '_' || character == '-';
      });
    }

    ProfileSliceScope ParseProfileSliceScope(const std::string &value) {
      if (value == "listed") {
        return ProfileSliceScope::kListed;
      }
      if (value == "fit_selection") {
        return ProfileSliceScope::kFitSelection;
      }
      throw ConfigError("fit.profile_likelihood.slice_scope must be 'listed' or 'fit_selection'.");
    }

    ProfileRetryStrategy ParseProfileRetryStrategy(const std::string &value) {
      if (value == "reference_and_bidirectional_neighbors") {
        return ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors;
      }
      if (value == "reference_only") {
        return ProfileRetryStrategy::kReferenceOnly;
      }
      throw ConfigError("Unsupported fit.profile_likelihood.retry_strategy: " + value);
    }

    ProfileExecutionMode ParseProfileExecutionMode(const std::string &value) {
      if (value == "alongside_fit") return ProfileExecutionMode::kAlongsideFit;
      if (value == "profile_only") return ProfileExecutionMode::kProfileOnly;
      throw ConfigError("fit.profile_likelihood.execution_mode must be 'alongside_fit' or 'profile_only'.");
    }

    ProfileParallelBackend ParseProfileParallelBackend(const std::string &value) {
      if (value == "serial") return ProfileParallelBackend::kSerial;
      if (value == "process") return ProfileParallelBackend::kProcess;
      if (value == "thread") return ProfileParallelBackend::kThread;
      throw ConfigError("Unsupported fit.profile_likelihood.parallel_backend: " + value);
    }

    ProfileMinimizerBackend ParseProfileMinimizerBackend(const std::string &value) {
      if (value == "legacy_tminuit") return ProfileMinimizerBackend::kLegacyTMinuit;
      if (value == "minuit2") return ProfileMinimizerBackend::kMinuit2;
      throw ConfigError("Unsupported fit.profile_likelihood.minimizer_backend: " + value);
    }

    ProfileHesseStrategy ParseProfileHesseStrategy(const std::string &value) {
      if (value == "all_attempts") return ProfileHesseStrategy::kAllAttempts;
      if (value == "none") return ProfileHesseStrategy::kNone;
      throw ConfigError("Unsupported fit.profile_likelihood.hesse_strategy: " + value);
    }

    ProfileLikelihoodConfig ParseProfileLikelihoodConfig(const toml::table *table) {
      ProfileLikelihoodConfig config;
      if (table == nullptr) {
        return config;
      }
      const std::set<std::string> allowed_profile_fields = {"enabled", "slice_scope", "slice_ids", "retry_strategy",
                                                            "write_likelihood_slice", "contour_levels", "scans",
                                                            "execution_mode", "parallel_backend", "workers",
                                                            "minimizer_backend", "hesse_strategy", "checkpoint"};
      for (const auto &[raw_key, node] : *table) {
        (void)node;
        if (allowed_profile_fields.count(std::string(raw_key.str())) == 0U) {
          throw ConfigError("Unsupported field in fit.profile_likelihood: " + std::string(raw_key.str()));
        }
      }
      config.enabled = ReadOptionalBool(*table, "enabled", config.enabled);
      config.slice_scope = ParseProfileSliceScope(
          ReadOptionalString(*table, "slice_scope", "listed"));
      config.retry_strategy = ParseProfileRetryStrategy(ReadOptionalString(
          *table, "retry_strategy", "reference_and_bidirectional_neighbors"));
      config.write_likelihood_slice = ReadOptionalBool(*table, "write_likelihood_slice", true);
      config.execution_mode = ParseProfileExecutionMode(
          ReadOptionalString(*table, "execution_mode", "alongside_fit"));
      config.parallel_backend = ParseProfileParallelBackend(
          ReadOptionalString(*table, "parallel_backend", "serial"));
      config.minimizer_backend = ParseProfileMinimizerBackend(
          ReadOptionalString(*table, "minimizer_backend", "legacy_tminuit"));
      config.hesse_strategy = ParseProfileHesseStrategy(
          ReadOptionalString(*table, "hesse_strategy", "all_attempts"));
      if (table->contains("workers")) {
        const auto workers = (*table)["workers"].value<std::int64_t>();
        if (!workers.has_value() || *workers < 1 || *workers > std::numeric_limits<int>::max()) {
          throw ConfigError("fit.profile_likelihood.workers must be a positive integer.");
        }
        config.workers = static_cast<int>(*workers);
      }
      if (const auto *checkpoint = (*table)["checkpoint"].as_table(); checkpoint != nullptr) {
        const std::set<std::string> allowed_checkpoint_fields = {"enabled", "resume", "run_id", "directory"};
        for (const auto &[raw_key, node] : *checkpoint) {
          (void)node;
          if (allowed_checkpoint_fields.count(std::string(raw_key.str())) == 0U) {
            throw ConfigError("Unsupported field in fit.profile_likelihood.checkpoint: " + std::string(raw_key.str()));
          }
        }
        config.checkpoint.enabled = ReadOptionalBool(*checkpoint, "enabled", false);
        config.checkpoint.resume = ReadOptionalBool(*checkpoint, "resume", false);
        config.checkpoint.run_id = ReadOptionalString(*checkpoint, "run_id", "");
        config.checkpoint.directory = ReadOptionalString(*checkpoint, "directory", "");
      } else if (table->contains("checkpoint")) {
        throw ConfigError("fit.profile_likelihood.checkpoint must be a table.");
      }
      if (table->contains("slice_ids")) {
        config.slice_ids = ReadStringArray(*table, "slice_ids", "fit.profile_likelihood");
      }
      if (table->contains("contour_levels")) {
        config.contour_levels = ReadFiniteDoubleArray(*table, "contour_levels", "fit.profile_likelihood");
      }
      const auto *scans = (*table)["scans"].as_array();
      if (table->contains("scans") && scans == nullptr) {
        throw ConfigError("fit.profile_likelihood.scans must be an array of tables.");
      }
      if (scans != nullptr) {
        for (const toml::node &node : *scans) {
          const auto *scan_table = node.as_table();
          if (scan_table == nullptr) {
            throw ConfigError("Each [[fit.profile_likelihood.scans]] entry must be a table.");
          }
          const std::set<std::string> allowed_scan_fields = {"id", "parameters", "points", "min", "max", "refine",
                                                              "refinement_points"};
          for (const auto &[raw_key, child] : *scan_table) {
            (void)child;
            if (allowed_scan_fields.count(std::string(raw_key.str())) == 0U) {
              throw ConfigError("Unsupported field in fit.profile_likelihood.scans: " + std::string(raw_key.str()));
            }
          }
          ProfileScanConfig scan;
          scan.id = ReadRequiredString(*scan_table, "id", "fit.profile_likelihood.scans");
          scan.parameters = ReadStringArray(*scan_table, "parameters", "fit.profile_likelihood.scans");
          scan.points = ReadPointArray(*scan_table, "points", "fit.profile_likelihood.scans");
          const bool has_min = scan_table->contains("min");
          const bool has_max = scan_table->contains("max");
          if (has_min != has_max) {
            throw ConfigError("Profile scan min and max must be specified together.");
          }
          if (has_min) {
            scan.min = ReadFiniteDoubleArray(*scan_table, "min", "fit.profile_likelihood.scans");
            scan.max = ReadFiniteDoubleArray(*scan_table, "max", "fit.profile_likelihood.scans");
          }
          scan.refine = ReadOptionalBool(*scan_table, "refine", false);
          if (scan_table->contains("refinement_points")) {
            scan.refinement_points =
                ReadPointArray(*scan_table, "refinement_points", "fit.profile_likelihood.scans");
          }
          config.scans.push_back(std::move(scan));
        }
      }
      return config;
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

    bool IsAllowedParameterField(const std::string &field) {
      return field == "initial" || field == "min" || field == "max" || field == "fixed_value";
    }

    bool IsFixableLevyParameter(const std::string &parameter_name) {
      return parameter_name == "lambda" || parameter_name == "alpha";
    }

    bool HasAnyOverride(const LevyFitParameterOverride &parameter) {
      return parameter.initial.has_value() || parameter.min.has_value() || parameter.max.has_value()
             || parameter.fixed_value.has_value();
    }

    void ValidateFixedValueInsideEffectiveBounds(const std::string &context,
                                                 const LevyFitParameterOverride &parameter,
                                                 const double default_min,
                                                 const double default_max) {
      if (!parameter.fixed_value.has_value()) {
        return;
      }
      const double effective_min = parameter.min.value_or(default_min);
      const double effective_max = parameter.max.value_or(default_max);
      if (*parameter.fixed_value < effective_min || *parameter.fixed_value > effective_max) {
        throw ConfigError(context + ".fixed_value must be inside the effective min/max range.");
      }
    }

    LevyFitParameterOverride *FindParameterOverride(LevyFitParameterOverrides &parameters,
                                                    const std::string &parameter_name) {
      if (parameter_name == "norm") {
        return &parameters.norm;
      }
      if (parameter_name == "lambda") {
        return &parameters.lambda;
      }
      if (parameter_name == "rout2") {
        return &parameters.rout2;
      }
      if (parameter_name == "rside2") {
        return &parameters.rside2;
      }
      if (parameter_name == "rlong2") {
        return &parameters.rlong2;
      }
      if (parameter_name == "routside2") {
        return &parameters.routside2;
      }
      if (parameter_name == "routlong2") {
        return &parameters.routlong2;
      }
      if (parameter_name == "rsidelong2") {
        return &parameters.rsidelong2;
      }
      if (parameter_name == "alpha") {
        return &parameters.alpha;
      }
      if (parameter_name == "baseline_q2") {
        return &parameters.baseline_q2;
      }
      return nullptr;
    }

    LevyFitParameterOverride ParseLevyFitParameterOverride(const toml::table &table,
                                                           const std::string &parameter_name) {
      const std::string context = "fit.parameters." + parameter_name;
      for (const auto &[raw_key, node] : table) {
        (void)node;
        const std::string key(raw_key.str());
        if (!IsAllowedParameterField(key)) {
          throw ConfigError("Unsupported field '" + key + "' in " + context + ".");
        }
      }

      LevyFitParameterOverride parameter;
      parameter.initial = ReadOptionalFiniteDouble(table, "initial", context);
      parameter.min = ReadOptionalFiniteDouble(table, "min", context);
      parameter.max = ReadOptionalFiniteDouble(table, "max", context);
      parameter.fixed_value = ReadOptionalFiniteDouble(table, "fixed_value", context);

      if (parameter.min.has_value() != parameter.max.has_value()) {
        throw ConfigError(context + " must define both min and max when overriding parameter limits.");
      }
      if (parameter.min.has_value() && *parameter.min >= *parameter.max) {
        throw ConfigError(context + " must satisfy min < max.");
      }
      if (parameter.fixed_value.has_value() && !IsFixableLevyParameter(parameter_name)) {
        throw ConfigError(context + ".fixed_value is only supported for lambda and alpha.");
      }
      if (parameter.fixed_value.has_value() && parameter.min.has_value()
          && (*parameter.fixed_value < *parameter.min || *parameter.fixed_value > *parameter.max)) {
        throw ConfigError(context + ".fixed_value must be inside the configured min/max range.");
      }
      return parameter;
    }

    // Map supported Levy parameter tables into one validated contract shared by both fit objectives.
    LevyFitParameterOverrides ParseLevyFitParameterOverrides(const toml::table *parameters_table) {
      LevyFitParameterOverrides parameters;
      if (parameters_table == nullptr) {
        return parameters;
      }

      for (const auto &[raw_key, node] : *parameters_table) {
        const std::string parameter_name(raw_key.str());
        LevyFitParameterOverride *target = FindParameterOverride(parameters, parameter_name);
        if (target == nullptr) {
          throw ConfigError("Unsupported Levy fit parameter: fit.parameters." + parameter_name + ".");
        }
        const auto *parameter_table = node.as_table();
        if (parameter_table == nullptr) {
          throw ConfigError("Expected table for fit.parameters." + parameter_name + ".");
        }
        *target = ParseLevyFitParameterOverride(*parameter_table, parameter_name);
      }
      return parameters;
    }

    // Rebin tables are optional for legacy compatibility, but complete and strongly typed when present.
    AxisRebinConfig ParseAxisRebinConfig(const toml::table *rebin_table, const std::string &axis_name) {
      AxisRebinConfig config;
      if (rebin_table == nullptr) {
        return config;
      }

      const toml::node *axis_node = rebin_table->get(axis_name);
      if (axis_node == nullptr) {
        return config;
      }
      const auto *axis_table = axis_node->as_table();
      if (axis_table == nullptr) {
        throw ConfigError("build.rebin." + axis_name + " must be a table.");
      }

      config.configured = true;
      const auto enabled = (*axis_table)["enabled"].value<bool>();
      if (!enabled.has_value()) {
        throw ConfigError("build.rebin." + axis_name + ".enabled must be specified as a boolean.");
      }
      config.enabled = *enabled;

      if (axis_table->contains("factor")) {
        const auto factor = (*axis_table)["factor"].value<std::int64_t>();
        if (!factor.has_value() || *factor < std::numeric_limits<int>::min()
            || *factor > std::numeric_limits<int>::max()) {
          throw ConfigError("build.rebin." + axis_name + ".factor must be an integer.");
        }
        config.factor = static_cast<int>(*factor);
      }
      if (axis_table->contains("min")) {
        const auto min = (*axis_table)["min"].value<double>();
        if (!min.has_value()) {
          throw ConfigError("build.rebin." + axis_name + ".min must be numeric.");
        }
        config.min = *min;
      }
      if (axis_table->contains("max")) {
        const auto max = (*axis_table)["max"].value<double>();
        if (!max.has_value()) {
          throw ConfigError("build.rebin." + axis_name + ".max must be numeric.");
        }
        config.max = *max;
      }
      return config;
    }

    // Resolve the configuration-only mode; physical edge and factor checks wait for the ROOT axis.
    void ValidateAxisRebinConfig(const std::string &axis_name,
                                 AxisRebinConfig &rebin,
                                 const std::vector<RangeBin> &ranges,
                                 const bool allow_legacy_ranges) {
      const std::string context = "build.rebin." + axis_name;
      if (!rebin.configured) {
        if (!ranges.empty() && !allow_legacy_ranges) {
          throw ConfigError("[[bins." + axis_name + "]] requires an explicit " + context + " table.");
        }
        rebin.enabled = false;
        rebin.mode = ranges.empty() ? RebinMode::kNative : RebinMode::kLegacyRanges;
        return;
      }

      const bool has_factor = rebin.factor.has_value();
      const bool has_ranges = !ranges.empty();
      const bool has_min = rebin.min.has_value();
      const bool has_max = rebin.max.has_value();
      if (has_min != has_max) {
        throw ConfigError(context + ".min and .max must be specified together.");
      }
      if (has_min && (!std::isfinite(*rebin.min) || !std::isfinite(*rebin.max) || *rebin.max <= *rebin.min)) {
        throw ConfigError(context + ".min/.max must define a finite increasing range.");
      }

      if (!rebin.enabled) {
        if (has_factor || has_ranges || has_min) {
          throw ConfigError(context + " is disabled and cannot contain factor, min/max, or [[bins." + axis_name
                            + "]] ranges.");
        }
        rebin.mode = RebinMode::kNative;
        return;
      }

      if (has_factor == has_ranges) {
        throw ConfigError(context + " enabled=true requires exactly one of factor or [[bins." + axis_name + "]].");
      }
      if (has_factor) {
        if (*rebin.factor < 2) {
          throw ConfigError(context + ".factor must be >= 2.");
        }
        rebin.mode = RebinMode::kFactor;
        return;
      }
      if (has_min) {
        throw ConfigError(context + ".min/.max are only valid in factor mode.");
      }
      rebin.mode = RebinMode::kRanges;
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
    config.output.fit_report_directory =
        ReadOptionalString(output, "fit_report_directory", config.output.fit_report_directory);
    config.output.fit_report_root_name =
        ReadOptionalString(output, "fit_report_root_name", config.output.fit_report_root_name);
    config.output.profile_root_name =
        ReadOptionalString(output, "profile_root_name", config.output.profile_root_name);
    config.output.log_level = ParseLogLevel(ReadOptionalString(output, "log_level", ToString(config.output.log_level)));

    config.build.map_pair_phi_to_symmetric_range =
        ReadOptionalBool(build, "map_pair_phi_to_symmetric_range", config.build.map_pair_phi_to_symmetric_range);
    config.build.write_normalized_se_me_1d_projections = ReadOptionalBool(
        build, "write_normalized_se_me_1d_projections", config.build.write_normalized_se_me_1d_projections);
    config.build.reopen_output_file_per_slice =
        ReadOptionalBool(build, "reopen_output_file_per_slice", config.build.reopen_output_file_per_slice);
    config.build.split_mixed_event_by_phi =
        ReadOptionalBool(build, "split_mixed_event_by_phi", config.build.split_mixed_event_by_phi);
    config.build.split_same_event_by_qn =
        ReadOptionalBool(build, "split_same_event_by_qn", config.build.split_same_event_by_qn);
    config.build.split_mixed_event_by_qn =
        ReadOptionalBool(build, "split_mixed_event_by_qn", config.build.split_mixed_event_by_qn);
    config.build.progress = ReadOptionalProgressMode(build, "progress", config.build.progress);

    const toml::table *rebin = build["rebin"].as_table();
    if (build.contains("rebin") && rebin == nullptr) {
      throw ConfigError("build.rebin must be a table.");
    }
    config.build.phi_rebin = ParseAxisRebinConfig(rebin, "phi");
    config.build.mt_rebin = ParseAxisRebinConfig(rebin, "mt");

    config.fit.model = ParseFitModel(ReadOptionalString(fit, "model", ToString(config.fit.model)));
    const std::optional<bool> legacy_use_coulomb = ReadOptionalNullableBool(fit, "use_coulomb");
    const auto explicit_coulomb_mode = fit["coulomb_mode"].value<std::string>();
    if (explicit_coulomb_mode.has_value()) {
      config.fit.options.coulomb_mode = ParseCoulombMode(*explicit_coulomb_mode);
      if (legacy_use_coulomb.has_value()) {
        const CoulombMode legacy_mode = *legacy_use_coulomb ? CoulombMode::kGamow : CoulombMode::kNone;
        if (legacy_mode != config.fit.options.coulomb_mode) {
          throw ConfigError("Conflicting fit.use_coulomb and fit.coulomb_mode values.");
        }
      }
    } else if (legacy_use_coulomb.has_value()) {
      config.fit.options.coulomb_mode = *legacy_use_coulomb ? CoulombMode::kGamow : CoulombMode::kNone;
    }

    const auto explicit_finite_source_mode = fit["finite_source_mode"].value<std::string>();
    if (explicit_finite_source_mode.has_value()) {
      config.fit.options.finite_source_mode = ParseFiniteSourceMode(*explicit_finite_source_mode);
      if (config.fit.options.coulomb_mode != CoulombMode::kFiniteSource) {
        throw ConfigError("fit.finite_source_mode is only valid with fit.coulomb_mode = \"finite_source\".");
      }
    }
    config.fit.options.use_core_halo_lambda =
        ReadOptionalBool(fit, "use_core_halo_lambda", config.fit.options.use_core_halo_lambda);
    config.fit.options.use_q2_baseline = ReadOptionalBool(fit, "use_q2_baseline", config.fit.options.use_q2_baseline);
    config.fit.options.use_pml = ReadOptionalBool(fit, "use_pml", config.fit.options.use_pml);
    config.fit.options.fit_q_max = ReadOptionalDouble(fit, "fit_q_max", config.fit.options.fit_q_max);
    const toml::table *fit_parameters = fit["parameters"].as_table();
    if (fit.contains("parameters") && fit_parameters == nullptr) {
      throw ConfigError("fit.parameters must be a table.");
    }
    config.fit.options.parameters = ParseLevyFitParameterOverrides(fit_parameters);
    config.fit.map_pair_phi_to_symmetric_range = ReadOptionalNullableBool(fit, "map_pair_phi_to_symmetric_range");
    config.fit.reopen_output_file_per_slice =
        ReadOptionalBool(fit, "reopen_output_file_per_slice", config.fit.reopen_output_file_per_slice);
    config.fit.progress = ReadOptionalProgressMode(fit, "progress", config.fit.progress);
    const toml::table *profile_likelihood = fit["profile_likelihood"].as_table();
    if (fit.contains("profile_likelihood") && profile_likelihood == nullptr) {
      throw ConfigError("fit.profile_likelihood must be a table.");
    }
    config.fit.profile_likelihood = ParseProfileLikelihoodConfig(profile_likelihood);

    if (const auto *bins = root["bins"].as_table(); bins != nullptr) {
      config.centrality_bins = ParseRangeBinArray(GetOptionalArray(*bins, "centrality"), "bins.centrality");
      config.mt_bins = ParseRangeBinArray(GetOptionalArray(*bins, "mt"), "bins.mt");
      config.phi_bins = ParseRangeBinArray(GetOptionalArray(*bins, "phi"), "bins.phi");
      config.qn_bins = ParseRangeBinArray(GetOptionalArray(*bins, "qn"), "bins.qn");
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
    // Keep older configs valid while exposing an independently configurable report path.
    if (config.output.fit_report_directory.empty()) {
      config.output.fit_report_directory = config.output.output_directory;
    }
    if (config.output.fit_report_root_name.empty()) {
      throw ConfigError("output.fit_report_root_name is required.");
    }
    if (config.output.profile_root_name.empty()) {
      throw ConfigError("output.profile_root_name is required.");
    }
    if (config.fit.options.fit_q_max <= 0.0 || !std::isfinite(config.fit.options.fit_q_max)) {
      throw ConfigError("fit.fit_q_max must be finite and positive.");
    }
    if (!config.fit.options.use_core_halo_lambda && HasAnyOverride(config.fit.options.parameters.lambda)) {
      throw ConfigError("fit.parameters.lambda requires fit.use_core_halo_lambda = true.");
    }
    if (!config.fit.options.use_q2_baseline && HasAnyOverride(config.fit.options.parameters.baseline_q2)) {
      throw ConfigError("fit.parameters.baseline_q2 requires fit.use_q2_baseline = true.");
    }
    ValidateFixedValueInsideEffectiveBounds("fit.parameters.lambda", config.fit.options.parameters.lambda, 0.0, 1.0);
    ValidateFixedValueInsideEffectiveBounds("fit.parameters.alpha", config.fit.options.parameters.alpha, 0.5, 2.0);

    ValidateRangeCollection("centrality", config.centrality_bins, true);
    if (!config.mt_bins.empty()) {
      ValidateRangeCollection("mt", config.mt_bins, true);
    }
    if (!config.phi_bins.empty()) {
      ValidateRangeCollection("phi", config.phi_bins, true);
    }
    ValidateAxisRebinConfig("mt", config.build.mt_rebin, config.mt_bins, true);
    ValidateAxisRebinConfig("phi", config.build.phi_rebin, config.phi_bins, false);
    if (!config.build.mt_rebin.configured && config.mt_bins.empty()) {
      throw ConfigError("Missing [[bins.mt]] for legacy mode; configure build.rebin.mt for native/factor mode.");
    }
    if (config.build.split_mixed_event_by_qn && !config.build.split_same_event_by_qn) {
      throw ConfigError("build.split_mixed_event_by_qn requires build.split_same_event_by_qn = true.");
    }
    if (config.build.split_same_event_by_qn && config.qn_bins.empty()) {
      throw ConfigError("build.split_same_event_by_qn requires at least one [[bins.qn]] bin.");
    }
    for (std::size_t index = 0; index < config.qn_bins.size(); ++index) {
      if (config.qn_bins[index].label.empty()) {
        config.qn_bins[index].label = "qn" + std::to_string(index + 1U);
      }
    }
    if (!config.qn_bins.empty()) {
      ValidateRangeCollection("qn", config.qn_bins, false);
    }

    if (config.fit_centrality_bins.empty()) {
      config.fit_centrality_bins = config.centrality_bins;
    }
    ValidateRangeCollection("fit_selection.centrality", config.fit_centrality_bins, true);
    if (!config.fit_mt_bins.empty()) {
      ValidateRangeCollection("fit_selection.mt", config.fit_mt_bins, true);
    }

    for (const RangeBin &bin : config.fit_centrality_bins) {
      if (!ContainsExactRange(config.centrality_bins, bin)) {
        throw ConfigError("Each fit_selection.centrality bin must exactly match a build bin.");
      }
    }
    if (config.build.mt_rebin.mode == RebinMode::kRanges
        || config.build.mt_rebin.mode == RebinMode::kLegacyRanges) {
      for (const RangeBin &bin : config.fit_mt_bins) {
        if (!ContainsExactRange(config.mt_bins, bin)) {
          throw ConfigError("Each fit_selection.mt bin must exactly match a configured build bin.");
        }
      }
    }

    EnsureExtension(config.output.cf_root_name, ".root");
    EnsureExtension(config.output.fit_root_name, ".root");
    EnsureExtension(config.output.fit_summary_name, ".tsv");
    EnsureExtension(config.output.fit_report_root_name, ".root");
    EnsureExtension(config.output.profile_root_name, ".root");

    const ProfileLikelihoodConfig &profile = config.fit.profile_likelihood;
    if (!profile.enabled) {
      return;
    }
    if (!config.fit.options.use_pml) {
      throw ConfigError("fit.profile_likelihood.enabled requires fit.use_pml = true.");
    }
    if (profile.execution_mode == ProfileExecutionMode::kProfileOnly
        && profile.slice_scope != ProfileSliceScope::kListed) {
      throw ConfigError("fit.profile_likelihood.execution_mode = 'profile_only' requires slice_scope = 'listed'.");
    }
    if (profile.parallel_backend == ProfileParallelBackend::kProcess
        && profile.minimizer_backend != ProfileMinimizerBackend::kLegacyTMinuit) {
      throw ConfigError("fit.profile_likelihood.parallel_backend = 'process' requires minimizer_backend = 'legacy_tminuit'.");
    }
    if (profile.parallel_backend == ProfileParallelBackend::kProcess
        && profile.execution_mode != ProfileExecutionMode::kProfileOnly) {
      throw ConfigError("fit.profile_likelihood.parallel_backend = 'process' currently requires execution_mode = 'profile_only'.");
    }
    if (profile.parallel_backend == ProfileParallelBackend::kThread
        && profile.minimizer_backend != ProfileMinimizerBackend::kMinuit2) {
      throw ConfigError("fit.profile_likelihood.parallel_backend = 'thread' requires minimizer_backend = 'minuit2'.");
    }
    if (profile.workers < 1) {
      throw ConfigError("fit.profile_likelihood.workers must be positive.");
    }
    if (profile.checkpoint.resume && !profile.checkpoint.enabled) {
      throw ConfigError("fit.profile_likelihood.checkpoint.resume requires checkpoint.enabled = true.");
    }
    if (profile.checkpoint.enabled
        && (profile.checkpoint.run_id.empty() || profile.checkpoint.directory.empty())) {
      throw ConfigError("Enabled fit.profile_likelihood.checkpoint requires non-empty run_id and directory.");
    }
    if (profile.scans.empty()) {
      throw ConfigError("fit.profile_likelihood.enabled requires at least one [[fit.profile_likelihood.scans]] entry.");
    }
    if (profile.slice_scope == ProfileSliceScope::kListed) {
      if (profile.slice_ids.empty()) {
        throw ConfigError("fit.profile_likelihood.slice_ids is required when slice_scope = 'listed'.");
      }
      std::set<std::string> unique_slice_ids(profile.slice_ids.begin(), profile.slice_ids.end());
      if (unique_slice_ids.size() != profile.slice_ids.size()) {
        throw ConfigError("fit.profile_likelihood.slice_ids must be unique.");
      }
    } else if (!profile.slice_ids.empty()) {
      throw ConfigError("fit.profile_likelihood.slice_ids is forbidden when slice_scope = 'fit_selection'.");
    }
    if (profile.contour_levels.empty()
        || std::any_of(profile.contour_levels.begin(), profile.contour_levels.end(), [](const double value) {
             return !std::isfinite(value) || value <= 0.0;
           })) {
      throw ConfigError("fit.profile_likelihood.contour_levels must contain finite positive values.");
    }
    std::set<std::string> scan_ids;
    for (const ProfileScanConfig &scan : profile.scans) {
      if (!IsSafeProfileScanId(scan.id) || !scan_ids.insert(scan.id).second) {
        throw ConfigError("Profile scan id must be unique and contain only letters, digits, '_' or '-'.");
      }
      if (scan.parameters.empty() || scan.parameters.size() > 2U || scan.points.size() != scan.parameters.size()) {
        throw ConfigError("Profile scans require matching one- or two-dimensional parameters and points arrays.");
      }
      if (scan.parameters.size() == 2U && scan.parameters[0] == scan.parameters[1]) {
        throw ConfigError("A two-dimensional profile scan cannot repeat a parameter target.");
      }
      if (!scan.min.empty()
          && (scan.min.size() != scan.parameters.size() || scan.max.size() != scan.parameters.size())) {
        throw ConfigError("Profile scan min/max arrays must match parameters dimensionality.");
      }
      for (std::size_t axis = 0; axis < scan.min.size(); ++axis) {
        if (!(scan.max[axis] > scan.min[axis])) {
          throw ConfigError("Profile scan min/max must be finite increasing ranges.");
        }
      }
      if (scan.refine && scan.refinement_points.empty()) {
        throw ConfigError("refine = true requires refinement_points.");
      }
      if (scan.refine && scan.refinement_points.size() != scan.parameters.size()) {
        throw ConfigError("Profile scan refinement_points must match parameters dimensionality.");
      }
      if (!scan.refine && !scan.refinement_points.empty()) {
        throw ConfigError("refinement_points requires refine = true.");
      }
    }
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

  std::string ToString(const CoulombMode mode) {
    switch (mode) {
      case CoulombMode::kNone:
        return "none";
      case CoulombMode::kGamow:
        return "gamow";
      case CoulombMode::kFiniteSource:
        return "finite_source";
    }
    return "none";
  }

  CoulombMode ParseCoulombMode(const std::string &token) {
    const std::string lowered = ToLower(token);
    if (lowered == "none") {
      return CoulombMode::kNone;
    }
    if (lowered == "gamow") {
      return CoulombMode::kGamow;
    }
    if (lowered == "finite_source" || lowered == "finite-source") {
      return CoulombMode::kFiniteSource;
    }
    throw ConfigError("Unsupported Coulomb mode: " + token);
  }

  std::string ToString(const FiniteSourceMode mode) {
    switch (mode) {
      case FiniteSourceMode::kFixed1D:
        return "fixed_1d";
      case FiniteSourceMode::kIterative1D:
        return "iterative_1d";
    }
    return "fixed_1d";
  }

  FiniteSourceMode ParseFiniteSourceMode(const std::string &token) {
    const std::string lowered = ToLower(token);
    if (lowered == "fixed_1d" || lowered == "fixed-1d") {
      return FiniteSourceMode::kFixed1D;
    }
    if (lowered == "iterative_1d" || lowered == "iterative-1d") {
      return FiniteSourceMode::kIterative1D;
    }
    throw ConfigError("Unsupported finite-source mode: " + token);
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

  std::string ToString(const RebinMode mode) {
    switch (mode) {
      case RebinMode::kNative:
        return "native";
      case RebinMode::kFactor:
        return "factor";
      case RebinMode::kRanges:
        return "ranges";
      case RebinMode::kLegacyRanges:
        return "legacy";
    }
    return "native";
  }

}  // namespace exp_femto_3d
