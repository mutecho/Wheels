#pragma once

#include "femto3d/AnalysisConfig.h"

#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>

namespace femto3d {

class ConfigError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

struct CliOptions {
  std::string config_path;
  std::optional<std::string> input_root_override;
  std::optional<std::string> output_root_override;
  std::optional<InputSchema> input_schema_override;
};

[[nodiscard]] CliOptions ParseCliArgs(int argc, char** argv);

void PrintUsage(std::ostream& stream);

[[nodiscard]] ApplicationConfig LoadApplicationConfig(const std::string& path);

void ValidateApplicationConfig(ApplicationConfig& config);

void ApplyCliOverrides(const CliOptions& cli_options, ApplicationConfig& config);

[[nodiscard]] std::string ToString(InputSchema schema);

[[nodiscard]] InputSchema ParseInputSchema(const std::string& token);

[[nodiscard]] std::string ToString(EventPlaneWeightMode mode);

[[nodiscard]] EventPlaneWeightMode ParseEventPlaneWeightMode(
    const std::string& token);

}  // namespace femto3d
