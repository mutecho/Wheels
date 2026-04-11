#pragma once

#include "exp_femto_3d/Types.h"

#include <string>

namespace exp_femto_3d {

[[nodiscard]] ApplicationConfig LoadApplicationConfig(const std::string& path);
void ValidateApplicationConfig(ApplicationConfig& config);

[[nodiscard]] std::string ToString(LogLevel level);
[[nodiscard]] LogLevel ParseLogLevel(const std::string& token);
[[nodiscard]] std::string ToString(FitModel model);
[[nodiscard]] FitModel ParseFitModel(const std::string& token);

}  // namespace exp_femto_3d
