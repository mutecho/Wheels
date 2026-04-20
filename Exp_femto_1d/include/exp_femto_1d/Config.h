#pragma once

#include <string>

#include "exp_femto_1d/Types.h"

namespace exp_femto_1d {

  [[nodiscard]] ApplicationConfig LoadApplicationConfig(const std::string &path);
  void ValidateApplicationConfig(ApplicationConfig &config);

  [[nodiscard]] std::string ToString(LogLevel level);
  [[nodiscard]] LogLevel ParseLogLevel(const std::string &token);
  [[nodiscard]] std::string ToString(ProgressMode mode);
  [[nodiscard]] ProgressMode ParseProgressMode(const std::string &token);
  [[nodiscard]] std::string ToString(RegionKind kind);

}  // namespace exp_femto_1d
