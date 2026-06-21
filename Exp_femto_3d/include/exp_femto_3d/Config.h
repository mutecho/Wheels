#pragma once

#include <string>

#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  [[nodiscard]] ApplicationConfig LoadApplicationConfig(const std::string &path);
  void ValidateApplicationConfig(ApplicationConfig &config);

  [[nodiscard]] std::string ToString(LogLevel level);
  [[nodiscard]] LogLevel ParseLogLevel(const std::string &token);
  [[nodiscard]] std::string ToString(FitModel model);
  [[nodiscard]] FitModel ParseFitModel(const std::string &token);
  [[nodiscard]] std::string ToString(CoulombMode mode);
  [[nodiscard]] CoulombMode ParseCoulombMode(const std::string &token);
  [[nodiscard]] std::string ToString(FiniteSourceMode mode);
  [[nodiscard]] FiniteSourceMode ParseFiniteSourceMode(const std::string &token);
  [[nodiscard]] std::string ToString(ProgressMode mode);
  [[nodiscard]] ProgressMode ParseProgressMode(const std::string &token);

}  // namespace exp_femto_3d
