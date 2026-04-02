#pragma once

#include "femto3d/AnalysisConfig.h"
#include "femto3d/AnalysisTypes.h"

#include <array>
#include <limits>
#include <string>

namespace femto3d {

inline void ApplyAlphaErrorValidityMask(const bool error_valid,
                                        double* error) {
  if (!error_valid && error != nullptr) {
    *error = std::numeric_limits<double>::quiet_NaN();
  }
}

inline void ApplyHbtErrorValidityMask(const std::array<bool, 6>& valid_flags,
                                      double errors[6]) {
  for (std::size_t index = 0; index < valid_flags.size(); ++index) {
    if (!valid_flags[index]) {
      errors[index] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

inline void ApplyDirectionalErrorValidityMask(
    const std::array<bool, 6>& valid_flags,
    double errors[6]) {
  for (std::size_t index = 0; index < valid_flags.size(); ++index) {
    if (!valid_flags[index]) {
      errors[index] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

[[nodiscard]] AnalysisStatistics RunAnalysis(const AnalysisConfig& config,
                                             const std::string& input_root_path,
                                             const std::string& output_root_path);

}  // namespace femto3d
