#pragma once

#include "femto3d/AnalysisConfig.h"
#include "femto3d/AnalysisTypes.h"

#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace femto3d {

struct R2SummaryPointDecision {
  bool write_point = false;
  double error = 0.0;
  bool skipped_invalid_hbt_error = false;
};

// R2 summary plots keep the relaxed "central-value-only" fallback separate
// from fit_summary masking so callers can decide whether missing errors should
// hide a point or only remove its error bar.
inline R2SummaryPointDecision DecideHbtR2SummaryPoint(
    const ProjectionFitResult& hbt_result,
    const bool hbt_error_valid,
    const ProjectionFitConfig& fit_config) {
  R2SummaryPointDecision decision;
  if (!hbt_result.success || !std::isfinite(hbt_result.r2)) {
    return decision;
  }

  if (hbt_error_valid) {
    decision.write_point = true;
    decision.error = hbt_result.r2_error;
    return decision;
  }

  if (fit_config.accept_hbt_central_value_only_for_summary) {
    decision.write_point = true;
    decision.error = 0.0;
    return decision;
  }

  decision.skipped_invalid_hbt_error = true;
  return decision;
}

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

[[nodiscard]] AnalysisStatistics RunAnalysis(const ApplicationConfig& config);

}  // namespace femto3d
