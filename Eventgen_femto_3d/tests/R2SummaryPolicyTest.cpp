#include "femto3d/Workflow.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace {

using femto3d::DecideHbtR2SummaryPoint;
using femto3d::ProjectionFitConfig;
using femto3d::ProjectionFitResult;

ProjectionFitResult MakeResult(const double r2,
                               const double r2_error,
                               const bool success) {
  ProjectionFitResult result;
  result.projection_name = "Rout2";
  result.r2 = r2;
  result.r2_error = r2_error;
  result.success = success;
  return result;
}

}  // namespace

int main() {
  ProjectionFitConfig accept_config;
  const ProjectionFitResult central_only_result =
      MakeResult(2.5, std::numeric_limits<double>::quiet_NaN(), true);
  const femto3d::R2SummaryPointDecision accepted_decision =
      DecideHbtR2SummaryPoint(central_only_result, false, accept_config);
  if (!accepted_decision.write_point || accepted_decision.error != 0.0 ||
      accepted_decision.skipped_invalid_hbt_error) {
    std::cerr << "Expected enabled summary fallback to keep a central-value-only "
                 "point with zero error bar.\n";
    return 1;
  }

  ProjectionFitConfig strict_config;
  strict_config.accept_hbt_central_value_only_for_summary = false;
  const femto3d::R2SummaryPointDecision rejected_decision =
      DecideHbtR2SummaryPoint(central_only_result, false, strict_config);
  if (rejected_decision.write_point || rejected_decision.error != 0.0 ||
      !rejected_decision.skipped_invalid_hbt_error) {
    std::cerr << "Expected strict summary policy to skip central-value-only "
                 "points and mark the skip reason.\n";
    return 2;
  }

  const ProjectionFitResult full_result = MakeResult(1.25, 0.15, true);
  const femto3d::R2SummaryPointDecision valid_error_decision =
      DecideHbtR2SummaryPoint(full_result, true, strict_config);
  if (!valid_error_decision.write_point ||
      std::abs(valid_error_decision.error - 0.15) > 1.0e-12 ||
      valid_error_decision.skipped_invalid_hbt_error) {
    std::cerr << "Expected valid HBT errors to remain visible regardless of "
                 "summary fallback policy.\n";
    return 3;
  }

  const ProjectionFitResult failed_result = MakeResult(3.0, 0.2, false);
  const femto3d::R2SummaryPointDecision failed_decision =
      DecideHbtR2SummaryPoint(failed_result, true, accept_config);
  if (failed_decision.write_point || failed_decision.skipped_invalid_hbt_error) {
    std::cerr << "Expected failed HBT results to stay out of summary graphs.\n";
    return 4;
  }

  const ProjectionFitResult nan_result =
      MakeResult(std::numeric_limits<double>::quiet_NaN(), 0.2, true);
  const femto3d::R2SummaryPointDecision nan_decision =
      DecideHbtR2SummaryPoint(nan_result, true, accept_config);
  if (nan_decision.write_point || nan_decision.skipped_invalid_hbt_error) {
    std::cerr << "Expected non-finite HBT central values to stay out of "
                 "summary graphs.\n";
    return 5;
  }

  return 0;
}
