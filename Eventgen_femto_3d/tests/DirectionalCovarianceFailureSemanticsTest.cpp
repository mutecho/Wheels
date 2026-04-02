#include <cmath>
#include <iostream>
#include <limits>

#include "../src/ProjectionFit.cpp"

int main() {
  using femto3d::BuildDirectionalResult;
  using femto3d::BuildSuccessfulResult;
  using femto3d::BuildSuccessfulCentralValueOnlyResult;

  const femto3d::ProjectionFitResult valid_result =
      BuildSuccessfulResult("out", 1.25, 0.15);
  if (!valid_result.success || std::isnan(valid_result.r2_error) ||
      std::abs(valid_result.r2_error - 0.15) > 1.0e-12) {
    std::cerr << "Expected valid directional result to preserve finite error.\n";
    return 1;
  }

  const femto3d::ProjectionFitResult central_only_result =
      BuildSuccessfulCentralValueOnlyResult("out_side", 2.50);
  if (!central_only_result.success || !std::isnan(central_only_result.r2_error) ||
      std::abs(central_only_result.r2 - 2.50) > 1.0e-12) {
    std::cerr << "Expected covariance-failure fallback to keep center value and "
                 "mark error invalid via NaN.\n";
    return 2;
  }

  femto3d::ProjectionFitConfig strict_failure_config;
  strict_failure_config.fail_directional_results_when_error_invalid = true;
  const femto3d::ProjectionFitResult strict_failure_result =
      BuildDirectionalResult("out_long",
                             3.75,
                             false,
                             std::numeric_limits<double>::quiet_NaN(),
                             strict_failure_config);
  if (strict_failure_result.success || !std::isnan(strict_failure_result.r2_error)) {
    std::cerr << "Expected strict covariance-failure mode to mark directional "
                 "result as failed with invalid error.\n";
    return 3;
  }

  return 0;
}
