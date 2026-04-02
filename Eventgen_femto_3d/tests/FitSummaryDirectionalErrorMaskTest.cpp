#include <cmath>
#include <iostream>

#include "../src/ProjectionFit.cpp"
#include "femto3d/Workflow.h"

namespace femto3d {

int RunFitSummaryDirectionalErrorMaskTest() {
  std::array<ProjectionFitResult, 6> results{};
  std::array<bool, 6> directional_error_valid{};
  std::array<bool, 6> hbt_error_valid{};

  results[0] = BuildSuccessfulResult("out", 1.0, 0.2);
  directional_error_valid[0] = true;

  results[1] = BuildSuccessfulCentralValueOnlyResult("side", 2.0);
  directional_error_valid[1] = false;

  // Simulate a future regression where an invalid directional point
  // accidentally carries a default zero error into fit_summary.
  results[2] = BuildSuccessfulResult("long", 3.0, 0.0);
  directional_error_valid[2] = false;

  double values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double errors[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int success_flags[6] = {0, 0, 0, 0, 0, 0};
  int error_valid_flags[6] = {0, 0, 0, 0, 0, 0};

  for (std::size_t index = 0; index < results.size(); ++index) {
    values[index] = results[index].r2;
    errors[index] = results[index].r2_error;
    success_flags[index] = results[index].success ? 1 : 0;
    error_valid_flags[index] = directional_error_valid[index] ? 1 : 0;
  }
  ApplyDirectionalErrorValidityMask(directional_error_valid, errors);

  if (success_flags[0] != 1 || error_valid_flags[0] != 1 ||
      std::abs(values[0] - 1.0) > 1.0e-12 ||
      std::abs(errors[0] - 0.2) > 1.0e-12) {
    std::cerr << "Expected valid directional point to remain unchanged in fit_summary.\n";
    return 1;
  }

  if (success_flags[1] != 1 || error_valid_flags[1] != 0 ||
      std::abs(values[1] - 2.0) > 1.0e-12 || !std::isnan(errors[1])) {
    std::cerr << "Expected center-only directional point to keep value and mask error to NaN.\n";
    return 2;
  }

  if (success_flags[2] != 1 || error_valid_flags[2] != 0 ||
      std::abs(values[2] - 3.0) > 1.0e-12 || !std::isnan(errors[2])) {
    std::cerr << "Expected fit_summary mask to prevent invalid directional zero error from leaking through.\n";
    return 3;
  }

  double alpha_error = 0.0;
  ApplyAlphaErrorValidityMask(false, &alpha_error);
  if (!std::isnan(alpha_error)) {
    std::cerr << "Expected alpha fit_summary mask to replace invalid error with NaN.\n";
    return 4;
  }

  double alpha_error_valid_value = 0.11;
  ApplyAlphaErrorValidityMask(true, &alpha_error_valid_value);
  if (std::abs(alpha_error_valid_value - 0.11) > 1.0e-12) {
    std::cerr << "Expected alpha fit_summary mask to preserve valid error.\n";
    return 5;
  }

  double hbt_errors[6] = {0.4, 0.0, 0.6, 0.0, 0.0, 0.0};
  hbt_error_valid[0] = true;
  hbt_error_valid[1] = false;
  hbt_error_valid[2] = true;
  ApplyHbtErrorValidityMask(hbt_error_valid, hbt_errors);
  if (std::abs(hbt_errors[0] - 0.4) > 1.0e-12 ||
      !std::isnan(hbt_errors[1]) ||
      std::abs(hbt_errors[2] - 0.6) > 1.0e-12) {
    std::cerr << "Expected HBT fit_summary mask to preserve valid errors and mask invalid ones.\n";
    return 6;
  }

  return 0;
}

}  // namespace femto3d

int main() {
  return femto3d::RunFitSummaryDirectionalErrorMaskTest();
}
