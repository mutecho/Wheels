#include <cmath>
#include <iostream>
#include <limits>

#include "../src/ProjectionFit.cpp"

int main() {
  using femto3d::BuildAlphaResult;
  using femto3d::BuildHbtResult;

  femto3d::ProjectionFitConfig default_config;

  const auto alpha_valid = BuildAlphaResult(1.35, true, 0.08, default_config);
  if (!alpha_valid.success || !alpha_valid.error_valid ||
      std::abs(alpha_valid.value - 1.35) > 1.0e-12 ||
      std::abs(alpha_valid.error - 0.08) > 1.0e-12) {
    std::cerr << "Expected valid alpha result to keep finite value and error.\n";
    return 1;
  }

  const auto alpha_loose_invalid =
      BuildAlphaResult(1.42, false, 0.05, default_config);
  if (!alpha_loose_invalid.success || alpha_loose_invalid.error_valid ||
      std::abs(alpha_loose_invalid.value - 1.42) > 1.0e-12 ||
      !std::isnan(alpha_loose_invalid.error)) {
    std::cerr << "Expected loose alpha mode to keep center value and mask error.\n";
    return 2;
  }

  femto3d::ProjectionFitConfig strict_alpha_config;
  strict_alpha_config.fail_alpha_result_when_error_invalid = true;
  const auto alpha_strict_invalid =
      BuildAlphaResult(1.57, false, 0.04, strict_alpha_config);
  if (alpha_strict_invalid.success || alpha_strict_invalid.error_valid ||
      !std::isnan(alpha_strict_invalid.value) ||
      !std::isnan(alpha_strict_invalid.error)) {
    std::cerr << "Expected strict alpha mode to fail the result and mask value/error.\n";
    return 3;
  }

  const femto3d::ProjectionFitResult hbt_valid =
      BuildHbtResult("Rout2", 4.25, true, 0.30, default_config);
  if (!hbt_valid.success || std::isnan(hbt_valid.r2_error) ||
      std::abs(hbt_valid.r2 - 4.25) > 1.0e-12 ||
      std::abs(hbt_valid.r2_error - 0.30) > 1.0e-12) {
    std::cerr << "Expected valid HBT result to keep finite value and error.\n";
    return 4;
  }

  const femto3d::ProjectionFitResult hbt_loose_invalid =
      BuildHbtResult("Rside2", 5.50, false, 0.20, default_config);
  if (!hbt_loose_invalid.success || !std::isnan(hbt_loose_invalid.r2_error) ||
      std::abs(hbt_loose_invalid.r2 - 5.50) > 1.0e-12) {
    std::cerr << "Expected loose HBT mode to keep center value and mask error.\n";
    return 5;
  }

  femto3d::ProjectionFitConfig strict_hbt_config;
  strict_hbt_config.fail_hbt_results_when_error_invalid = true;
  const femto3d::ProjectionFitResult hbt_strict_invalid =
      BuildHbtResult("Rlong2", 6.75, false, 0.25, strict_hbt_config);
  if (hbt_strict_invalid.success || !std::isnan(hbt_strict_invalid.r2) ||
      !std::isnan(hbt_strict_invalid.r2_error)) {
    std::cerr << "Expected strict HBT mode to fail the result and mask value/error.\n";
    return 6;
  }

  return 0;
}
