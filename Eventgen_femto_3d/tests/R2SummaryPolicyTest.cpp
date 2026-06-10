#include "femto3d/ProgressReporter.h"
#include "femto3d/Workflow.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

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

bool Contains(const std::string& text, const std::string& expected) {
  return text.find(expected) != std::string::npos;
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

  const femto3d::ProgressRenderSnapshot progress_snapshot{
      10U, 5U, 1, std::chrono::seconds(20)};
  const std::string progress_line =
      femto3d::FormatProgressLine(progress_snapshot);
  if (!Contains(progress_line, "50%") ||
      !Contains(progress_line, "ETA 00:00:20")) {
    std::cerr << "Expected progress line to include integer percent and ETA.\n";
    return 6;
  }

  const femto3d::RangeBin mt_bin{0.2, 0.4, "mt_test"};
  std::vector<femto3d::R2SummaryPoint> rside_points;
  for (const double phi :
       {-3.0 * femto3d::kPi / 8.0,
        -1.0 * femto3d::kPi / 8.0,
        1.0 * femto3d::kPi / 8.0,
        3.0 * femto3d::kPi / 8.0}) {
    femto3d::R2SummaryPoint point;
    point.phi_center = phi;
    point.phi_error = femto3d::kPi / 8.0;
    point.value = 10.0 + 2.0 * std::cos(2.0 * phi);
    point.error = 0.5;
    point.valid = true;
    rside_points.push_back(point);
  }

  const std::optional<femto3d::EpsfSummaryPoint> epsf_point =
      femto3d::ComputeEpsfFromRsideSummaryPoints(rside_points, mt_bin);
  if (!epsf_point.has_value() || !epsf_point->valid ||
      std::abs(epsf_point->value - 0.2) > 1.0e-12 ||
      std::abs(epsf_point->mt_center - 0.3) > 1.0e-12 ||
      epsf_point->error <= 0.0) {
    std::cerr << "Expected epsf extraction to recover side-radius second "
                 "harmonic over femto mT.\n";
    return 7;
  }

  rside_points.resize(1U);
  if (femto3d::ComputeEpsfFromRsideSummaryPoints(rside_points, mt_bin)
          .has_value()) {
    std::cerr << "Expected epsf extraction to require at least two usable "
                 "phi points.\n";
    return 8;
  }

  return 0;
}
