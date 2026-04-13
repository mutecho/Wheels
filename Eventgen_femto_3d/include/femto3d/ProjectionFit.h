#pragma once

#include "femto3d/AnalysisTypes.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"

#include <array>
#include <memory>
#include <string>

namespace femto3d {

struct ProjectionFitConfig {
  bool use_adaptive_integration = false;
  bool accept_forced_posdef_covariance_as_valid = false;
  bool fail_alpha_result_when_error_invalid = false;
  bool fail_hbt_results_when_error_invalid = false;
  bool fail_directional_results_when_error_invalid = false;
  bool accept_hbt_central_value_only_for_summary = true;
};

struct SliceFitProducts {
  std::array<ProjectionFitResult, 6> directional_results;
  std::array<ProjectionFitResult, 6> hbt_radii_results;
  std::array<bool, 6> directional_error_valid{};
  std::array<bool, 6> hbt_error_valid{};
  std::array<std::array<double, 6>, 6> hbt_covariance{};
  std::array<std::array<double, 6>, 6> directional_covariance{};
  std::array<std::unique_ptr<TF1>, 6> fit_functions;
  std::unique_ptr<TCanvas> fit_canvas;
  double alpha = 0.0;
  double alpha_error = 0.0;
  bool alpha_success = false;
  bool alpha_error_valid = false;
  int covariance_status = 0;
  bool covariance_valid = false;
};

[[nodiscard]] SliceFitProducts FitSliceHistograms(
    const std::string& slice_name,
    std::array<std::unique_ptr<TH1D>, 6>& projection_histograms,
    const std::array<ProjectionDefinition, 6>& directions,
    ProjectionFitConfig fit_config = {});

}  // namespace femto3d
