#include "femto3d/ProjectionFit.h"

#include "femto3d/AnalysisTypes.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace femto3d {

namespace {

constexpr double kMinimumPositiveR2 = 1.0e-6;
constexpr double kMinimumAlpha = 0.20;
constexpr double kMaximumAlpha = 2.00;
constexpr double kLargePenalty = 1.0e12;
constexpr double kLevyIntegrationUpperBound = 60.0;
constexpr int kLevyIntegrationSteps = 1200;
constexpr std::size_t kLevyIntegrationPointCount =
    static_cast<std::size_t>(kLevyIntegrationSteps) + 1U;
constexpr double kAdaptiveAbsoluteTolerance = 1.0e-9;
constexpr double kAdaptiveRelativeTolerance = 1.0e-7;
constexpr int kAdaptiveMaxDepth = 20;
constexpr double kComparisonTolerance = 1.0e-12;
constexpr std::size_t kDirectionalCount = 6U;
constexpr std::size_t kFitParameterCount = 7U;
constexpr std::size_t kHbtParameterOffset = 1U;

using CovarianceMatrix6 =
    std::array<std::array<double, kDirectionalCount>, kDirectionalCount>;

struct FixedIntegrationWorkspace {
  std::array<double, kLevyIntegrationPointCount> u_grid{};
  std::array<double, kLevyIntegrationPointCount> simpson_weights{};
  double simpson_scale = 0.0;
  double inverse_pi = 1.0 / kPi;
};

struct RawProjectionSample {
  double abs_x = 0.0;
  double y = 0.0;
  double sigma = 0.0;
};

struct ProjectionHistogramCache {
  std::vector<std::size_t> abs_x_indices;
  std::vector<double> y_values;
  std::vector<double> sigma_values;
};

struct ProjectionFitCache {
  std::array<ProjectionHistogramCache, 6> histogram_caches;
  std::vector<double> unique_abs_x_grid;
  int used_points = 0;
};

struct ProjectedCurveSet {
  std::array<std::vector<double>, 6> curves;
  std::array<std::size_t, 6> direction_curve_indices{};
  std::array<double, 6> unique_scales{};
  std::size_t unique_scale_count = 0;
};

struct CovariancePropagationResult {
  CovarianceMatrix6 hbt_covariance{};
  CovarianceMatrix6 directional_covariance{};
  std::array<double, kDirectionalCount> directional_errors{};
  int covariance_status = 0;
  bool covariance_valid = false;
  bool directional_valid = false;
};

constexpr CovarianceMatrix6 kDirectionalTransformMatrix = {{
    {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}},
    {{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}},
    {{0.5, 0.5, 0.0, 1.0, 0.0, 0.0}},
    {{0.5, 0.0, 0.5, 0.0, 1.0, 0.0}},
    {{0.0, 0.5, 0.5, 0.0, 0.0, 1.0}},
}};

bool NearlyEqual(const double left, const double right) {
  return std::abs(left - right) <=
         kComparisonTolerance *
             std::max({1.0, std::abs(left), std::abs(right)});
}

const FixedIntegrationWorkspace& GetFixedIntegrationWorkspace() {
  static const FixedIntegrationWorkspace workspace = [] {
    FixedIntegrationWorkspace result;
    const double step =
        kLevyIntegrationUpperBound / static_cast<double>(kLevyIntegrationSteps);
    result.simpson_scale = step / 3.0;

    for (std::size_t index = 0; index < kLevyIntegrationPointCount; ++index) {
      result.u_grid[index] = static_cast<double>(index) * step;
      result.simpson_weights[index] =
          (index == 0U || index + 1U == kLevyIntegrationPointCount)
              ? 1.0
              : ((index % 2U == 0U) ? 2.0 : 4.0);
    }

    return result;
  }();

  return workspace;
}

std::vector<double> BuildUniqueAbsXGrid(std::vector<double> abs_x_values) {
  if (abs_x_values.empty()) {
    return {};
  }

  std::sort(abs_x_values.begin(), abs_x_values.end());

  std::vector<double> unique_abs_x_grid;
  unique_abs_x_grid.reserve(abs_x_values.size());
  for (const double abs_x : abs_x_values) {
    if (unique_abs_x_grid.empty() ||
        !NearlyEqual(unique_abs_x_grid.back(), abs_x)) {
      unique_abs_x_grid.push_back(abs_x);
    }
  }

  return unique_abs_x_grid;
}

std::size_t FindGridIndex(const std::vector<double>& unique_abs_x_grid,
                          const double abs_x) {
  const auto iterator =
      std::lower_bound(unique_abs_x_grid.begin(), unique_abs_x_grid.end(), abs_x);
  if (iterator != unique_abs_x_grid.end() && NearlyEqual(*iterator, abs_x)) {
    return static_cast<std::size_t>(
        std::distance(unique_abs_x_grid.begin(), iterator));
  }

  if (iterator != unique_abs_x_grid.begin()) {
    const auto previous = iterator - 1;
    if (NearlyEqual(*previous, abs_x)) {
      return static_cast<std::size_t>(
          std::distance(unique_abs_x_grid.begin(), previous));
    }
  }

  throw std::logic_error("Failed to locate cached |x| grid point.");
}

ProjectionFitCache BuildProjectionFitCache(
    const std::array<std::unique_ptr<TH1D>, 6>& normalized_histograms) {
  std::array<std::vector<RawProjectionSample>, 6> raw_samples;
  std::vector<double> abs_x_values;
  int used_points = 0;

  for (std::size_t histogram_index = 0; histogram_index < normalized_histograms.size();
       ++histogram_index) {
    const TH1D* histogram = normalized_histograms[histogram_index].get();
    if (histogram == nullptr) {
      continue;
    }

    auto& samples = raw_samples[histogram_index];
    samples.reserve(static_cast<std::size_t>(histogram->GetNbinsX()));

    for (int bin_index = 1; bin_index <= histogram->GetNbinsX(); ++bin_index) {
      const double sigma = histogram->GetBinError(bin_index);
      if (sigma <= 0.0) {
        continue;
      }

      const double abs_x = std::abs(histogram->GetBinCenter(bin_index));
      samples.push_back(
          {abs_x, histogram->GetBinContent(bin_index), sigma});
      abs_x_values.push_back(abs_x);
      ++used_points;
    }
  }

  ProjectionFitCache cache;
  cache.unique_abs_x_grid = BuildUniqueAbsXGrid(std::move(abs_x_values));
  cache.used_points = used_points;

  for (std::size_t histogram_index = 0; histogram_index < raw_samples.size();
       ++histogram_index) {
    const auto& samples = raw_samples[histogram_index];
    auto& histogram_cache = cache.histogram_caches[histogram_index];
    histogram_cache.abs_x_indices.reserve(samples.size());
    histogram_cache.y_values.reserve(samples.size());
    histogram_cache.sigma_values.reserve(samples.size());

    for (const RawProjectionSample& sample : samples) {
      histogram_cache.abs_x_indices.push_back(
          FindGridIndex(cache.unique_abs_x_grid, sample.abs_x));
      histogram_cache.y_values.push_back(sample.y);
      histogram_cache.sigma_values.push_back(sample.sigma);
    }
  }

  return cache;
}

ProjectionFitResult BuildFailedResult(const std::string& projection_name) {
  ProjectionFitResult result;
  result.projection_name = projection_name;
  result.r2 = std::numeric_limits<double>::quiet_NaN();
  result.r2_error = std::numeric_limits<double>::quiet_NaN();
  result.success = false;
  return result;
}

ProjectionFitResult BuildSuccessfulResult(const std::string& projection_name,
                                          const double r2,
                                          const double r2_error) {
  ProjectionFitResult result;
  result.projection_name = projection_name;
  result.amplitude = 1.0;
  result.r2 = r2;
  result.r2_error = r2_error;
  result.success = true;
  return result;
}

ProjectionFitResult BuildSuccessfulCentralValueOnlyResult(
    const std::string& projection_name,
    const double r2) {
  ProjectionFitResult result;
  result.projection_name = projection_name;
  result.amplitude = 1.0;
  result.r2 = r2;
  result.r2_error = std::numeric_limits<double>::quiet_NaN();
  result.success = true;
  return result;
}

ProjectionFitResult BuildNamedResult(const std::string& projection_name,
                                     const double r2,
                                     const bool error_valid,
                                     const double r2_error,
                                     const bool fail_result_when_error_invalid) {
  if (error_valid) {
    return BuildSuccessfulResult(projection_name, r2, r2_error);
  }

  if (fail_result_when_error_invalid) {
    return BuildFailedResult(projection_name);
  }

  return BuildSuccessfulCentralValueOnlyResult(projection_name, r2);
}

ProjectionFitResult BuildDirectionalResult(
    const std::string& projection_name,
    const double r2,
    const bool error_valid,
    const double r2_error,
    const ProjectionFitConfig& fit_config) {
  return BuildNamedResult(projection_name,
                          r2,
                          error_valid,
                          r2_error,
                          fit_config.fail_directional_results_when_error_invalid);
}

struct ScalarFitValue {
  double value = std::numeric_limits<double>::quiet_NaN();
  double error = std::numeric_limits<double>::quiet_NaN();
  bool success = false;
  bool error_valid = false;
};

ScalarFitValue BuildScalarFitValue(const double value,
                                   const bool error_valid,
                                   const double error,
                                   const bool fail_result_when_error_invalid) {
  ScalarFitValue result;
  if (error_valid) {
    result.value = value;
    result.error = error;
    result.success = true;
    result.error_valid = true;
    return result;
  }

  if (fail_result_when_error_invalid) {
    return result;
  }

  result.value = value;
  result.success = true;
  return result;
}

ScalarFitValue BuildAlphaResult(const double alpha,
                                const bool error_valid,
                                const double alpha_error,
                                const ProjectionFitConfig& fit_config) {
  return BuildScalarFitValue(alpha,
                             error_valid,
                             alpha_error,
                             fit_config.fail_alpha_result_when_error_invalid);
}

ProjectionFitResult BuildHbtResult(const std::string& projection_name,
                                   const double r2,
                                   const bool error_valid,
                                   const double r2_error,
                                   const ProjectionFitConfig& fit_config) {
  return BuildNamedResult(projection_name,
                          r2,
                          error_valid,
                          r2_error,
                          fit_config.fail_hbt_results_when_error_invalid);
}

bool IsValidLevyParameters(const double alpha, const double scale) {
  return alpha >= kMinimumAlpha && alpha <= kMaximumAlpha && scale > 0.0;
}

std::array<double, kLevyIntegrationPointCount> BuildFixedWeightedAlphaKernel(
    const double alpha) {
  std::array<double, kLevyIntegrationPointCount> weighted_alpha_kernel{};
  const FixedIntegrationWorkspace& workspace = GetFixedIntegrationWorkspace();

  for (std::size_t index = 0; index < kLevyIntegrationPointCount; ++index) {
    weighted_alpha_kernel[index] =
        workspace.simpson_weights[index] *
        std::exp(-0.5 * std::pow(workspace.u_grid[index], alpha));
  }

  return weighted_alpha_kernel;
}

double EvaluateFixedValueWithKernel(
    const double abs_x,
    const double scale,
    const std::array<double, kLevyIntegrationPointCount>& weighted_alpha_kernel) {
  const FixedIntegrationWorkspace& workspace = GetFixedIntegrationWorkspace();
  const double beta = abs_x / scale;

  double weighted_sum = 0.0;
  for (std::size_t index = 0; index < kLevyIntegrationPointCount; ++index) {
    weighted_sum +=
        weighted_alpha_kernel[index] * std::cos(beta * workspace.u_grid[index]);
  }

  return workspace.simpson_scale * weighted_sum * workspace.inverse_pi / scale;
}

std::vector<double> EvaluateFixedCurveWithKernel(
    const std::vector<double>& unique_abs_x_grid,
    const double scale,
    const std::array<double, kLevyIntegrationPointCount>& weighted_alpha_kernel) {
  std::vector<double> curve(unique_abs_x_grid.size(), 0.0);
  if (scale <= 0.0) {
    return curve;
  }

  for (std::size_t index = 0; index < unique_abs_x_grid.size(); ++index) {
    curve[index] =
        EvaluateFixedValueWithKernel(unique_abs_x_grid[index],
                                     scale,
                                     weighted_alpha_kernel);
  }

  return curve;
}

double SimpsonEstimate(const double interval_width,
                       const double left,
                       const double middle,
                       const double right) {
  return interval_width * (left + 4.0 * middle + right) / 6.0;
}

template <typename Integrand>
double AdaptiveSimpsonRecursive(const Integrand& integrand,
                                const double left_edge,
                                const double right_edge,
                                const double left_value,
                                const double middle_value,
                                const double right_value,
                                const double whole_interval_estimate,
                                const double absolute_tolerance,
                                const int depth_remaining) {
  const double midpoint = 0.5 * (left_edge + right_edge);
  const double left_midpoint = 0.5 * (left_edge + midpoint);
  const double right_midpoint = 0.5 * (midpoint + right_edge);
  const double left_midpoint_value = integrand(left_midpoint);
  const double right_midpoint_value = integrand(right_midpoint);

  const double left_estimate = SimpsonEstimate(midpoint - left_edge,
                                               left_value,
                                               left_midpoint_value,
                                               middle_value);
  const double right_estimate = SimpsonEstimate(right_edge - midpoint,
                                                middle_value,
                                                right_midpoint_value,
                                                right_value);
  const double refined_estimate = left_estimate + right_estimate;
  const double error_estimate =
      std::abs(refined_estimate - whole_interval_estimate);
  const double local_tolerance =
      std::max(absolute_tolerance,
               kAdaptiveRelativeTolerance * std::abs(refined_estimate));

  if (depth_remaining <= 0 || error_estimate <= 15.0 * local_tolerance) {
    return refined_estimate +
           (refined_estimate - whole_interval_estimate) / 15.0;
  }

  return AdaptiveSimpsonRecursive(integrand,
                                  left_edge,
                                  midpoint,
                                  left_value,
                                  left_midpoint_value,
                                  middle_value,
                                  left_estimate,
                                  0.5 * absolute_tolerance,
                                  depth_remaining - 1) +
         AdaptiveSimpsonRecursive(integrand,
                                  midpoint,
                                  right_edge,
                                  middle_value,
                                  right_midpoint_value,
                                  right_value,
                                  right_estimate,
                                  0.5 * absolute_tolerance,
                                  depth_remaining - 1);
}

double LevyStable1DValueFixed(const double rho,
                              const double alpha,
                              const double scale) {
  if (!IsValidLevyParameters(alpha, scale)) {
    return 0.0;
  }

  const auto weighted_alpha_kernel = BuildFixedWeightedAlphaKernel(alpha);
  return EvaluateFixedValueWithKernel(std::abs(rho), scale, weighted_alpha_kernel);
}

double LevyStable1DValueAdaptive(const double rho,
                                 const double alpha,
                                 const double scale) {
  if (!IsValidLevyParameters(alpha, scale)) {
    return 0.0;
  }

  const double beta = std::abs(rho) / scale;
  const auto integrand = [beta, alpha](const double u) {
    return std::cos(beta * u) * std::exp(-0.5 * std::pow(u, alpha));
  };

  const double left_edge = 0.0;
  const double right_edge = kLevyIntegrationUpperBound;
  const double midpoint = 0.5 * (left_edge + right_edge);
  const double left_value = integrand(left_edge);
  const double middle_value = integrand(midpoint);
  const double right_value = integrand(right_edge);
  const double whole_interval_estimate = SimpsonEstimate(right_edge - left_edge,
                                                         left_value,
                                                         middle_value,
                                                         right_value);

  const double integral = AdaptiveSimpsonRecursive(integrand,
                                                   left_edge,
                                                   right_edge,
                                                   left_value,
                                                   middle_value,
                                                   right_value,
                                                   whole_interval_estimate,
                                                   kAdaptiveAbsoluteTolerance,
                                                   kAdaptiveMaxDepth);
  return integral / (kPi * scale);
}

std::vector<double> EvaluateAdaptiveCurve(const std::vector<double>& unique_abs_x_grid,
                                          const double alpha,
                                          const double scale) {
  std::vector<double> curve(unique_abs_x_grid.size(), 0.0);
  if (!IsValidLevyParameters(alpha, scale)) {
    return curve;
  }

  for (std::size_t index = 0; index < unique_abs_x_grid.size(); ++index) {
    curve[index] = LevyStable1DValueAdaptive(unique_abs_x_grid[index],
                                             alpha,
                                             scale);
  }

  return curve;
}

std::vector<double> EvaluateLevyCurveOnGrid(
    const std::vector<double>& unique_abs_x_grid,
    const double alpha,
    const double scale,
    const ProjectionFitConfig& fit_config,
    const std::array<double, kLevyIntegrationPointCount>* fixed_weighted_alpha_kernel =
        nullptr) {
  if (fit_config.use_adaptive_integration) {
    return EvaluateAdaptiveCurve(unique_abs_x_grid, alpha, scale);
  }

  if (!IsValidLevyParameters(alpha, scale)) {
    return std::vector<double>(unique_abs_x_grid.size(), 0.0);
  }

  if (fixed_weighted_alpha_kernel != nullptr) {
    return EvaluateFixedCurveWithKernel(unique_abs_x_grid,
                                        scale,
                                        *fixed_weighted_alpha_kernel);
  }

  const auto weighted_alpha_kernel = BuildFixedWeightedAlphaKernel(alpha);
  return EvaluateFixedCurveWithKernel(unique_abs_x_grid,
                                      scale,
                                      weighted_alpha_kernel);
}

std::size_t FindMatchingScaleIndex(
    const std::array<double, 6>& unique_scales,
    const std::size_t unique_scale_count,
    const double scale) {
  for (std::size_t index = 0; index < unique_scale_count; ++index) {
    if (NearlyEqual(unique_scales[index], scale)) {
      return index;
    }
  }

  return unique_scale_count;
}

ProjectedCurveSet BuildProjectedCurves(
    const std::vector<double>& unique_abs_x_grid,
    const double alpha,
    const std::array<double, 6>& projected_r2,
    const ProjectionFitConfig& fit_config) {
  ProjectedCurveSet curve_set;
  std::array<double, kLevyIntegrationPointCount> fixed_weighted_alpha_kernel{};
  const std::array<double, kLevyIntegrationPointCount>* kernel_pointer = nullptr;

  if (!fit_config.use_adaptive_integration) {
    fixed_weighted_alpha_kernel = BuildFixedWeightedAlphaKernel(alpha);
    kernel_pointer = &fixed_weighted_alpha_kernel;
  }

  for (std::size_t direction_index = 0; direction_index < projected_r2.size();
       ++direction_index) {
    const double scale = std::sqrt(projected_r2[direction_index]);
    const std::size_t existing_curve_index =
        FindMatchingScaleIndex(curve_set.unique_scales,
                               curve_set.unique_scale_count,
                               scale);

    if (existing_curve_index < curve_set.unique_scale_count) {
      curve_set.direction_curve_indices[direction_index] = existing_curve_index;
      continue;
    }

    const std::size_t curve_index = curve_set.unique_scale_count;
    curve_set.unique_scales[curve_index] = scale;
    curve_set.curves[curve_index] =
        EvaluateLevyCurveOnGrid(unique_abs_x_grid,
                                alpha,
                                scale,
                                fit_config,
                                kernel_pointer);
    curve_set.direction_curve_indices[direction_index] = curve_index;
    ++curve_set.unique_scale_count;
  }

  return curve_set;
}

double LevyStableTf1EvaluatorFixed(double* x, double* parameters) {
  return LevyStable1DValueFixed(x[0], parameters[0], parameters[1]);
}

double LevyStableTf1EvaluatorAdaptive(double* x, double* parameters) {
  return LevyStable1DValueAdaptive(x[0], parameters[0], parameters[1]);
}

std::array<double, 6> ComputeProjectedR2(
    const std::array<double, 6>& hbt_radii_r2) {
  return {{
      hbt_radii_r2[0],
      hbt_radii_r2[1],
      hbt_radii_r2[2],
      0.5 * (hbt_radii_r2[0] + hbt_radii_r2[1]) + hbt_radii_r2[3],
      0.5 * (hbt_radii_r2[0] + hbt_radii_r2[2]) + hbt_radii_r2[4],
      0.5 * (hbt_radii_r2[1] + hbt_radii_r2[2]) + hbt_radii_r2[5],
  }};
}

bool IsFiniteMatrix(const CovarianceMatrix6& matrix) {
  for (const auto& row : matrix) {
    for (const double value : row) {
      if (!std::isfinite(value)) {
        return false;
      }
    }
  }

  return true;
}

bool HasNonNegativeDiagonal(const CovarianceMatrix6& matrix) {
  for (std::size_t index = 0; index < matrix.size(); ++index) {
    const double diagonal = matrix[index][index];
    if (!std::isfinite(diagonal) || diagonal < 0.0) {
      return false;
    }
  }

  return true;
}

CovarianceMatrix6 PropagateDirectionalCovariance(
    const CovarianceMatrix6& hbt_covariance) {
  CovarianceMatrix6 directional_covariance{};

  for (std::size_t row = 0; row < kDirectionalCount; ++row) {
    for (std::size_t column = 0; column < kDirectionalCount; ++column) {
      double propagated_value = 0.0;
      for (std::size_t left_index = 0; left_index < kDirectionalCount;
           ++left_index) {
        for (std::size_t right_index = 0; right_index < kDirectionalCount;
             ++right_index) {
          propagated_value +=
              kDirectionalTransformMatrix[row][left_index] *
              hbt_covariance[left_index][right_index] *
              kDirectionalTransformMatrix[column][right_index];
        }
      }
      directional_covariance[row][column] = propagated_value;
    }
  }

  return directional_covariance;
}

CovariancePropagationResult BuildCovariancePropagationResult(
    const ROOT::Math::Minimizer& minimizer,
    const bool minimize_success,
    const bool hesse_success,
    const ProjectionFitConfig& fit_config) {
  CovariancePropagationResult result;
  result.covariance_status = minimizer.CovMatrixStatus();
  if (!minimize_success || !hesse_success ||
      minimizer.NDim() < kFitParameterCount ||
      result.covariance_status <= 0) {
    return result;
  }

  std::array<double, kFitParameterCount * kFitParameterCount> full_covariance{};
  if (!minimizer.GetCovMatrix(full_covariance.data())) {
    return result;
  }

  for (std::size_t row = 0; row < kDirectionalCount; ++row) {
    for (std::size_t column = 0; column < kDirectionalCount; ++column) {
      result.hbt_covariance[row][column] =
          full_covariance[(row + kHbtParameterOffset) * kFitParameterCount +
                          (column + kHbtParameterOffset)];
    }
  }

  if (!IsFiniteMatrix(result.hbt_covariance) ||
      !HasNonNegativeDiagonal(result.hbt_covariance)) {
    return result;
  }

  result.covariance_valid =
      (result.covariance_status == 3) ||
      (fit_config.accept_forced_posdef_covariance_as_valid &&
       result.covariance_status == 2);

  result.directional_covariance =
      PropagateDirectionalCovariance(result.hbt_covariance);
  if (!IsFiniteMatrix(result.directional_covariance) ||
      !HasNonNegativeDiagonal(result.directional_covariance)) {
    result.directional_covariance = {};
    return result;
  }

  for (std::size_t index = 0; index < kDirectionalCount; ++index) {
    result.directional_errors[index] =
        std::sqrt(result.directional_covariance[index][index]);
  }
  result.directional_valid = result.covariance_valid;
  return result;
}

std::unique_ptr<TH1D> MakeNormalizedClone(const TH1D& histogram,
                                          const std::string& clone_name) {
  auto clone =
      std::unique_ptr<TH1D>(static_cast<TH1D*>(histogram.Clone(clone_name.c_str())));
  clone->SetDirectory(nullptr);

  const double normalization = clone->Integral("width");
  if (normalization > 0.0) {
    clone->Scale(1.0 / normalization);
  }

  return clone;
}

struct SimultaneousLevyChi2 {
  std::shared_ptr<const ProjectionFitCache> cache;
  ProjectionFitConfig fit_config;

  double operator()(const double* parameters) const {
    if (cache == nullptr || cache->used_points <= 0) {
      return kLargePenalty;
    }

    const double alpha = parameters[0];
    if (alpha < kMinimumAlpha || alpha > kMaximumAlpha) {
      return kLargePenalty;
    }

    const std::array<double, 6> hbt_radii_r2 = {
        parameters[1], parameters[2], parameters[3],
        parameters[4], parameters[5], parameters[6]};
    const std::array<double, 6> projected_r2 = ComputeProjectedR2(hbt_radii_r2);

    for (const double projected_value : projected_r2) {
      if (projected_value <= kMinimumPositiveR2) {
        return kLargePenalty;
      }
    }

    const ProjectedCurveSet curve_set =
        BuildProjectedCurves(cache->unique_abs_x_grid,
                             alpha,
                             projected_r2,
                             fit_config);
    double chi2 = 0.0;
    for (std::size_t histogram_index = 0;
         histogram_index < cache->histogram_caches.size();
         ++histogram_index) {
      const ProjectionHistogramCache& histogram_cache =
          cache->histogram_caches[histogram_index];
      const std::vector<double>& model_curve =
          curve_set.curves[curve_set.direction_curve_indices[histogram_index]];

      for (std::size_t sample_index = 0;
           sample_index < histogram_cache.y_values.size();
           ++sample_index) {
        const double model =
            model_curve[histogram_cache.abs_x_indices[sample_index]];
        const double residual =
            (histogram_cache.y_values[sample_index] - model) /
            histogram_cache.sigma_values[sample_index];
        chi2 += residual * residual;
      }
    }

    return chi2;
  }
};

}  // namespace

SliceFitProducts FitSliceHistograms(
    const std::string& slice_name,
    std::array<std::unique_ptr<TH1D>, 6>& projection_histograms,
    const std::array<ProjectionDefinition, 6>& directions,
    ProjectionFitConfig fit_config) {
  SliceFitProducts products;

  std::array<std::unique_ptr<TH1D>, 6> normalized_histograms;

  bool has_enough_statistics = true;
  for (std::size_t index = 0; index < directions.size(); ++index) {
    TH1D* histogram = projection_histograms[index].get();
    if (histogram == nullptr || histogram->GetEntries() < 5.0 ||
        histogram->Integral("width") <= 0.0) {
      has_enough_statistics = false;
      break;
    }

    normalized_histograms[index] = MakeNormalizedClone(
        *histogram, histogram->GetName() + std::string("_normalized"));
  }

  if (has_enough_statistics) {
    const auto fit_cache =
        std::make_shared<ProjectionFitCache>(BuildProjectionFitCache(
            normalized_histograms));

    std::array<double, 6> initial_hbt_radii_r2 = {0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
    for (std::size_t index = 0; index < 3U; ++index) {
      const double rms = projection_histograms[index]->GetRMS();
      initial_hbt_radii_r2[index] =
          std::max(0.5 * rms * rms, 0.05);
    }

    const double cross_limit =
        5.0 * std::max({initial_hbt_radii_r2[0],
                        initial_hbt_radii_r2[1],
                        initial_hbt_radii_r2[2],
                        1.0});

    SimultaneousLevyChi2 chi2{fit_cache, fit_config};
    ROOT::Math::Functor objective(chi2, 7);
    std::unique_ptr<ROOT::Math::Minimizer> minimizer(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

    if (minimizer != nullptr) {
      minimizer->SetMaxFunctionCalls(20000);
      minimizer->SetMaxIterations(5000);
      minimizer->SetTolerance(1.0e-3);
      minimizer->SetFunction(objective);

      minimizer->SetLimitedVariable(0, "alpha", 1.8, 0.05, kMinimumAlpha, kMaximumAlpha);
      minimizer->SetLimitedVariable(1,
                                    "Rout2",
                                    initial_hbt_radii_r2[0],
                                    0.05,
                                    kMinimumPositiveR2,
                                    1.0e3);
      minimizer->SetLimitedVariable(2,
                                    "Rside2",
                                    initial_hbt_radii_r2[1],
                                    0.05,
                                    kMinimumPositiveR2,
                                    1.0e3);
      minimizer->SetLimitedVariable(3,
                                    "Rlong2",
                                    initial_hbt_radii_r2[2],
                                    0.05,
                                    kMinimumPositiveR2,
                                    1.0e3);
      minimizer->SetLimitedVariable(4, "Ros2", 0.0, 0.02, -cross_limit, cross_limit);
      minimizer->SetLimitedVariable(5, "Rol2", 0.0, 0.02, -cross_limit, cross_limit);
      minimizer->SetLimitedVariable(6, "Rsl2", 0.0, 0.02, -cross_limit, cross_limit);

      const bool minimize_success = minimizer->Minimize();
      const bool hesse_success = minimize_success ? minimizer->Hesse() : false;
      if (minimize_success) {
        const double* best_fit_parameters = minimizer->X();
        const double* best_fit_errors = minimizer->Errors();
        const CovariancePropagationResult covariance_result =
            BuildCovariancePropagationResult(*minimizer,
                                             minimize_success,
                                             hesse_success,
                                             fit_config);

        const std::array<double, 6> hbt_radii_r2 = {
            best_fit_parameters[1], best_fit_parameters[2], best_fit_parameters[3],
            best_fit_parameters[4], best_fit_parameters[5], best_fit_parameters[6]};
        const std::array<double, 6> hbt_radii_r2_errors = {
            best_fit_errors[1], best_fit_errors[2], best_fit_errors[3],
            best_fit_errors[4], best_fit_errors[5], best_fit_errors[6]};
        const std::array<double, 6> projected_r2 = ComputeProjectedR2(hbt_radii_r2);
        std::array<double, 6> projected_r2_errors{};
        if (covariance_result.directional_valid) {
          projected_r2_errors = covariance_result.directional_errors;
        }
        const ScalarFitValue alpha_result =
            BuildAlphaResult(best_fit_parameters[0],
                             covariance_result.covariance_valid,
                             best_fit_errors[0],
                             fit_config);

        products.alpha = alpha_result.value;
        products.alpha_error = alpha_result.error;
        products.alpha_success = alpha_result.success;
        products.alpha_error_valid = alpha_result.error_valid;
        products.covariance_status = covariance_result.covariance_status;
        products.hbt_covariance = covariance_result.hbt_covariance;
        products.directional_covariance = covariance_result.directional_covariance;
        products.covariance_valid = covariance_result.covariance_valid;

        const std::array<std::string, 6> hbt_names = {
            "Rout2", "Rside2", "Rlong2", "Ros2", "Rol2", "Rsl2"};
        for (std::size_t index = 0; index < hbt_names.size(); ++index) {
          products.hbt_error_valid[index] = covariance_result.covariance_valid;
          products.hbt_radii_results[index] = BuildHbtResult(
              hbt_names[index],
              hbt_radii_r2[index],
              covariance_result.covariance_valid,
              hbt_radii_r2_errors[index],
              fit_config);
        }

        for (std::size_t index = 0; index < directions.size(); ++index) {
          products.directional_error_valid[index] =
              covariance_result.directional_valid;
          products.directional_results[index] = BuildDirectionalResult(
              directions[index].short_name,
              projected_r2[index],
              covariance_result.directional_valid,
              projected_r2_errors[index],
              fit_config);

          TH1D* normalized_histogram = normalized_histograms[index].get();
          const std::string fit_name = "fit_" + directions[index].short_name;
          products.fit_functions[index] = std::make_unique<TF1>(
              fit_name.c_str(),
              fit_config.use_adaptive_integration
                  ? LevyStableTf1EvaluatorAdaptive
                  : LevyStableTf1EvaluatorFixed,
              normalized_histogram->GetXaxis()->GetXmin(),
              normalized_histogram->GetXaxis()->GetXmax(),
              2);
          products.fit_functions[index]->SetParNames("alpha", "R");
          products.fit_functions[index]->SetParameters(best_fit_parameters[0],
                                                      std::sqrt(projected_r2[index]));
          products.fit_functions[index]->SetLineColor(kRed + 1);
          products.fit_functions[index]->SetNpx(400);
        }
      }
    }
  }

  for (std::size_t index = 0; index < directions.size(); ++index) {
    if (!products.directional_results[index].success) {
      products.directional_results[index] =
          BuildFailedResult(directions[index].short_name);
    }
    if (!products.hbt_radii_results[index].success) {
      const std::array<std::string, 6> hbt_names = {
          "Rout2", "Rside2", "Rlong2", "Ros2", "Rol2", "Rsl2"};
      products.hbt_radii_results[index] = BuildFailedResult(hbt_names[index]);
    }
  }

  const std::string canvas_name = slice_name + "_fit_canvas";
  products.fit_canvas =
      std::make_unique<TCanvas>(canvas_name.c_str(), slice_name.c_str(), 1800, 900);
  products.fit_canvas->Divide(3, 2);

  for (std::size_t index = 0; index < directions.size(); ++index) {
    products.fit_canvas->cd(static_cast<int>(index + 1U));
    if (normalized_histograms[index] == nullptr) {
      continue;
    }

    normalized_histograms[index]->SetTitle(directions[index].title.c_str());
    normalized_histograms[index]->GetYaxis()->SetTitle("Normalized counts");
    normalized_histograms[index]->DrawCopy("E");
    if (products.directional_results[index].success &&
        products.fit_functions[index] != nullptr) {
      products.fit_functions[index]->Draw("SAME");
    }
  }

  return products;
}

}  // namespace femto3d
