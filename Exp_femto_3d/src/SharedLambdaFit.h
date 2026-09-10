#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace exp_femto_3d::shared_lambda {

  /** One local fit parameter before shared/free parameters are flattened for TMinuit. */
  struct LocalParameter {
    std::string name;
    double initial = 0.0;
    double step = 0.01;
    double lower = 0.0;
    double upper = 0.0;
    bool fixed = false;
  };

  /** The stable free-parameter ordering used by both Minuit and the saved covariance matrix. */
  struct ParameterLayout {
    std::vector<std::string> labels;
    std::vector<std::string> slice_ids;
    std::vector<int> local_indices;
    std::vector<double> initial;
    std::vector<double> steps;
    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<int> has_limits;
    std::vector<std::vector<int>> local_to_global;
  };

  /**
   * Build a layout with lambda (local index 1) represented exactly once.
   * Fixed local parameters map to -1; every other free nuisance is owned by its slice.
   */
  ParameterLayout BuildParameterLayout(const std::vector<std::string> &slice_ids,
                                       const std::vector<std::vector<LocalParameter>> &parameters,
                                       int lambda_local_index = 1);

  /** Expand a free vector back to one complete local vector, restoring fixed values. */
  std::vector<double> ExpandLocalValues(const ParameterLayout &layout,
                                        std::size_t member_index,
                                        const std::vector<LocalParameter> &parameters,
                                        const double *free_values);

  /** True only for a finite, symmetric, strictly positive-definite row-major matrix. */
  bool IsPositiveDefinite(const std::vector<double> &matrix,
                          std::size_t dimension,
                          double relative_tolerance = 1.0e-12);

  /** GLS result for y(phi)=A+2B*cos(2phi) or A+2B*sin(2phi). */
  struct HarmonicResult {
    double intercept = 0.0;
    double harmonic = 0.0;
    double intercept_variance = 0.0;
    double harmonic_variance = 0.0;
    double covariance = 0.0;
    bool valid = false;
    std::string failure_reason;
  };

  /**
   * Fit a second harmonic with the complete cross-phi covariance.
   * The covariance is solved through Cholesky factors; no explicit inverse is formed.
   */
  HarmonicResult FitSecondHarmonicGLS(const std::vector<double> &phi,
                                     const std::vector<double> &values,
                                     const std::vector<double> &covariance,
                                     bool sine_form);

  struct ConstantResult {
    double value = 0.0;
    double variance = 0.0;
    bool valid = false;
    std::string failure_reason;
  };

  /** Fit one common value with a complete covariance, without repeated-point shortcuts. */
  ConstantResult FitConstantGLS(const std::vector<double> &values,
                                const std::vector<double> &covariance);

}  // namespace exp_femto_3d::shared_lambda
