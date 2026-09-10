#include "SharedLambdaFit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace exp_femto_3d::shared_lambda {

  namespace {

    bool Cholesky(const std::vector<double> &matrix,
                  const std::size_t dimension,
                  std::vector<double> &lower,
                  const double relative_tolerance) {
      if (matrix.size() != dimension * dimension || dimension == 0U) return false;
      lower.assign(matrix.size(), 0.0);
      double diagonal_scale = 0.0;
      for (std::size_t row = 0; row < dimension; ++row) {
        const double diagonal = matrix[row * dimension + row];
        if (!std::isfinite(diagonal)) return false;
        diagonal_scale = std::max(diagonal_scale, std::abs(diagonal));
        for (std::size_t column = 0; column < dimension; ++column) {
          const double value = matrix[row * dimension + column];
          const double transpose = matrix[column * dimension + row];
          if (!std::isfinite(value) || !std::isfinite(transpose)
              || std::abs(value - transpose) > relative_tolerance * (1.0 + std::abs(value))) {
            return false;
          }
        }
      }
      const double floor = relative_tolerance * std::max(1.0, diagonal_scale);
      for (std::size_t row = 0; row < dimension; ++row) {
        for (std::size_t column = 0; column <= row; ++column) {
          double value = matrix[row * dimension + column];
          for (std::size_t inner = 0; inner < column; ++inner) {
            value -= lower[row * dimension + inner] * lower[column * dimension + inner];
          }
          if (row == column) {
            if (!std::isfinite(value) || value <= floor) return false;
            lower[row * dimension + column] = std::sqrt(value);
          } else {
            lower[row * dimension + column] = value / lower[column * dimension + column];
          }
        }
      }
      return true;
    }

    std::vector<double> SolveCholesky(const std::vector<double> &lower,
                                      const std::size_t dimension,
                                      const std::vector<double> &rhs) {
      std::vector<double> intermediate(dimension, 0.0);
      std::vector<double> solution(dimension, 0.0);
      for (std::size_t row = 0; row < dimension; ++row) {
        double value = rhs[row];
        for (std::size_t column = 0; column < row; ++column) {
          value -= lower[row * dimension + column] * intermediate[column];
        }
        intermediate[row] = value / lower[row * dimension + row];
      }
      for (std::size_t reverse = 0; reverse < dimension; ++reverse) {
        const std::size_t row = dimension - reverse - 1U;
        double value = intermediate[row];
        for (std::size_t column = row + 1U; column < dimension; ++column) {
          value -= lower[column * dimension + row] * solution[column];
        }
        solution[row] = value / lower[row * dimension + row];
      }
      return solution;
    }

  }  // namespace

  ParameterLayout BuildParameterLayout(const std::vector<std::string> &slice_ids,
                                       const std::vector<std::vector<LocalParameter>> &parameters,
                                       const int lambda_local_index) {
    if (slice_ids.empty() || slice_ids.size() != parameters.size()) {
      throw std::invalid_argument("shared-lambda layout requires one parameter list per non-empty member list");
    }
    ParameterLayout layout;
    layout.local_to_global.resize(slice_ids.size());
    auto append = [&](const std::string &label,
                            const std::string &slice_id,
                            const int local_index,
                            const LocalParameter &parameter) {
      layout.labels.push_back(label);
      layout.slice_ids.push_back(slice_id);
      layout.local_indices.push_back(local_index);
      layout.initial.push_back(parameter.initial);
      layout.steps.push_back(parameter.step);
      layout.lower.push_back(parameter.lower);
      layout.upper.push_back(parameter.upper);
      layout.has_limits.push_back(parameter.lower < parameter.upper ? 1 : 0);
      return static_cast<int>(layout.labels.size() - 1U);
    };

    if (lambda_local_index < 0
        || static_cast<std::size_t>(lambda_local_index) >= parameters.front().size()
        || parameters.front()[static_cast<std::size_t>(lambda_local_index)].fixed) {
      throw std::invalid_argument("shared lambda must be a free local parameter");
    }
    const LocalParameter &lambda_parameter = parameters.front()[static_cast<std::size_t>(lambda_local_index)];
    const int shared_lambda_index = append("lambda", "<shared>", lambda_local_index, lambda_parameter);
    for (std::size_t member = 0; member < parameters.size(); ++member) {
      layout.local_to_global[member].assign(parameters[member].size(), -1);
      for (std::size_t local = 0; local < parameters[member].size(); ++local) {
        const LocalParameter &parameter = parameters[member][local];
        if (parameter.fixed) continue;
        if (static_cast<int>(local) == lambda_local_index) {
          layout.local_to_global[member][local] = shared_lambda_index;
          continue;
        }
        layout.local_to_global[member][local] = append(
            slice_ids[member] + "::" + parameter.name, slice_ids[member], static_cast<int>(local), parameter);
      }
    }
    return layout;
  }

  std::vector<double> ExpandLocalValues(const ParameterLayout &layout,
                                        const std::size_t member_index,
                                        const std::vector<LocalParameter> &parameters,
                                        const double *free_values) {
    if (member_index >= layout.local_to_global.size()
        || layout.local_to_global[member_index].size() != parameters.size() || free_values == nullptr) {
      throw std::invalid_argument("invalid shared-lambda local expansion request");
    }
    std::vector<double> result(parameters.size(), 0.0);
    for (std::size_t local = 0; local < parameters.size(); ++local) {
      const int global = layout.local_to_global[member_index][local];
      result[local] = global >= 0 ? free_values[global] : parameters[local].initial;
    }
    return result;
  }

  bool IsPositiveDefinite(const std::vector<double> &matrix,
                          const std::size_t dimension,
                          const double relative_tolerance) {
    std::vector<double> lower;
    return Cholesky(matrix, dimension, lower, relative_tolerance);
  }

  HarmonicResult FitSecondHarmonicGLS(const std::vector<double> &phi,
                                     const std::vector<double> &values,
                                     const std::vector<double> &covariance,
                                     const bool sine_form) {
    HarmonicResult result;
    const std::size_t count = phi.size();
    if (count < 2U || values.size() != count || covariance.size() != count * count) {
      result.failure_reason = "insufficient points or inconsistent matrix dimensions";
      return result;
    }
    if (std::any_of(phi.begin(), phi.end(), [](const double value) { return !std::isfinite(value); })
        || std::any_of(values.begin(), values.end(), [](const double value) { return !std::isfinite(value); })) {
      result.failure_reason = "non-finite input point";
      return result;
    }
    std::vector<double> lower;
    if (!Cholesky(covariance, count, lower, 1.0e-12)) {
      result.failure_reason = "covariance is not finite positive definite";
      return result;
    }
    std::vector<double> column0(count, 1.0);
    std::vector<double> column1(count, 0.0);
    for (std::size_t index = 0; index < count; ++index) {
      // The fitted coefficient is B in A+2B*cos(2phi) or A+2B*sin(2phi).
      column1[index] = 2.0 * (sine_form ? std::sin(2.0 * phi[index]) : std::cos(2.0 * phi[index]));
    }
    const std::vector<double> vinv_x0 = SolveCholesky(lower, count, column0);
    const std::vector<double> vinv_x1 = SolveCholesky(lower, count, column1);
    const std::vector<double> vinv_y = SolveCholesky(lower, count, values);
    double normal00 = 0.0;
    double normal01 = 0.0;
    double normal11 = 0.0;
    double rhs0 = 0.0;
    double rhs1 = 0.0;
    for (std::size_t index = 0; index < count; ++index) {
      normal00 += column0[index] * vinv_x0[index];
      normal01 += column0[index] * vinv_x1[index];
      normal11 += column1[index] * vinv_x1[index];
      rhs0 += column0[index] * vinv_y[index];
      rhs1 += column1[index] * vinv_y[index];
    }
    const double determinant = normal00 * normal11 - normal01 * normal01;
    if (!std::isfinite(determinant) || determinant <= 1.0e-20) {
      result.failure_reason = "harmonic design matrix is rank deficient";
      return result;
    }
    result.intercept = (normal11 * rhs0 - normal01 * rhs1) / determinant;
    result.harmonic = (normal00 * rhs1 - normal01 * rhs0) / determinant;
    result.intercept_variance = normal11 / determinant;
    result.harmonic_variance = normal00 / determinant;
    result.covariance = -normal01 / determinant;
    result.valid = std::isfinite(result.intercept) && std::isfinite(result.harmonic)
                   && result.intercept_variance >= 0.0 && result.harmonic_variance >= 0.0;
    if (!result.valid) result.failure_reason = "non-finite generalized least-squares result";
    return result;
  }

  ConstantResult FitConstantGLS(const std::vector<double> &values,
                                const std::vector<double> &covariance) {
    ConstantResult result;
    const std::size_t count = values.size();
    if (count == 0U || covariance.size() != count * count
        || std::any_of(values.begin(), values.end(), [](const double value) { return !std::isfinite(value); })) {
      result.failure_reason = "empty, non-finite, or dimensionally inconsistent constant input";
      return result;
    }
    std::vector<double> lower;
    if (!Cholesky(covariance, count, lower, 1.0e-12)) {
      result.failure_reason = "covariance is not finite positive definite";
      return result;
    }
    const std::vector<double> ones(count, 1.0);
    const std::vector<double> vinv_ones = SolveCholesky(lower, count, ones);
    const std::vector<double> vinv_values = SolveCholesky(lower, count, values);
    double denominator = 0.0;
    double numerator = 0.0;
    for (std::size_t index = 0; index < count; ++index) {
      denominator += vinv_ones[index];
      numerator += vinv_values[index];
    }
    if (!std::isfinite(denominator) || denominator <= 0.0) {
      result.failure_reason = "constant design matrix is rank deficient";
      return result;
    }
    result.value = numerator / denominator;
    result.variance = 1.0 / denominator;
    result.valid = std::isfinite(result.value) && std::isfinite(result.variance) && result.variance >= 0.0;
    if (!result.valid) result.failure_reason = "non-finite constant GLS result";
    return result;
  }

}  // namespace exp_femto_3d::shared_lambda
