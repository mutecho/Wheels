#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "exp_femto_3d/Types.h"

namespace exp_femto_3d::profile_likelihood {

  enum class PointStatus {
    kValid,
    kNonconverged,
    kMinimizerError,
    kModelDomainInvalid,
    kObjectiveInvalid,
    kNoValidAttempt,
  };

  enum class SeedOrigin {
    kNominal,
    kForwardNeighbor,
    kReverseNeighbor,
  };

  struct MinimizationResult {
    std::vector<double> values;
    std::vector<double> errors;
    double objective = std::numeric_limits<double>::quiet_NaN();
    double edm = std::numeric_limits<double>::quiet_NaN();
    int migrad_status = -1;
    int hesse_status = -1;
    int minuit_istat = -1;
    bool model_domain_valid = false;
    bool objective_valid = false;
    bool minimizer_error = false;
    std::size_t fcn_calls = 0;
    double total_wall_ms = std::numeric_limits<double>::quiet_NaN();
    double migrad_wall_ms = std::numeric_limits<double>::quiet_NaN();
    double hesse_wall_ms = std::numeric_limits<double>::quiet_NaN();
    bool hesse_ran = false;
    bool parameter_errors_valid = false;
    std::vector<int> at_lower_bound;
    std::vector<int> at_upper_bound;
  };

  struct AttemptRecord {
    int attempt_index = 0;
    int point_index = 0;
    int stage = 0;
    int ix = 0;
    int iy = -1;
    SeedOrigin seed_origin = SeedOrigin::kNominal;
    std::vector<double> coordinates;
    MinimizationResult result;
    PointStatus status = PointStatus::kNoValidAttempt;
  };

  struct PointRecord {
    int point_index = 0;
    int stage = 0;
    int ix = 0;
    int iy = -1;
    std::vector<double> coordinates;
    MinimizationResult winner;
    PointStatus status = PointStatus::kNoValidAttempt;
    SeedOrigin winner_seed = SeedOrigin::kNominal;
    int attempt_count = 0;
    int valid_attempt_count = 0;
    double objective_spread = std::numeric_limits<double>::quiet_NaN();
    double likelihood_slice_objective = std::numeric_limits<double>::quiet_NaN();
    bool likelihood_slice_objective_valid = false;
  };

  struct ScanResult {
    std::vector<PointRecord> points;
    std::vector<AttemptRecord> attempts;
    bool refinement_performed = false;
    std::string refinement_reason;
  };

  struct SliceEvaluation {
    SliceEvaluation() = default;
    // Preserve source compatibility with existing tests and callers that return only an objective scalar.
    SliceEvaluation(const double objective_value)
        : objective(objective_value),
          model_domain_valid(std::isfinite(objective_value)),
          objective_valid(std::isfinite(objective_value)) {}

    double objective = std::numeric_limits<double>::quiet_NaN();
    bool model_domain_valid = false;
    bool objective_valid = false;
  };

  struct ReferenceMinimum {
    double objective = std::numeric_limits<double>::quiet_NaN();
    std::string source = "none";
  };

  using MinimizeFunction = std::function<MinimizationResult(const std::vector<double>& seed,
                                                             const std::vector<int>& fixed_indices,
                                                             const std::vector<double>& fixed_values)>;
  using SliceFunction = std::function<SliceEvaluation(const std::vector<int>& fixed_indices,
                                                      const std::vector<double>& fixed_values)>;

  /**
   * @brief Serial deterministic driver for one 1D/2D profile grid.
   *
   * The driver knows nothing about ROOT or the PML evaluator.  It always attempts a
   * nominal seed and, when enabled, axis-adjacent forward and reverse warm starts.
   * A point is valid only when its supplied minimizer result says that the objective
   * and physical model domain are valid and MIGRAD succeeded.
   */
  ScanResult RunProfileScan(const ProfileScanConfig& config,
                            const std::vector<int>& target_indices,
                            const std::vector<double>& effective_min,
                            const std::vector<double>& effective_max,
                            const std::vector<double>& nominal_values,
                            ProfileRetryStrategy retry_strategy,
                            const MinimizeFunction& minimize,
                            const SliceFunction& evaluate_slice);

  // Select one reference only after every scan (including refinement) is complete.
  ReferenceMinimum SelectGlobalReference(
      double nominal_objective,
      bool nominal_valid,
      const std::vector<std::pair<std::string, const ScanResult*>>& scans);

  [[nodiscard]] bool IsValid(const MinimizationResult& result);
  [[nodiscard]] PointStatus StatusFor(const MinimizationResult& result);
  [[nodiscard]] const char* ToString(PointStatus status);
  [[nodiscard]] const char* ToString(SeedOrigin origin);

}  // namespace exp_femto_3d::profile_likelihood
