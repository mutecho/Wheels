#pragma once

#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

#include "ProfileLikelihood.h"

namespace exp_femto_3d::combined_profile {

  struct MemberAttempt {
    int stage = 0;
    int point_index = -1;
    int member_index = -1;
    double lambda = std::numeric_limits<double>::quiet_NaN();
    profile_likelihood::SeedOrigin seed_origin = profile_likelihood::SeedOrigin::kNominal;
    profile_likelihood::MinimizationResult result;
  };

  struct MemberPoint {
    int member_index = -1;
    profile_likelihood::MinimizationResult winner;
    profile_likelihood::PointStatus status = profile_likelihood::PointStatus::kNoValidAttempt;
    profile_likelihood::SeedOrigin winner_seed = profile_likelihood::SeedOrigin::kNominal;
    int attempt_count = 0;
    int valid_attempt_count = 0;
  };

  struct CombinedPoint {
    int stage = 0;
    int point_index = -1;
    double lambda = std::numeric_limits<double>::quiet_NaN();
    std::vector<MemberPoint> members;
    double objective = std::numeric_limits<double>::quiet_NaN();
    double delta_objective = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    std::string failure_reason;
  };

  struct CandidateInterval {
    double lower = std::numeric_limits<double>::quiet_NaN();
    double upper = std::numeric_limits<double>::quiet_NaN();
    bool touches_lower_bound = false;
    bool touches_upper_bound = false;
    bool contains_failure_gap = false;
  };

  struct StageRecord {
    int stage = 0;
    std::vector<double> coordinates;
    std::vector<CandidateInterval> candidates;
    double minimum_objective = std::numeric_limits<double>::quiet_NaN();
    double minimum_lambda = std::numeric_limits<double>::quiet_NaN();
    bool complete = false;
  };

  struct Result {
    std::vector<StageRecord> stages;
    std::vector<CombinedPoint> points;
    std::vector<MemberAttempt> attempts;
    std::vector<profile_likelihood::MinimizationResult> final_members;
    double lambda_hat = std::numeric_limits<double>::quiet_NaN();
    double objective_hat = std::numeric_limits<double>::quiet_NaN();
    bool scan_coverage_complete = false;
    bool search_converged = false;
    bool unresolved_gap = false;
    bool at_boundary = false;
    std::string failure_reason;
  };

  using MinimizeMemberFunction = std::function<profile_likelihood::MinimizationResult(
      std::size_t member_index, double lambda, const std::vector<double>& seed)>;
  using StageCommitFunction = std::function<void(const Result&)>;

  /**
   * Profile every member at common lambda coordinates and minimize their unweighted sum.
   * Seeds are selected only from the deterministic nominal/coordinate-neighbor order;
   * invalid member points remain in the returned numerical record.
   */
  Result Run(const CombinedProfileConfig& config,
             double lambda_lower,
             double lambda_upper,
             const std::vector<std::vector<double>>& nominal_member_values,
             const MinimizeMemberFunction& minimize_member,
             const Result* resume = nullptr,
             const StageCommitFunction& commit_stage = {});

  /** Binary checkpoint format for complete stages; malformed or truncated input throws. */
  void SaveCheckpoint(std::ostream& output, const Result& result);
  Result LoadCheckpoint(std::istream& input);

}  // namespace exp_femto_3d::combined_profile
