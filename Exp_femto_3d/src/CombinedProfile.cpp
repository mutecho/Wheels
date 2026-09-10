#include "CombinedProfile.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <istream>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace exp_femto_3d::combined_profile {
  namespace {

    std::vector<double> Grid(const double lower, const double upper, const int count) {
      std::vector<double> coordinates(static_cast<std::size_t>(count));
      for (int index = 0; index < count; ++index) {
        coordinates[static_cast<std::size_t>(index)] =
            lower + (upper - lower) * static_cast<double>(index) / static_cast<double>(count - 1);
      }
      return coordinates;
    }

    bool Better(const profile_likelihood::MinimizationResult& candidate,
                const profile_likelihood::MinimizationResult& current) {
      return std::isfinite(candidate.objective)
             && (!std::isfinite(current.objective) || candidate.objective < current.objective);
    }

    std::vector<CandidateInterval> FindCandidates(const std::vector<CombinedPoint*>& ordered,
                                                   const double lower_bound,
                                                   const double upper_bound) {
      std::vector<CandidateInterval> intervals;
      for (std::size_t index = 0; index < ordered.size(); ++index) {
        const CombinedPoint& point = *ordered[index];
        if (!point.valid) continue;
        const CombinedPoint* left = index > 0U ? ordered[index - 1U] : nullptr;
        const CombinedPoint* right = index + 1U < ordered.size() ? ordered[index + 1U] : nullptr;
        const bool no_higher_than_left = left == nullptr || !left->valid || point.objective <= left->objective;
        const bool no_higher_than_right = right == nullptr || !right->valid || point.objective <= right->objective;
        if (!no_higher_than_left || !no_higher_than_right) continue;
        CandidateInterval interval;
        interval.lower = left == nullptr ? lower_bound : left->lambda;
        interval.upper = right == nullptr ? upper_bound : right->lambda;
        interval.touches_lower_bound = left == nullptr || point.lambda == lower_bound;
        interval.touches_upper_bound = right == nullptr || point.lambda == upper_bound;
        interval.contains_failure_gap = (left != nullptr && !left->valid) || (right != nullptr && !right->valid);
        intervals.push_back(interval);
      }
      std::sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        return left.lower < right.lower || (left.lower == right.lower && left.upper < right.upper);
      });
      std::vector<CandidateInterval> merged;
      for (const CandidateInterval& interval : intervals) {
        if (merged.empty() || interval.lower > merged.back().upper) {
          merged.push_back(interval);
        } else {
          merged.back().upper = std::max(merged.back().upper, interval.upper);
          merged.back().touches_lower_bound = merged.back().touches_lower_bound || interval.touches_lower_bound;
          merged.back().touches_upper_bound = merged.back().touches_upper_bound || interval.touches_upper_bound;
          merged.back().contains_failure_gap = merged.back().contains_failure_gap || interval.contains_failure_gap;
        }
      }
      return merged;
    }

  }  // namespace

  Result Run(const CombinedProfileConfig& config,
             const double lambda_lower,
             const double lambda_upper,
             const std::vector<std::vector<double>>& nominal_member_values,
             const MinimizeMemberFunction& minimize_member,
             const Result* resume,
             const StageCommitFunction& commit_stage) {
    if (!(lambda_lower < lambda_upper) || config.coarse_points < 3 || config.refinement_points < 3
        || config.max_refinement_rounds < 0 || nominal_member_values.empty() || !minimize_member) {
      throw std::invalid_argument("invalid combined-profile driver contract");
    }
    Result output = resume == nullptr ? Result{} : *resume;
    output.final_members.clear();
    output.lambda_hat = std::numeric_limits<double>::quiet_NaN();
    output.objective_hat = std::numeric_limits<double>::quiet_NaN();
    output.failure_reason.clear();
    std::map<double, std::size_t> point_by_coordinate;
    for (std::size_t index = 0; index < output.points.size(); ++index) {
      point_by_coordinate[output.points[index].lambda] = index;
    }
    std::vector<std::size_t> previous_forward_points;

    auto run_stage = [&](const int stage, std::vector<double> coordinates) {
      std::sort(coordinates.begin(), coordinates.end());
      coordinates.erase(std::unique(coordinates.begin(), coordinates.end()), coordinates.end());
      StageRecord record;
      record.stage = stage;
      record.coordinates = coordinates;
      std::vector<std::size_t> stage_points;
      previous_forward_points.clear();
      for (std::size_t coordinate_index = 0; coordinate_index < coordinates.size(); ++coordinate_index) {
        const double lambda = coordinates[coordinate_index];
        const auto existing = point_by_coordinate.find(lambda);
        if (existing != point_by_coordinate.end()) {
          stage_points.push_back(existing->second);
          previous_forward_points.assign(nominal_member_values.size(), static_cast<std::size_t>(-1));
          for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
            if (profile_likelihood::IsValid(output.points[existing->second].members[member].winner)) {
              previous_forward_points[member] = existing->second;
            }
          }
          continue;
        }
        CombinedPoint point;
        point.stage = stage;
        point.point_index = static_cast<int>(output.points.size());
        point.lambda = lambda;
        point.members.resize(nominal_member_values.size());
        for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
          MemberPoint& member_point = point.members[member];
          member_point.member_index = static_cast<int>(member);
          auto attempt = [&](const std::vector<double>& seed,
                             const profile_likelihood::SeedOrigin origin) {
            MemberAttempt saved;
            saved.stage = stage;
            saved.point_index = point.point_index;
            saved.member_index = static_cast<int>(member);
            saved.lambda = lambda;
            saved.seed_origin = origin;
            saved.result = minimize_member(member, lambda, seed);
            ++member_point.attempt_count;
            if (profile_likelihood::IsValid(saved.result)) {
              ++member_point.valid_attempt_count;
              if (!profile_likelihood::IsValid(member_point.winner)
                  || Better(saved.result, member_point.winner)) {
                member_point.winner = saved.result;
                member_point.winner_seed = origin;
              }
            } else if (member_point.valid_attempt_count == 0
                       && Better(saved.result, member_point.winner)) {
              member_point.winner = saved.result;
              member_point.winner_seed = origin;
            }
            output.attempts.push_back(std::move(saved));
          };
          attempt(nominal_member_values[member], profile_likelihood::SeedOrigin::kNominal);
          if (config.retry_strategy == ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors
              && previous_forward_points.size() == nominal_member_values.size()
              && previous_forward_points[member] != static_cast<std::size_t>(-1)) {
            attempt(output.points[previous_forward_points[member]].members[member].winner.values,
                    profile_likelihood::SeedOrigin::kForwardNeighbor);
          }
          member_point.status = profile_likelihood::StatusFor(member_point.winner);
        }
        output.points.push_back(std::move(point));
        const std::size_t saved_index = output.points.size() - 1U;
        point_by_coordinate[lambda] = saved_index;
        stage_points.push_back(saved_index);
        previous_forward_points.resize(nominal_member_values.size(), saved_index);
        for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
          if (!profile_likelihood::IsValid(output.points[saved_index].members[member].winner)) {
            previous_forward_points[member] = static_cast<std::size_t>(-1);
          }
        }
      }

      if (config.retry_strategy == ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors) {
        std::vector<std::size_t> next_valid(nominal_member_values.size(), static_cast<std::size_t>(-1));
        for (std::size_t reverse = stage_points.size(); reverse-- > 0U;) {
          CombinedPoint& point = output.points[stage_points[reverse]];
          if (point.stage != stage) {
            for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
              if (profile_likelihood::IsValid(point.members[member].winner)) {
                next_valid[member] = stage_points[reverse];
              }
            }
            continue;
          }
          for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
            if (next_valid[member] != static_cast<std::size_t>(-1)) {
              MemberAttempt saved;
              saved.stage = stage;
              saved.point_index = point.point_index;
              saved.member_index = static_cast<int>(member);
              saved.lambda = point.lambda;
              saved.seed_origin = profile_likelihood::SeedOrigin::kReverseNeighbor;
              saved.result = minimize_member(member, point.lambda,
                  output.points[next_valid[member]].members[member].winner.values);
              MemberPoint& member_point = point.members[member];
              ++member_point.attempt_count;
              if (profile_likelihood::IsValid(saved.result)) {
                ++member_point.valid_attempt_count;
                if (!profile_likelihood::IsValid(member_point.winner)
                    || Better(saved.result, member_point.winner)) {
                  member_point.winner = saved.result;
                  member_point.winner_seed = saved.seed_origin;
                }
              }
              output.attempts.push_back(std::move(saved));
            }
            point.members[member].status = profile_likelihood::StatusFor(point.members[member].winner);
            if (profile_likelihood::IsValid(point.members[member].winner)) next_valid[member] = stage_points[reverse];
          }
        }
      }

      for (const std::size_t index : stage_points) {
        CombinedPoint& point = output.points[index];
        point.valid = std::all_of(point.members.begin(), point.members.end(), [](const MemberPoint& member) {
          return profile_likelihood::IsValid(member.winner);
        });
        if (point.valid) {
          point.objective = 0.0;
          for (const MemberPoint& member : point.members) point.objective += member.winner.objective;
          if (!std::isfinite(point.objective)) point.valid = false;
        }
        if (!point.valid) {
          std::ostringstream reason;
          reason << "invalid members:";
          for (const MemberPoint& member : point.members) {
            if (!profile_likelihood::IsValid(member.winner)) {
              reason << ' ' << member.member_index << '=' << profile_likelihood::ToString(member.status);
            }
          }
          point.failure_reason = reason.str();
        } else if (!std::isfinite(record.minimum_objective) || point.objective < record.minimum_objective) {
          record.minimum_objective = point.objective;
          record.minimum_lambda = point.lambda;
        }
      }
      std::vector<CombinedPoint*> ordered;
      ordered.reserve(stage_points.size());
      for (const std::size_t index : stage_points) ordered.push_back(&output.points[index]);
      record.candidates = FindCandidates(ordered, lambda_lower, lambda_upper);
      record.complete = true;
      output.stages.push_back(std::move(record));
      if (commit_stage) commit_stage(output);
    };

    if (output.stages.empty()) {
      run_stage(0, Grid(lambda_lower, lambda_upper, config.coarse_points));
    }
    double previous_minimum = output.stages.back().minimum_objective;
    const int first_round = output.stages.back().stage + 1;
    for (int round = first_round; round <= config.max_refinement_rounds; ++round) {
      const auto& candidates = output.stages.back().candidates;
      if (candidates.empty()) break;
      double widest = 0.0;
      std::vector<double> coordinates;
      for (const CandidateInterval& candidate : candidates) {
        widest = std::max(widest, candidate.upper - candidate.lower);
        const auto grid = Grid(candidate.lower, candidate.upper, config.refinement_points);
        coordinates.insert(coordinates.end(), grid.begin(), grid.end());
      }
      if (widest <= config.lambda_tolerance) {
        output.search_converged = true;
        break;
      }
      run_stage(round, std::move(coordinates));
      const double current_minimum = output.stages.back().minimum_objective;
      double current_widest = 0.0;
      for (const CandidateInterval& candidate : output.stages.back().candidates) {
        current_widest = std::max(current_widest, candidate.upper - candidate.lower);
      }
      if (current_widest <= config.lambda_tolerance && std::isfinite(previous_minimum)
          && std::isfinite(current_minimum)
          && std::abs(current_minimum - previous_minimum) <= config.objective_tolerance) {
        output.search_converged = true;
        break;
      }
      previous_minimum = current_minimum;
    }

    const CombinedPoint* best = nullptr;
    for (const CombinedPoint& point : output.points) {
      if (point.valid && (best == nullptr || point.objective < best->objective)) best = &point;
    }
    output.scan_coverage_complete = !output.stages.empty() && output.stages.front().complete;
    output.unresolved_gap = std::any_of(output.points.begin(), output.points.end(), [](const CombinedPoint& point) {
      return !point.valid;
    });
    if (best == nullptr) {
      output.failure_reason = "all combined-profile coordinates are invalid";
      return output;
    }
    output.lambda_hat = best->lambda;
    output.objective_hat = best->objective;
    output.at_boundary = best->lambda == lambda_lower || best->lambda == lambda_upper;
    for (std::size_t member = 0; member < nominal_member_values.size(); ++member) {
      output.final_members.push_back(minimize_member(member, output.lambda_hat, best->members[member].winner.values));
      if (!profile_likelihood::IsValid(output.final_members.back())) {
        output.failure_reason = "final conditional minimization failed for member " + std::to_string(member);
      }
    }
    double final_objective = 0.0;
    for (const auto& member : output.final_members) final_objective += member.objective;
    if (output.failure_reason.empty() && std::isfinite(final_objective)) output.objective_hat = final_objective;
    for (CombinedPoint& point : output.points) {
      if (point.valid) point.delta_objective = point.objective - output.objective_hat;
    }
    if (!output.search_converged && output.failure_reason.empty()) {
      output.failure_reason = "combined-profile refinement budget exhausted before convergence";
    }
    return output;
  }

  namespace {
    constexpr std::uint64_t kCheckpointMagic = 0x4350524f46494c31ULL;

    template <typename T>
    void WriteScalar(std::ostream& output, const T& value) {
      output.write(reinterpret_cast<const char*>(&value), sizeof(T));
      if (!output) throw std::runtime_error("cannot write combined-profile checkpoint");
    }

    template <typename T>
    T ReadScalar(std::istream& input) {
      T value{};
      input.read(reinterpret_cast<char*>(&value), sizeof(T));
      if (!input) throw std::runtime_error("truncated combined-profile checkpoint");
      return value;
    }

    void WriteString(std::ostream& output, const std::string& value) {
      WriteScalar(output, static_cast<std::uint64_t>(value.size()));
      output.write(value.data(), static_cast<std::streamsize>(value.size()));
      if (!output) throw std::runtime_error("cannot write combined-profile checkpoint string");
    }

    std::string ReadString(std::istream& input) {
      const std::uint64_t size = ReadScalar<std::uint64_t>(input);
      if (size > (1ULL << 30U)) throw std::runtime_error("invalid combined-profile checkpoint string size");
      std::string value(static_cast<std::size_t>(size), '\0');
      input.read(value.data(), static_cast<std::streamsize>(size));
      if (!input) throw std::runtime_error("truncated combined-profile checkpoint string");
      return value;
    }

    template <typename T>
    void WriteVector(std::ostream& output, const std::vector<T>& values) {
      WriteScalar(output, static_cast<std::uint64_t>(values.size()));
      for (const T& value : values) WriteScalar(output, value);
    }

    template <typename T>
    std::vector<T> ReadVector(std::istream& input) {
      const std::uint64_t size = ReadScalar<std::uint64_t>(input);
      if (size > (1ULL << 30U)) throw std::runtime_error("invalid combined-profile checkpoint vector size");
      std::vector<T> values(static_cast<std::size_t>(size));
      for (T& value : values) value = ReadScalar<T>(input);
      return values;
    }

    void WriteMinimization(std::ostream& output, const profile_likelihood::MinimizationResult& value) {
      WriteVector(output, value.values);
      WriteVector(output, value.errors);
      WriteScalar(output, value.objective);
      WriteScalar(output, value.edm);
      WriteScalar(output, value.migrad_status);
      WriteScalar(output, value.hesse_status);
      WriteScalar(output, value.minuit_istat);
      WriteScalar(output, value.model_domain_valid);
      WriteScalar(output, value.objective_valid);
      WriteScalar(output, value.minimizer_error);
      WriteScalar(output, static_cast<std::uint64_t>(value.fcn_calls));
      WriteScalar(output, value.total_wall_ms);
      WriteScalar(output, value.migrad_wall_ms);
      WriteScalar(output, value.hesse_wall_ms);
      WriteScalar(output, value.hesse_ran);
      WriteScalar(output, value.parameter_errors_valid);
      WriteVector(output, value.at_lower_bound);
      WriteVector(output, value.at_upper_bound);
    }

    profile_likelihood::MinimizationResult ReadMinimization(std::istream& input) {
      profile_likelihood::MinimizationResult value;
      value.values = ReadVector<double>(input);
      value.errors = ReadVector<double>(input);
      value.objective = ReadScalar<double>(input);
      value.edm = ReadScalar<double>(input);
      value.migrad_status = ReadScalar<int>(input);
      value.hesse_status = ReadScalar<int>(input);
      value.minuit_istat = ReadScalar<int>(input);
      value.model_domain_valid = ReadScalar<bool>(input);
      value.objective_valid = ReadScalar<bool>(input);
      value.minimizer_error = ReadScalar<bool>(input);
      value.fcn_calls = static_cast<std::size_t>(ReadScalar<std::uint64_t>(input));
      value.total_wall_ms = ReadScalar<double>(input);
      value.migrad_wall_ms = ReadScalar<double>(input);
      value.hesse_wall_ms = ReadScalar<double>(input);
      value.hesse_ran = ReadScalar<bool>(input);
      value.parameter_errors_valid = ReadScalar<bool>(input);
      value.at_lower_bound = ReadVector<int>(input);
      value.at_upper_bound = ReadVector<int>(input);
      return value;
    }
  }

  void SaveCheckpoint(std::ostream& output, const Result& result) {
    WriteScalar(output, kCheckpointMagic);
    WriteScalar(output, static_cast<std::uint64_t>(result.stages.size()));
    for (const StageRecord& stage : result.stages) {
      WriteScalar(output, stage.stage);
      WriteVector(output, stage.coordinates);
      WriteScalar(output, static_cast<std::uint64_t>(stage.candidates.size()));
      for (const CandidateInterval& candidate : stage.candidates) {
        WriteScalar(output, candidate.lower);
        WriteScalar(output, candidate.upper);
        WriteScalar(output, candidate.touches_lower_bound);
        WriteScalar(output, candidate.touches_upper_bound);
        WriteScalar(output, candidate.contains_failure_gap);
      }
      WriteScalar(output, stage.minimum_objective);
      WriteScalar(output, stage.minimum_lambda);
      WriteScalar(output, stage.complete);
    }
    WriteScalar(output, static_cast<std::uint64_t>(result.points.size()));
    for (const CombinedPoint& point : result.points) {
      WriteScalar(output, point.stage);
      WriteScalar(output, point.point_index);
      WriteScalar(output, point.lambda);
      WriteScalar(output, static_cast<std::uint64_t>(point.members.size()));
      for (const MemberPoint& member : point.members) {
        WriteScalar(output, member.member_index);
        WriteMinimization(output, member.winner);
        WriteScalar(output, static_cast<int>(member.status));
        WriteScalar(output, static_cast<int>(member.winner_seed));
        WriteScalar(output, member.attempt_count);
        WriteScalar(output, member.valid_attempt_count);
      }
      WriteScalar(output, point.objective);
      WriteScalar(output, point.delta_objective);
      WriteScalar(output, point.valid);
      WriteString(output, point.failure_reason);
    }
    WriteScalar(output, static_cast<std::uint64_t>(result.attempts.size()));
    for (const MemberAttempt& attempt : result.attempts) {
      WriteScalar(output, attempt.stage);
      WriteScalar(output, attempt.point_index);
      WriteScalar(output, attempt.member_index);
      WriteScalar(output, attempt.lambda);
      WriteScalar(output, static_cast<int>(attempt.seed_origin));
      WriteMinimization(output, attempt.result);
    }
    WriteScalar(output, result.scan_coverage_complete);
    WriteScalar(output, result.search_converged);
    WriteScalar(output, result.unresolved_gap);
    WriteScalar(output, result.at_boundary);
  }

  Result LoadCheckpoint(std::istream& input) {
    if (ReadScalar<std::uint64_t>(input) != kCheckpointMagic) {
      throw std::runtime_error("unsupported combined-profile checkpoint format");
    }
    Result result;
    const std::uint64_t stage_count = ReadScalar<std::uint64_t>(input);
    for (std::uint64_t index = 0; index < stage_count; ++index) {
      StageRecord stage;
      stage.stage = ReadScalar<int>(input);
      stage.coordinates = ReadVector<double>(input);
      const std::uint64_t candidate_count = ReadScalar<std::uint64_t>(input);
      for (std::uint64_t candidate_index = 0; candidate_index < candidate_count; ++candidate_index) {
        CandidateInterval candidate;
        candidate.lower = ReadScalar<double>(input);
        candidate.upper = ReadScalar<double>(input);
        candidate.touches_lower_bound = ReadScalar<bool>(input);
        candidate.touches_upper_bound = ReadScalar<bool>(input);
        candidate.contains_failure_gap = ReadScalar<bool>(input);
        stage.candidates.push_back(candidate);
      }
      stage.minimum_objective = ReadScalar<double>(input);
      stage.minimum_lambda = ReadScalar<double>(input);
      stage.complete = ReadScalar<bool>(input);
      if (!stage.complete) throw std::runtime_error("combined-profile checkpoint contains an incomplete stage");
      result.stages.push_back(std::move(stage));
    }
    const std::uint64_t point_count = ReadScalar<std::uint64_t>(input);
    for (std::uint64_t index = 0; index < point_count; ++index) {
      CombinedPoint point;
      point.stage = ReadScalar<int>(input);
      point.point_index = ReadScalar<int>(input);
      point.lambda = ReadScalar<double>(input);
      const std::uint64_t member_count = ReadScalar<std::uint64_t>(input);
      for (std::uint64_t member_index = 0; member_index < member_count; ++member_index) {
        MemberPoint member;
        member.member_index = ReadScalar<int>(input);
        member.winner = ReadMinimization(input);
        member.status = static_cast<profile_likelihood::PointStatus>(ReadScalar<int>(input));
        member.winner_seed = static_cast<profile_likelihood::SeedOrigin>(ReadScalar<int>(input));
        member.attempt_count = ReadScalar<int>(input);
        member.valid_attempt_count = ReadScalar<int>(input);
        point.members.push_back(std::move(member));
      }
      point.objective = ReadScalar<double>(input);
      point.delta_objective = ReadScalar<double>(input);
      point.valid = ReadScalar<bool>(input);
      point.failure_reason = ReadString(input);
      result.points.push_back(std::move(point));
    }
    const std::uint64_t attempt_count = ReadScalar<std::uint64_t>(input);
    for (std::uint64_t index = 0; index < attempt_count; ++index) {
      MemberAttempt attempt;
      attempt.stage = ReadScalar<int>(input);
      attempt.point_index = ReadScalar<int>(input);
      attempt.member_index = ReadScalar<int>(input);
      attempt.lambda = ReadScalar<double>(input);
      attempt.seed_origin = static_cast<profile_likelihood::SeedOrigin>(ReadScalar<int>(input));
      attempt.result = ReadMinimization(input);
      result.attempts.push_back(std::move(attempt));
    }
    result.scan_coverage_complete = ReadScalar<bool>(input);
    result.search_converged = ReadScalar<bool>(input);
    result.unresolved_gap = ReadScalar<bool>(input);
    result.at_boundary = ReadScalar<bool>(input);
    return result;
  }

}  // namespace exp_femto_3d::combined_profile
