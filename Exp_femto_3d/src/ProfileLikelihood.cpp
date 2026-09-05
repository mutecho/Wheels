#include "ProfileLikelihood.h"

#include <algorithm>
#include <cmath>

namespace exp_femto_3d::profile_likelihood {

  namespace {
    std::vector<double> MakeGrid(const double lower, const double upper, const int points) {
      std::vector<double> grid(static_cast<std::size_t>(points));
      for (int index = 0; index < points; ++index) {
        grid[static_cast<std::size_t>(index)] =
            lower + (upper - lower) * static_cast<double>(index) / static_cast<double>(points - 1);
      }
      return grid;
    }

    bool Better(const MinimizationResult& candidate, const MinimizationResult& current) {
      return std::isfinite(candidate.objective)
             && (!std::isfinite(current.objective) || candidate.objective < current.objective);
    }

    bool IsValidSliceEvaluation(const SliceEvaluation& evaluation) {
      return evaluation.model_domain_valid && evaluation.objective_valid && std::isfinite(evaluation.objective);
    }

    struct GridSpec {
      std::vector<std::vector<double>> axes;
      int stage = 0;
    };

    void RunGrid(const GridSpec& grid,
                 const std::vector<int>& target_indices,
                 const std::vector<double>& nominal_values,
                 const ProfileRetryStrategy retry_strategy,
                 const MinimizeFunction& minimize,
                 const SliceFunction& evaluate_slice,
                 ScanResult& output) {
      const int nx = static_cast<int>(grid.axes[0].size());
      const int ny = grid.axes.size() == 2U ? static_cast<int>(grid.axes[1].size()) : 1;
      std::vector<int> forward_indices(static_cast<std::size_t>(nx * ny), -1);
      const auto flat = [nx](const int ix, const int iy) { return iy * nx + ix; };

      const auto finalize = [&](PointRecord& point) {
        MinimizationResult best_valid;
        MinimizationResult best_failed;
        bool have_failed_attempt = false;
        point.valid_attempt_count = 0;
        double min_objective = std::numeric_limits<double>::infinity();
        double max_objective = -std::numeric_limits<double>::infinity();
        for (const AttemptRecord& attempt : output.attempts) {
          if (attempt.stage != point.stage || attempt.ix != point.ix || attempt.iy != point.iy) continue;
          if (IsValid(attempt.result)) {
            ++point.valid_attempt_count;
            if (!IsValid(best_valid) || Better(attempt.result, best_valid)) {
              best_valid = attempt.result;
              point.winner_seed = attempt.seed_origin;
            }
          } else if (!have_failed_attempt || Better(attempt.result, best_failed)) {
            // Preserve the best numerical failure, including finite penalty-valued
            // domain/objective failures, so the point keeps its diagnostic cause.
            best_failed = attempt.result;
            point.winner_seed = attempt.seed_origin;
            have_failed_attempt = true;
          }
          if (std::isfinite(attempt.result.objective)) {
            min_objective = std::min(min_objective, attempt.result.objective);
            max_objective = std::max(max_objective, attempt.result.objective);
          }
        }
        point.winner = IsValid(best_valid) ? best_valid : best_failed;
        point.status = IsValid(best_valid) ? PointStatus::kValid : StatusFor(point.winner);
        if (!IsValid(best_valid) && (!have_failed_attempt || !std::isfinite(point.winner.objective))) {
          point.status = PointStatus::kNoValidAttempt;
        }
        if (std::isfinite(min_objective) && std::isfinite(max_objective)) point.objective_spread = max_objective - min_objective;
      };
      const auto solve_sweep = [&](const bool reverse) {
        for (int outer = 0; outer < ny; ++outer) {
          const int iy = reverse ? ny - 1 - outer : outer;
          for (int inner = 0; inner < nx; ++inner) {
            const int ix = reverse ? nx - 1 - inner : inner;
            std::vector<double> coordinates{grid.axes[0][static_cast<std::size_t>(ix)]};
            if (grid.axes.size() == 2U) {
              coordinates.push_back(grid.axes[1][static_cast<std::size_t>(iy)]);
            }
            const int flat_index = flat(ix, iy);
            if (!reverse) {
              PointRecord point;
              point.point_index = static_cast<int>(output.points.size());
              point.stage = grid.stage;
              point.ix = ix;
              point.iy = grid.axes.size() == 2U ? iy : -1;
              point.coordinates = coordinates;
              const SliceEvaluation slice_evaluation = evaluate_slice(target_indices, coordinates);
              point.likelihood_slice_objective = slice_evaluation.objective;
              point.likelihood_slice_objective_valid = IsValidSliceEvaluation(slice_evaluation);
              output.points.push_back(std::move(point));
              forward_indices[static_cast<std::size_t>(flat_index)] = static_cast<int>(output.points.size() - 1U);
            }
            PointRecord& point = output.points[static_cast<std::size_t>(forward_indices[static_cast<std::size_t>(flat_index)])];

            const auto record_attempt = [&](const std::vector<double>& seed, const SeedOrigin origin) {
              AttemptRecord attempt;
              attempt.attempt_index = static_cast<int>(output.attempts.size());
              attempt.point_index = point.point_index;
              attempt.stage = point.stage;
              attempt.ix = point.ix;
              attempt.iy = point.iy;
              attempt.coordinates = coordinates;
              attempt.seed_origin = origin;
              attempt.result = minimize(seed, target_indices, coordinates);
              attempt.status = StatusFor(attempt.result);
              ++point.attempt_count;
              if (IsValid(attempt.result)) {
                ++point.valid_attempt_count;
              }
              output.attempts.push_back(std::move(attempt));
            };
            if (!reverse) record_attempt(nominal_values, SeedOrigin::kNominal);

            if (retry_strategy == ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors) {
              // Only use already-solved adjacent axis points: no row-wrap neighbors in 2D.
              int best_neighbor = -1;
              const auto consider = [&](const int candidate) {
                if (candidate < 0) {
                  return;
                }
                const PointRecord& candidate_point = output.points[static_cast<std::size_t>(candidate)];
                if (candidate_point.status != PointStatus::kValid) {
                  return;
                }
                if (best_neighbor < 0 || Better(candidate_point.winner, output.points[static_cast<std::size_t>(best_neighbor)].winner)) {
                  best_neighbor = candidate;
                }
              };
              if ((!reverse && ix > 0) || (reverse && ix + 1 < nx)) {
                consider(forward_indices[static_cast<std::size_t>(flat(reverse ? ix + 1 : ix - 1, iy))]);
              }
              if (grid.axes.size() == 2U && ((!reverse && iy > 0) || (reverse && iy + 1 < ny))) {
                consider(forward_indices[static_cast<std::size_t>(flat(ix, reverse ? iy + 1 : iy - 1))]);
              }
              if (best_neighbor >= 0) {
                record_attempt(output.points[static_cast<std::size_t>(best_neighbor)].winner.values,
                               reverse ? SeedOrigin::kReverseNeighbor : SeedOrigin::kForwardNeighbor);
              }
            }

            finalize(point);
          }
        }
      };
      solve_sweep(false);
      if (retry_strategy == ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors) {
        solve_sweep(true);
      }
    }
  }  // namespace

  bool IsValid(const MinimizationResult& result) {
    return result.model_domain_valid && result.objective_valid && std::isfinite(result.objective)
           && result.migrad_status == 0 && !result.minimizer_error;
  }

  PointStatus StatusFor(const MinimizationResult& result) {
    if (IsValid(result)) return PointStatus::kValid;
    // A setup/runtime minimizer failure has no trustworthy final model evaluation.
    // Preserve that provenance instead of letting default invalid flags mask it.
    if (result.minimizer_error) return PointStatus::kMinimizerError;
    if (!result.model_domain_valid) return PointStatus::kModelDomainInvalid;
    if (!result.objective_valid || !std::isfinite(result.objective)) return PointStatus::kObjectiveInvalid;
    if (result.migrad_status != 0) return PointStatus::kNonconverged;
    return PointStatus::kNoValidAttempt;
  }

  const char* ToString(const PointStatus status) {
    switch (status) {
      case PointStatus::kValid: return "valid";
      case PointStatus::kNonconverged: return "nonconverged";
      case PointStatus::kMinimizerError: return "minimizer_error";
      case PointStatus::kModelDomainInvalid: return "model_domain_invalid";
      case PointStatus::kObjectiveInvalid: return "objective_invalid";
      case PointStatus::kNoValidAttempt: return "no_valid_attempt";
    }
    return "no_valid_attempt";
  }

  const char* ToString(const SeedOrigin origin) {
    switch (origin) {
      case SeedOrigin::kNominal: return "nominal";
      case SeedOrigin::kForwardNeighbor: return "forward_neighbor";
      case SeedOrigin::kReverseNeighbor: return "reverse_neighbor";
    }
    return "nominal";
  }

  ScanResult RunProfileScan(const ProfileScanConfig& config,
                            const std::vector<int>& target_indices,
                            const std::vector<double>& effective_min,
                            const std::vector<double>& effective_max,
                            const std::vector<double>& nominal_values,
                            const ProfileRetryStrategy retry_strategy,
                            const MinimizeFunction& minimize,
                            const SliceFunction& evaluate_slice) {
    ScanResult output;
    GridSpec coarse;
    for (std::size_t axis = 0; axis < target_indices.size(); ++axis) {
      coarse.axes.push_back(MakeGrid(effective_min[axis], effective_max[axis], config.points[axis]));
    }
    RunGrid(coarse, target_indices, nominal_values, retry_strategy, minimize, evaluate_slice, output);

    int best_point = -1;
    for (std::size_t index = 0; index < output.points.size(); ++index) {
      const PointRecord& point = output.points[index];
      if (point.stage == 0 && point.status == PointStatus::kValid
          && (best_point < 0 || Better(point.winner, output.points[static_cast<std::size_t>(best_point)].winner))) {
        best_point = static_cast<int>(index);
      }
    }
    if (!config.refine) {
      output.refinement_reason = "disabled";
      return output;
    }
    if (best_point < 0) {
      output.refinement_reason = "no_valid_coarse_minimum";
      return output;
    }
    const PointRecord& best = output.points[static_cast<std::size_t>(best_point)];
    const int ny = coarse.axes.size() == 2U ? static_cast<int>(coarse.axes[1].size()) : 1;
    if (best.ix == 0 || best.ix + 1 == static_cast<int>(coarse.axes[0].size())
        || (ny > 1 && (best.iy == 0 || best.iy + 1 == ny))) {
      output.refinement_reason = "coarse_minimum_on_boundary";
      return output;
    }
    GridSpec refined;
    refined.stage = 1;
    for (std::size_t axis = 0; axis < target_indices.size(); ++axis) {
      const int center = axis == 0 ? best.ix : best.iy;
      const int refined_points = config.refinement_points.empty() ? config.points[axis] : config.refinement_points[axis];
      refined.axes.push_back(MakeGrid(coarse.axes[axis][static_cast<std::size_t>(center - 1)],
                                      coarse.axes[axis][static_cast<std::size_t>(center + 1)],
                                      refined_points));
    }
    output.refinement_performed = true;
    output.refinement_reason = "interior_coarse_minimum";
    RunGrid(refined, target_indices, nominal_values, retry_strategy, minimize, evaluate_slice, output);
    return output;
  }

  ReferenceMinimum SelectGlobalReference(
      const double nominal_objective,
      const bool nominal_valid,
      const std::vector<std::pair<std::string, const ScanResult*>>& scans) {
    ReferenceMinimum reference;
    if (nominal_valid && std::isfinite(nominal_objective)) {
      reference.objective = nominal_objective;
      reference.source = "nominal";
    }
    for (const auto& [scan_id, scan] : scans) {
      if (scan == nullptr) continue;
      for (const PointRecord& point : scan->points) {
        if (point.status == PointStatus::kValid && std::isfinite(point.winner.objective)
            && (!std::isfinite(reference.objective) || point.winner.objective < reference.objective)) {
          reference.objective = point.winner.objective;
          reference.source = "profile:" + scan_id;
        }
      }
    }
    return reference;
  }

}  // namespace exp_femto_3d::profile_likelihood
