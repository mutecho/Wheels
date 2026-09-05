#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ProfileLikelihood.h"

namespace {

  using exp_femto_3d::ProfileRetryStrategy;
  using exp_femto_3d::ProfileScanConfig;
  using namespace exp_femto_3d::profile_likelihood;

  void Expect(const bool condition, const std::string& message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  MinimizationResult ValidResult(const std::vector<double>& values, const double objective) {
    MinimizationResult result;
    result.values = values;
    result.errors.assign(values.size(), 0.1);
    result.objective = objective;
    result.edm = 1.0e-8;
    result.migrad_status = 0;
    result.hesse_status = 0;
    result.minuit_istat = 3;
    result.model_domain_valid = true;
    result.objective_valid = true;
    return result;
  }

  const PointRecord& PointAt(const ScanResult& result, const int stage, const int ix, const int iy = -1) {
    const auto found = std::find_if(result.points.begin(), result.points.end(), [&](const PointRecord& point) {
      return point.stage == stage && point.ix == ix && point.iy == iy;
    });
    Expect(found != result.points.end(), "requested profile point is missing");
    return *found;
  }

  bool HasAttempt(const ScanResult& result, const int ix, const SeedOrigin origin) {
    return std::any_of(result.attempts.begin(), result.attempts.end(), [&](const AttemptRecord& attempt) {
      return attempt.stage == 0 && attempt.ix == ix && attempt.seed_origin == origin;
    });
  }

}  // namespace

int main() {
  ProfileScanConfig one_dimensional;
  one_dimensional.id = "branch_arbitration";
  one_dimensional.parameters = {"rout2"};
  one_dimensional.points = {3};

  // The scripted objective deliberately makes the reverse-neighbor seed best at x=1.
  // This verifies that all attempts survive and that the winner is chosen by objective,
  // rather than scan order or a fixed seed preference.
  const auto scripted_minimize = [](const std::vector<double>& seed,
                                    const std::vector<int>& fixed_indices,
                                    const std::vector<double>& fixed_values) {
    Expect(fixed_indices == std::vector<int>{2}, "1D target must be fixed at the requested Minuit index");
    Expect(fixed_values.size() == 1U, "1D scan must pass exactly one fixed coordinate");
    const double x = fixed_values[0];
    double objective = x == 0.0 ? 4.0 : x == 2.0 ? 2.0 : 10.0;
    if (x == 1.0 && !seed.empty() && std::abs(seed[0] - 0.0) < 1.0e-12) objective = 3.0;
    if (x == 1.0 && !seed.empty() && std::abs(seed[0] - 2.0) < 1.0e-12) objective = 1.0;
    return ValidResult({x, 0.25, -0.5}, objective);
  };
  const auto scripted_slice = [](const std::vector<int>& fixed_indices, const std::vector<double>& fixed_values) {
    Expect(fixed_indices == std::vector<int>{2}, "slice target mismatch");
    return 20.0 + fixed_values[0] * fixed_values[0];
  };
  const ScanResult arbitration = RunProfileScan(one_dimensional,
                                                {2},
                                                {0.0},
                                                {2.0},
                                                {1.0, 0.0, 0.0},
                                                ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors,
                                                scripted_minimize,
                                                scripted_slice);
  Expect(arbitration.points.size() == 3U, "1D coarse grid must contain the requested exact point count");
  Expect(arbitration.attempts.size() > arbitration.points.size(), "neighbor retries must retain separate attempts");
  Expect(HasAttempt(arbitration, 1, SeedOrigin::kNominal), "nominal attempt missing");
  Expect(HasAttempt(arbitration, 1, SeedOrigin::kForwardNeighbor), "forward-neighbor attempt missing");
  Expect(HasAttempt(arbitration, 1, SeedOrigin::kReverseNeighbor), "reverse-neighbor attempt missing");
  const PointRecord& arbitration_middle = PointAt(arbitration, 0, 1);
  Expect(arbitration_middle.status == PointStatus::kValid, "arbitrated point should remain valid");
  Expect(arbitration_middle.winner_seed == SeedOrigin::kReverseNeighbor,
         "lowest valid reverse-neighbor objective must win arbitration");
  Expect(std::abs(arbitration_middle.winner.objective - 1.0) < 1.0e-12,
         "arbitrated objective mismatch");
  Expect(arbitration_middle.likelihood_slice_objective > arbitration_middle.winner.objective,
         "a fitted profile must improve over the fixed-nuisance likelihood slice in this toy model");
  const ReferenceMinimum lower_profile_reference =
      SelectGlobalReference(5.0, true, {{one_dimensional.id, &arbitration}});
  Expect(std::abs(lower_profile_reference.objective - 1.0) < 1.0e-12
             && lower_profile_reference.source == "profile:branch_arbitration",
         "a lower valid profile minimum must replace the diagnostic reference");
  const ReferenceMinimum lower_nominal_reference =
      SelectGlobalReference(0.5, true, {{one_dimensional.id, &arbitration}});
  Expect(std::abs(lower_nominal_reference.objective - 0.5) < 1.0e-12
             && lower_nominal_reference.source == "nominal",
         "the nominal fit must remain the reference when every profile point is higher");

  ProfileScanConfig two_dimensional;
  two_dimensional.id = "two_dimensional";
  two_dimensional.parameters = {"rout2", "rside2"};
  two_dimensional.points = {3, 3};
  const std::vector<std::vector<double>> expected_coordinates = {
      {0.0, 5.0}, {1.0, 5.0}, {2.0, 5.0}, {0.0, 6.0}, {1.0, 6.0},
      {2.0, 6.0}, {0.0, 7.0}, {1.0, 7.0}, {2.0, 7.0}};
  const auto quadratic_minimize = [](const std::vector<double>&,
                                     const std::vector<int>& fixed_indices,
                                     const std::vector<double>& fixed_values) {
    Expect(fixed_indices == std::vector<int>({2, 3}), "2D targets must be fixed in axis order");
    Expect(fixed_values.size() == 2U, "2D scan must pass two exact fixed coordinates");
    const double x = fixed_values[0];
    const double y = fixed_values[1];
    return ValidResult({x, y, x + y}, (x - 1.0) * (x - 1.0) + (y - 6.0) * (y - 6.0));
  };
  const ScanResult two_d = RunProfileScan(two_dimensional,
                                          {2, 3},
                                          {0.0, 5.0},
                                          {2.0, 7.0},
                                          {1.0, 6.0, 0.0},
                                          ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors,
                                          quadratic_minimize,
                                          [](const std::vector<int>&, const std::vector<double>& values) {
                                            return 5.0 + values[0] + values[1];
                                          });
  Expect(two_d.points.size() == expected_coordinates.size(), "2D grid must contain every requested coordinate");
  for (const auto& expected : expected_coordinates) {
    const auto found = std::find_if(two_d.points.begin(), two_d.points.end(), [&](const PointRecord& point) {
      return point.stage == 0 && point.coordinates == expected;
    });
    Expect(found != two_d.points.end(), "2D grid omitted an exact requested coordinate");
  }

  ProfileScanConfig refined = one_dimensional;
  refined.id = "interior_refinement";
  refined.refine = true;
  refined.refinement_points = {5};
  const ScanResult refined_result = RunProfileScan(refined,
                                                    {2},
                                                    {0.0},
                                                    {2.0},
                                                    {1.0, 0.0, 0.0},
                                                    ProfileRetryStrategy::kReferenceOnly,
                                                    [](const std::vector<double>&,
                                                       const std::vector<int>&,
                                                       const std::vector<double>& values) {
                                                      return ValidResult(values, (values[0] - 1.0) * (values[0] - 1.0));
                                                    },
                                                    [](const std::vector<int>&, const std::vector<double>&) { return 9.0; });
  Expect(refined_result.refinement_performed, "interior coarse minimum must trigger one refinement layer");
  Expect(refined_result.refinement_reason == "interior_coarse_minimum", "unexpected refinement reason");
  Expect(refined_result.points.size() == 8U, "coarse plus one refined 1D grid must retain all points");
  Expect(PointAt(refined_result, 1, 0).coordinates[0] == 0.0
             && PointAt(refined_result, 1, 4).coordinates[0] == 2.0,
         "refined grid must use the neighboring coarse-coordinate interval exactly");

  const ScanResult boundary_result = RunProfileScan(refined,
                                                     {2},
                                                     {0.0},
                                                     {2.0},
                                                     {1.0, 0.0, 0.0},
                                                     ProfileRetryStrategy::kReferenceOnly,
                                                     [](const std::vector<double>&,
                                                        const std::vector<int>&,
                                                        const std::vector<double>& values) {
                                                       return ValidResult(values, values[0]);
                                                     },
                                                     [](const std::vector<int>&, const std::vector<double>&) { return 1.0; });
  Expect(!boundary_result.refinement_performed
             && boundary_result.refinement_reason == "coarse_minimum_on_boundary",
         "boundary coarse minimum must skip refinement with its explicit reason");

  const auto failure_result = [](const PointStatus expected_status, const MinimizationResult& template_result) {
    ProfileScanConfig config;
    config.id = "failure";
    config.parameters = {"rout2"};
    config.points = {3};
    const ScanResult result = RunProfileScan(config,
                                              {2},
                                              {0.0},
                                              {2.0},
                                              {0.0, 0.0, 0.0},
                                              ProfileRetryStrategy::kReferenceOnly,
                                              [template_result](const std::vector<double>&,
                                                                const std::vector<int>&,
                                                                const std::vector<double>& values) {
                                                MinimizationResult result = template_result;
                                                result.values = values;
                                                return result;
                                              },
                                              [](const std::vector<int>&, const std::vector<double>&) { return 0.0; });
    Expect(result.attempts.size() == 3U, "failed attempts must still be stored for every coarse coordinate");
    for (const PointRecord& point : result.points) Expect(point.status == expected_status, "failure classification mismatch");
    return result;
  };

  MinimizationResult nonconverged = ValidResult({}, 4.0);
  nonconverged.migrad_status = 4;
  Expect(StatusFor(nonconverged) == PointStatus::kNonconverged, "MIGRAD failure classification mismatch");
  (void)failure_result(PointStatus::kNonconverged, nonconverged);
  MinimizationResult minimizer_error = ValidResult({}, 4.0);
  minimizer_error.minimizer_error = true;
  Expect(StatusFor(minimizer_error) == PointStatus::kMinimizerError, "minimizer error classification mismatch");
  (void)failure_result(PointStatus::kMinimizerError, minimizer_error);
  MinimizationResult model_invalid = ValidResult({}, 4.0);
  model_invalid.model_domain_valid = false;
  Expect(StatusFor(model_invalid) == PointStatus::kModelDomainInvalid, "model-domain classification mismatch");
  (void)failure_result(PointStatus::kModelDomainInvalid, model_invalid);
  MinimizationResult objective_invalid = ValidResult({}, std::numeric_limits<double>::quiet_NaN());
  objective_invalid.objective_valid = false;
  Expect(StatusFor(objective_invalid) == PointStatus::kObjectiveInvalid, "objective-invalid classification mismatch");
  const ScanResult no_valid = failure_result(PointStatus::kNoValidAttempt, objective_invalid);
  Expect(no_valid.refinement_reason == "disabled", "non-refined invalid scans must report disabled refinement");

  ProfileScanConfig invalid_refinement = refined;
  invalid_refinement.id = "no_valid_refinement";
  const ScanResult no_valid_refinement = RunProfileScan(invalid_refinement,
                                                         {2},
                                                         {0.0},
                                                         {2.0},
                                                         {0.0, 0.0, 0.0},
                                                         ProfileRetryStrategy::kReferenceOnly,
                                                         [objective_invalid](const std::vector<double>&,
                                                                             const std::vector<int>&,
                                                                             const std::vector<double>& values) {
                                                           MinimizationResult result = objective_invalid;
                                                           result.values = values;
                                                           return result;
                                                         },
                                                         [](const std::vector<int>&, const std::vector<double>&) { return 0.0; });
  Expect(!no_valid_refinement.refinement_performed
             && no_valid_refinement.refinement_reason == "no_valid_coarse_minimum",
         "an all-invalid coarse grid must skip refinement explicitly");

  std::cout << "profile_likelihood_driver_test passed\n";
  return 0;
}
