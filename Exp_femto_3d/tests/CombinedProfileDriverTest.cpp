#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

#include "CombinedProfile.h"

namespace {
  void Expect(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
  }

  exp_femto_3d::profile_likelihood::MinimizationResult ValidResult(
      const double lambda, const double nuisance, const double objective) {
    exp_femto_3d::profile_likelihood::MinimizationResult result;
    result.values = {lambda, nuisance};
    result.objective = objective;
    result.edm = 0.0;
    result.migrad_status = 0;
    result.model_domain_valid = true;
    result.objective_valid = true;
    return result;
  }
}

int main() {
  using namespace exp_femto_3d;
  CombinedProfileConfig config;
  config.coarse_points = 21;
  config.refinement_points = 21;
  config.max_refinement_rounds = 6;
  config.lambda_tolerance = 1.0e-5;
  config.objective_tolerance = 1.0e-8;
  const std::vector<std::vector<double>> nominal{{0.5, 0.0}, {0.5, 0.0}};
  const auto quadratic = [](const std::size_t member, const double lambda,
                            const std::vector<double>&) {
    const double target = member == 0U ? 0.3 : 0.7;
    const double weight = member == 0U ? 1.0 : 3.0;
    const double nuisance = (member == 0U ? 1.0 : -0.5) + 2.0 * lambda;
    const double reference_constant = member == 0U ? 17.0 : 103.0;
    return ValidResult(lambda, nuisance,
                       reference_constant + weight * (lambda - target) * (lambda - target));
  };
  const auto result = combined_profile::Run(config, 0.0, 1.0, nominal, quadratic);
  Expect(std::abs(result.lambda_hat - 0.6) <= config.lambda_tolerance,
         "combined profile must minimize the unweighted sum of member objectives");
  Expect(result.final_members.size() == 2U
             && std::abs(result.final_members[0].values[1] - 2.2) < 1.0e-5,
         "each nuisance must be reoptimized at the selected shared lambda");
  Expect(result.search_converged, "regular quadratic profile should satisfy both stop tolerances");
  const double minimum_delta = std::min_element(
      result.points.begin(), result.points.end(), [](const auto& left, const auto& right) {
        const double lhs = left.valid ? left.delta_objective : std::numeric_limits<double>::infinity();
        const double rhs = right.valid ? right.delta_objective : std::numeric_limits<double>::infinity();
        return lhs < rhs;
      })->delta_objective;
  Expect(std::abs(minimum_delta) < 1.0e-7,
         "combined delta must use the combined final conditional reference, not member minima");

  const auto asymmetric = combined_profile::Run(
      config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        const double offset = lambda - 0.37;
        return ValidResult(lambda, static_cast<double>(member),
                           (1.0 + lambda) * offset * offset + 11.0 * member);
      });
  Expect(std::abs(asymmetric.lambda_hat - 0.37) <= config.lambda_tolerance,
         "asymmetric profiles must be refined in objective coordinates without parabolic assumptions");

  const auto multimodal = combined_profile::Run(
      config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        const double product = (lambda - 0.25) * (lambda - 0.75);
        return ValidResult(lambda, static_cast<double>(member), product * product + 5.0 * member);
      });
  Expect(multimodal.stages.front().candidates.size() == 2U,
         "all separated coarse-grid minimum candidates must be retained");
  Expect(std::abs(multimodal.lambda_hat - 0.25) <= config.lambda_tolerance
             || std::abs(multimodal.lambda_hat - 0.75) <= config.lambda_tolerance,
         "multimodal refinement must finish on one of the equal valid branches");

  auto flat_config = config;
  flat_config.max_refinement_rounds = 2;
  const auto flat = combined_profile::Run(
      flat_config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        return ValidResult(lambda, static_cast<double>(member), 7.0 + member);
      });
  Expect(!flat.search_converged && !flat.stages.back().candidates.empty()
             && flat.failure_reason.find("budget exhausted") != std::string::npos,
         "a flat unresolved minimum must retain candidates and report refinement-budget exhaustion");

  auto no_refinement_config = config;
  no_refinement_config.coarse_points = 7;
  no_refinement_config.max_refinement_rounds = 0;
  const auto budget_exhausted = combined_profile::Run(
      no_refinement_config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        return ValidResult(lambda, static_cast<double>(member), (lambda - 0.43) * (lambda - 0.43));
      });
  Expect(!budget_exhausted.search_converged
             && budget_exhausted.failure_reason.find("budget exhausted") != std::string::npos,
         "a valid center must remain distinct from outer-search convergence when no refinement is allowed");

  const auto boundary = combined_profile::Run(
      config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        return ValidResult(lambda, static_cast<double>(member), (lambda + 0.2) * (lambda + 0.2));
      });
  Expect(boundary.at_boundary && std::abs(boundary.lambda_hat) <= config.lambda_tolerance,
         "one-sided refinement must retain a boundary minimum");

  auto gap_config = config;
  gap_config.max_refinement_rounds = 2;
  const auto gap = combined_profile::Run(
      gap_config, 0.0, 1.0, nominal,
      [](const std::size_t member, const double lambda, const std::vector<double>&) {
        auto result = ValidResult(lambda, 0.0, (lambda - 0.5) * (lambda - 0.5));
        if (member == 1U && lambda > 0.45 && lambda < 0.55) {
          result.migrad_status = 4;
        }
        return result;
      });
  Expect(gap.unresolved_gap, "invalid member coordinates must remain visible as an unresolved gap");
  Expect(std::any_of(gap.points.begin(), gap.points.end(), [](const auto& point) {
           return !point.valid && !point.failure_reason.empty();
         }), "failed combined points must retain their member-level cause");

  std::stringstream checkpoint(std::ios::in | std::ios::out | std::ios::binary);
  combined_profile::SaveCheckpoint(checkpoint, gap);
  checkpoint.seekg(0);
  const combined_profile::Result restored = combined_profile::LoadCheckpoint(checkpoint);
  Expect(restored.stages.size() == gap.stages.size() && restored.points.size() == gap.points.size()
             && restored.attempts.size() == gap.attempts.size(),
         "stage checkpoint must preserve stages, winners, and every attempt");
  bool rejected_truncated = false;
  try {
    std::stringstream truncated(std::string("bad"), std::ios::in | std::ios::binary);
    (void)combined_profile::LoadCheckpoint(truncated);
  } catch (const std::runtime_error&) {
    rejected_truncated = true;
  }
  Expect(rejected_truncated, "truncated checkpoint input must fail explicitly");
  std::stringstream encoded(std::ios::in | std::ios::out | std::ios::binary);
  combined_profile::SaveCheckpoint(encoded, gap);
  std::string corrupted = encoded.str();
  corrupted.front() ^= 0x1;
  bool rejected_corrupt_magic = false;
  try {
    std::stringstream corrupt_stream(corrupted, std::ios::in | std::ios::binary);
    (void)combined_profile::LoadCheckpoint(corrupt_stream);
  } catch (const std::runtime_error&) {
    rejected_corrupt_magic = true;
  }
  Expect(rejected_corrupt_magic, "checkpoint corruption must be rejected before any state is reused");
  int resumed_calls = 0;
  const auto resumed = combined_profile::Run(
      gap_config, 0.0, 1.0, nominal,
      [&](const std::size_t member, const double lambda, const std::vector<double>& seed) {
        ++resumed_calls;
        return quadratic(member, lambda, seed);
      }, &restored);
  Expect(resumed.stages.size() == restored.stages.size() && resumed_calls == 2,
         "resume after the last complete stage must run only final member conditionals");

  return 0;
}
