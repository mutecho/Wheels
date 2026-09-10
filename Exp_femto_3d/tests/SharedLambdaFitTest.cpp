#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "SharedLambdaFit.h"

namespace {

  void Expect(const bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
  }

  bool Near(const double left, const double right, const double tolerance = 1.0e-10) {
    return std::abs(left - right) <= tolerance * (1.0 + std::max(std::abs(left), std::abs(right)));
  }

}  // namespace

int main() {
  using namespace exp_femto_3d::shared_lambda;

  std::vector<std::vector<LocalParameter>> local(2);
  for (auto &member : local) {
    member = {{"norm", 1.0, 0.01, 0.5, 1.5, false},
              {"lambda", 0.5, 0.01, 0.0, 1.0, false},
              {"rout2", 25.0, 0.1, 0.01, 400.0, false},
              {"alpha", 2.0, 0.01, 0.0, 0.0, true}};
  }
  const ParameterLayout layout = BuildParameterLayout({"phi0", "phi1"}, local);
  Expect(layout.labels.size() == 5U, "layout must contain one lambda and two nuisances per slice");
  Expect(layout.labels.front() == "lambda" && layout.slice_ids.front() == "<shared>",
         "shared lambda must own the first covariance coordinate");
  Expect(layout.local_to_global[0][1] == 0 && layout.local_to_global[1][1] == 0,
         "every member must map lambda to the same global coordinate");
  Expect(layout.local_to_global[0][3] == -1 && layout.local_to_global[1][3] == -1,
         "fixed local parameters must stay outside the free covariance");
  const std::vector<double> free_values = {0.62, 1.1, 30.0, 0.9, 20.0};
  const std::vector<double> expanded = ExpandLocalValues(layout, 1U, local[1], free_values.data());
  Expect(Near(expanded[0], 0.9) && Near(expanded[1], 0.62) && Near(expanded[2], 20.0)
             && Near(expanded[3], 2.0),
         "local expansion must combine shared, per-slice, and fixed values");

  const std::vector<double> positive = {4.0, 1.0, 1.0, 2.0};
  const std::vector<double> singular = {1.0, 1.0, 1.0, 1.0};
  Expect(IsPositiveDefinite(positive, 2U), "known SPD covariance should pass");
  Expect(!IsPositiveDefinite(singular, 2U), "singular covariance must be rejected");

  const std::vector<double> phi = {0.0, 0.5, 1.0, 1.5};
  constexpr double intercept = 24.0;
  constexpr double harmonic = 1.75;
  std::vector<double> values;
  for (const double angle : phi) values.push_back(intercept + 2.0 * harmonic * std::cos(2.0 * angle));
  const std::vector<double> covariance = {
      0.25, 0.04, 0.01, 0.00,
      0.04, 0.36, 0.03, 0.01,
      0.01, 0.03, 0.49, 0.02,
      0.00, 0.01, 0.02, 0.64};
  const HarmonicResult fit = FitSecondHarmonicGLS(phi, values, covariance, false);
  Expect(fit.valid && Near(fit.intercept, intercept) && Near(fit.harmonic, harmonic),
         "correlated GLS must recover A and B in A+2B cos(2phi)");
  Expect(fit.intercept_variance > 0.0 && fit.harmonic_variance > 0.0,
         "GLS coefficient covariance must be positive");
  const HarmonicResult invalid = FitSecondHarmonicGLS({0.0, 1.0}, {1.0, 2.0}, singular, false);
  Expect(!invalid.valid && !invalid.failure_reason.empty(),
         "invalid covariance must not trigger an equal-weight fallback");

  return 0;
}
