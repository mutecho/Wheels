#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TMath.h"
#include "exp_femto_3d/Types.h"

#ifdef EXP_FEMTO_3D_HAS_CATS
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#endif

namespace {

  constexpr double kHbarC = 0.1973269804;
  constexpr double kPiPiLikeSignBohrRadiusFm = 387.5;
  constexpr double kChargedPionMassMeV = 139.57018;

  void Expect(bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  double ComputeKStarMeV(const double q_out, const double q_side, const double q_long) {
    return 0.5 * std::sqrt(q_out * q_out + q_side * q_side + q_long * q_long) * 1000.0;
  }

  double ComputeLikeSignPiPiGamowFactor(const double q_out, const double q_side, const double q_long) {
    const double q_magnitude = std::sqrt(q_out * q_out + q_side * q_side + q_long * q_long);
    const double k_star_fm = 0.5 * q_magnitude / kHbarC;
    if (k_star_fm <= 1.0e-12) {
      return 0.0;
    }
    const double eta = 1.0 / (k_star_fm * kPiPiLikeSignBohrRadiusFm);
    const double two_pi_eta = 2.0 * TMath::Pi() * eta;
    if (two_pi_eta > 700.0) {
      return 0.0;
    }
    const double denominator = std::exp(two_pi_eta) - 1.0;
    if (denominator <= 0.0) {
      return 0.0;
    }
    return std::max(0.0, std::min(two_pi_eta / denominator, 1.0));
  }

  double EvaluateKernelTable(const std::vector<double> &kstar_mev,
                             const std::vector<double> &values,
                             const double kstar) {
    if (kstar <= kstar_mev.front()) {
      return values.front();
    }
    if (kstar >= kstar_mev.back()) {
      return 1.0;
    }
    const auto upper = std::upper_bound(kstar_mev.begin(), kstar_mev.end(), kstar);
    const std::size_t hi = static_cast<std::size_t>(std::distance(kstar_mev.begin(), upper));
    const std::size_t lo = hi > 0 ? hi - 1 : 0;
    const double fraction = (kstar - kstar_mev[lo]) / (kstar_mev[hi] - kstar_mev[lo]);
    return values[lo] + fraction * (values[hi] - values[lo]);
  }

#ifdef EXP_FEMTO_3D_HAS_CATS
  double CatsSphericalGaussianSource(double *parameters) {
    const double radius = parameters[1];
    const double source_size = parameters[3];
    return 4.0 * TMath::Pi() * radius * radius
           * std::pow(4.0 * TMath::Pi() * source_size * source_size, -1.5)
           * std::exp(-(radius * radius) / (4.0 * source_size * source_size));
  }

  double EvaluateCatsPiPiKernel(const double radius_fm,
                                const int q1q2,
                                const short quantum_statistics,
                                const double kstar_mev) {
    CATS kitty;
    kitty.SetNotifications(CATS::nError);
    kitty.SetMomBins(250, 0.0, 250.0);
    CATSparameters source_parameters(CATSparameters::tSource, 1, true);
    source_parameters.SetParameter(0, radius_fm);
    kitty.SetAnaSource(CatsSphericalGaussianSource, source_parameters);
    kitty.SetUseAnalyticSource(true);
    kitty.SetNumChannels(1);
    kitty.SetNumPW(0, 0);
    kitty.SetSpin(0, 0);
    kitty.SetChannelWeight(0, 1.0);
    kitty.SetQ1Q2(q1q2);
    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(quantum_statistics);
    kitty.SetRedMass(0.5 * kChargedPionMassMeV);
    kitty.KillTheCat();
    return kitty.EvalCorrFun(kstar_mev);
  }
#endif

}  // namespace

int main() {
  using namespace exp_femto_3d;

  Expect(UsesCoulomb(CoulombMode::kNone) == false, "none should not use Coulomb");
  Expect(UsesCoulomb(CoulombMode::kGamow), "gamow should use Coulomb");
  Expect(UsesCoulomb(CoulombMode::kFiniteSource), "finite_source should use Coulomb");

  const double expected_kstar = 0.5 * std::sqrt(0.03 * 0.03 + 0.04 * 0.04 + 0.12 * 0.12) * 1000.0;
  Expect(std::abs(ComputeKStarMeV(0.03, 0.04, 0.12) - expected_kstar) < 1.0e-12,
         "q_inv GeV/c to k* MeV/c conversion mismatch");

  const std::vector<double> kstar_table = {1.0, 3.0, 5.0};
  const std::vector<double> value_table = {0.4, 0.7, 0.9};
  Expect(std::abs(EvaluateKernelTable(kstar_table, value_table, 0.2) - 0.4) < 1.0e-12,
         "low-k* table behavior should clamp to first finite bin");
  Expect(std::abs(EvaluateKernelTable(kstar_table, value_table, 4.0) - 0.8) < 1.0e-12,
         "in-table kernel interpolation mismatch");
  Expect(std::abs(EvaluateKernelTable(kstar_table, value_table, 6.0) - 1.0) < 1.0e-12,
         "high-k* table behavior should return unity outside support");

#ifdef EXP_FEMTO_3D_HAS_CATS
  const double no_fsi_no_qs = EvaluateCatsPiPiKernel(4.0, 0, 0, 20.0);
  Expect(std::abs(no_fsi_no_qs - 1.0) < 5.0e-3,
         "CATS no-FSI/no-QS reference should be unity for a spherical Gaussian source");

  const double kstar_mev = 20.0;
  const double q_inv_gev = 2.0 * kstar_mev / 1000.0;
  const double gamow = ComputeLikeSignPiPiGamowFactor(q_inv_gev, 0.0, 0.0);
  const double small_radius = EvaluateCatsPiPiKernel(0.05, 1, 0, kstar_mev);
  const double large_radius = EvaluateCatsPiPiKernel(8.0, 1, 0, kstar_mev);
  Expect(std::isfinite(small_radius) && std::isfinite(large_radius), "CATS finite-source values should be finite");
  Expect(std::abs(small_radius - gamow) < std::abs(large_radius - gamow),
         "small finite-source radius should be closer to point-source Gamow than a large source radius");
#endif

  std::cout << "coulomb_kernel_validation_test passed\n";
  return 0;
}
