#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "exp_femto_1d/CatsModel.h"

namespace {

  void Expect(const bool condition, const std::string &message) {
    if (!condition) {
      throw std::runtime_error(message);
    }
  }

  double BaselinePolynomial(const double kstar_gev, const double p0, const double p1, const double p2) {
    return p0 * (1.0 + p1 * kstar_gev + p2 * kstar_gev * kstar_gev);
  }

}  // namespace

int main() {
  using namespace exp_femto_1d;

  FitConfig fit_config;
  fit_config.fit_kstar_max = 0.20;
  fit_config.use_coulomb = false;
  fit_config.baseline_p0_init = 1.0;
  fit_config.baseline_p0_min = 0.9;
  fit_config.baseline_p0_max = 1.1;
  fit_config.baseline_p1_init = 0.0;
  fit_config.baseline_p1_min = -0.1;
  fit_config.baseline_p1_max = 0.1;
  fit_config.baseline_p2_init = 0.0;
  fit_config.baseline_p2_min = -0.5;
  fit_config.baseline_p2_max = 0.5;
  fit_config.baseline_p3_fixed = true;
  fit_config.baseline_p4_fixed = true;
  fit_config.source_size_init = 5.5;
  fit_config.source_size_min = 3.0;
  fit_config.source_size_max = 12.0;
  fit_config.cats_num_mom_bins = 300;
  fit_config.cats_kmin_mev = 0.0;
  fit_config.cats_kmax_mev = 250.0;

  PiPiCatsModel model(
      fit_config.cats_num_mom_bins, fit_config.cats_kmin_mev, fit_config.cats_kmax_mev, fit_config.use_coulomb);
  Expect(std::isfinite(model.Evaluate(50.0, 6.0)), "CATS model should evaluate to a finite number");

  constexpr double kTrueP0 = 1.02;
  constexpr double kTrueP1 = 0.015;
  constexpr double kTrueP2 = 0.04;
  constexpr double kTrueSourceSize = 6.2;

  TH1D data_histogram("toy_cf", "toy_cf; k* (GeV/c); C(k*)", 40, 0.0, fit_config.fit_kstar_max);
  for (int bin = 1; bin <= data_histogram.GetNbinsX(); ++bin) {
    const double kstar_gev = data_histogram.GetXaxis()->GetBinCenter(bin);
    const double femto = model.Evaluate(1000.0 * kstar_gev, kTrueSourceSize);
    const double baseline = BaselinePolynomial(kstar_gev, kTrueP0, kTrueP1, kTrueP2);
    data_histogram.SetBinContent(bin, baseline * femto);
    data_histogram.SetBinError(bin, 0.01);
  }

  std::unique_ptr<TF1> fit_function(model.BuildFitFunction("toy_fit", fit_config));
  const TFitResultPtr fit_result_ptr = data_histogram.Fit(fit_function.get(), "QRS0");
  const int fit_status = static_cast<int>(fit_result_ptr);
  const TFitResult *fit_result = fit_result_ptr.Get();
  const int covariance_status = fit_result != nullptr ? fit_result->CovMatrixStatus() : -1;

  Expect(fit_status == 0, "toy CATS fit should succeed");
  Expect(covariance_status >= 0, "fit should report a covariance status");
  Expect(std::abs(fit_function->GetParameter(0) - kTrueP0) < 0.05, "p0 should be recovered within tolerance");
  Expect(std::abs(fit_function->GetParameter(5) - kTrueSourceSize) < 1.0,
         "source size should be recovered within tolerance");

  return 0;
}
