#include "exp_femto_1d/CatsModel.h"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

#include "CATS.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "TF1.h"

namespace exp_femto_1d {

  namespace {

    constexpr double kChargedPionMassMeV = 139.57;

    // Keep the Gaussian source local so the model depends only on installed CATS headers.
    double GaussianSource(double *parameters) {
      const double &radius = parameters[1];
      const double &size = parameters[3];
      return 4.0 * kPi * radius * radius * std::pow(4.0 * kPi * size * size, -1.5)
             * std::exp(-(radius * radius) / (4.0 * size * size));
    }

    // Evaluate the fixed baseline polynomial in GeV/c as required by the public contract.
    double BaselinePolynomial(const double kstar_gev, const double *parameters) {
      return parameters[0]
             * (1.0 + parameters[1] * kstar_gev + parameters[2] * kstar_gev * kstar_gev
                + parameters[3] * std::pow(kstar_gev, 3) + parameters[4] * std::pow(kstar_gev, 4));
    }

  }  // namespace

  class PiPiCatsModel::Impl {
   public:
    Impl(const unsigned num_mom_bins, const double k_min_mev, const double k_max_mev, const bool use_coulomb)
        : source_parameters(CATSparameters::tSource, 1, true), use_coulomb_flag(use_coulomb) {
      source_parameters.SetParameter(0, 6.0);

      kitty.SetMomBins(num_mom_bins, k_min_mev, k_max_mev);
      kitty.SetAnaSource(GaussianSource, source_parameters);
      kitty.SetAnaSource(0, source_parameters.GetParameter(0));
      kitty.SetUseAnalyticSource(true);
      kitty.SetNumChannels(1);
      kitty.SetNumPW(0, 0);
      kitty.SetSpin(0, 0);
      kitty.SetChannelWeight(0, 1);
      kitty.SetPdgId(211, 211);
      kitty.SetQ1Q2(use_coulomb_flag ? 1 : 0);
      kitty.SetRedMass(0.5 * kChargedPionMassMeV);
      kitty.SetNotifications(CATS::nError);
      kitty.SetMaxNumThreads(1);
      kitty.KillTheCat();
    }

    // Recompute the correlation on every source update so TF1 callbacks stay self-contained.
    double Evaluate(const double kstar_mev, const double source_size_fm) {
      if (!std::isfinite(kstar_mev) || !std::isfinite(source_size_fm) || source_size_fm <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
      }

      kitty.SetAnaSource(0, source_size_fm, false);
      kitty.KillTheCat();
      return kitty.EvalCorrFun(kstar_mev);
    }

    CATS kitty;
    CATSparameters source_parameters;
    bool use_coulomb_flag = true;
  };

  PiPiCatsModel::PiPiCatsModel(const unsigned num_mom_bins,
                               const double k_min_mev,
                               const double k_max_mev,
                               const bool use_coulomb)
      : impl_(new Impl(num_mom_bins, k_min_mev, k_max_mev, use_coulomb)) {
  }

  PiPiCatsModel::~PiPiCatsModel() {
    delete impl_;
    impl_ = nullptr;
  }

  double PiPiCatsModel::Evaluate(const double kstar_mev, const double source_size_fm) {
    return impl_->Evaluate(kstar_mev, source_size_fm);
  }

  TF1 *PiPiCatsModel::BuildFitFunction(const std::string &name, const FitConfig &config) {
    auto *fit_function = new TF1(name.c_str(),
                                 [this](double *x, double *parameters) {
                                   const double kstar_gev = x[0];
                                   const double kstar_mev = 1000.0 * kstar_gev;
                                   const double source_size = parameters[5];
                                   return BaselinePolynomial(kstar_gev, parameters) * Evaluate(kstar_mev, source_size);
                                 },
                                 0.0,
                                 config.fit_kstar_max,
                                 6);
    fit_function->SetParName(0, "p0");
    fit_function->SetParName(1, "p1");
    fit_function->SetParName(2, "p2");
    fit_function->SetParName(3, "p3");
    fit_function->SetParName(4, "p4");
    fit_function->SetParName(5, "source_size");
    fit_function->SetParameters(config.baseline_p0_init,
                                config.baseline_p1_init,
                                config.baseline_p2_init,
                                config.baseline_p3_value,
                                config.baseline_p4_value,
                                config.source_size_init);
    fit_function->SetParLimits(0, config.baseline_p0_min, config.baseline_p0_max);
    fit_function->SetParLimits(1, config.baseline_p1_min, config.baseline_p1_max);
    fit_function->SetParLimits(2, config.baseline_p2_min, config.baseline_p2_max);
    fit_function->SetParLimits(5, config.source_size_min, config.source_size_max);
    if (config.baseline_p3_fixed) {
      fit_function->FixParameter(3, config.baseline_p3_value);
    } else {
      fit_function->SetParameter(3, config.baseline_p3_value);
    }
    if (config.baseline_p4_fixed) {
      fit_function->FixParameter(4, config.baseline_p4_value);
    } else {
      fit_function->SetParameter(4, config.baseline_p4_value);
    }
    return fit_function;
  }

}  // namespace exp_femto_1d
