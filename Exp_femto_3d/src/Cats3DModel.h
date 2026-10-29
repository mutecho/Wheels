#pragma once

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "exp_femto_3d/Types.h"

class TH1;
class TH3;
class TH3D;
namespace exp_femto_3d {
  class Logger;

  /** A trial source matrix is outside the Gaussian domain; never a valid CF value. */
  class Cats3DDomainError : public std::runtime_error {
   public:
    using std::runtime_error::runtime_error;
  };

  /** Immutable external pair weight. Histogram errors are retained, not fitted. */
  class Cats3DLambda {
   public:
    explicit Cats3DLambda(const Cats3DLambdaInput &input);
    ~Cats3DLambda();
    double Evaluate(double k_mev) const;
    const Cats3DLambdaInput &Input() const {
      return input_;
    }
    const TH1 *Histogram() const {
      return histogram_.get();
    }

   private:
    Cats3DLambdaInput input_;
    std::unique_ptr<TH1> histogram_;
    double scale_ = 1.0;
  };

  /**
   * One bounded, uniformly binned q grid and one CATS3D wave kernel. Momentum
   * coordinates enter in GeV; only the adapter converts them to MeV (no extra /2).
   * The matrix parameters are (oo,ss,ll,os,ol,sl) in fm^2. Kernel storage is reused
   * across sources, groups and lambda scans; only an exactly identical source
   * matrix can reuse a previous source/WF fold.
   */
  class Cats3DModel {
   public:
    Cats3DModel(const TH3 &reference, double fit_q_max, const Cats3DConfig &config, const Logger *logger = nullptr);
    ~Cats3DModel();
    Cats3DModel(const Cats3DModel &) = delete;
    Cats3DModel &operator=(const Cats3DModel &) = delete;
    static bool Available();
    // Pure preflight: checks binning, memory budget and lambda coverage before any
    // formal fit output is recreated. It does not build a wave-function kernel.
    static void ValidateInputGrid(const TH3 &reference,
                                  double qmax,
                                  const Cats3DConfig &config,
                                  const Cats3DLambda &input);
    bool Matches(const TH3 &reference, double fit_q_max) const;
    void SetLambda(std::shared_ptr<const Cats3DLambda> input);
    // The constant scan is a final lambda; scaling acts after purity conversion.
    void SetScan(std::optional<double> constant, double scale = 1.0);
    void Prepare(const double *parameters, bool full);
    double Evaluate(double qo, double qs, double ql, const double *parameters, bool full, bool q2_baseline);
    double ScalarLambda() const;
    std::unique_ptr<TH3D> CoreHistogram() const;
    std::unique_ptr<TH3D> FinalHistogram(const double *parameters, bool full, bool q2_baseline);
    std::unique_ptr<TH3D> LambdaHistogram() const;
    const Cats3DLambda &LambdaInput() const;
    // Human-readable provenance and automatic filename stem; neither is a
    // mandatory integrity token in a portable native kernel file.
    const std::string &Contract() const;
    const std::string &CacheKey() const;
    const std::string &KernelCachePath() const;
    std::size_t KernelBytes() const;
    /** Two complete float kernels during native cache loading, plus ROOT/input reserve. */
    double EstimatedPeakGiB() const;
    bool WithCoulomb() const;
    std::size_t FoldCount() const;
    int PmlMaxCalls() const;

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}  // namespace exp_femto_3d
