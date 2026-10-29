#pragma once

#include <array>
#include <memory>
#include <string>

#include "exp_femto_3d/Types.h"

class TAxis;
class TH3;
class TTree;
class TTreeReader;

namespace exp_femto_3d {

  constexpr double kCfNormalizationTolerance = 1.e-9;  // GeV/c
  std::string ToString(CfNormalizationMode mode);
  void ValidateCfNormalizationConfig(const CfNormalizationConfig &config);

  /** Resolve normal-bin inner/outer boxes using centers, inclusive shell endpoints, no clipping. */
  CfNormalizationResult ResolveCfNormalizationRegion(const CfNormalizationConfig &config,
                                                     const std::array<const TAxis *, 3> &axes);
  void ValidateCfNormalizationAxes(const std::array<const TAxis *, 3> &se, const std::array<const TAxis *, 3> &me);
  /** Preserve the legacy ROOT Integral("width") convention, including for nonuniform grids. */
  double CfNormalizationIntegral(const TH3 &hist, const CfNormalizationResult &region);
  CfNormalizationResult ComputeCfNormalization(const TH3 &se, const TH3 &me, const CfNormalizationConfig &config);
  void ValidateCfNormalizationResult(const CfNormalizationResult &result);
  /** Convert a normalized model to raw S/B. Only legacy files recompute full integrals. */
  double ResolvePMLRawScale(const TH3 &se, const TH3 &me, const CfNormalizationResult &stored);

  /** Shared fixed-schema ROOT branches for CF, fit, and diagnostic provenance. */
  class CfNormalizationTreeWriter {
   public:
    explicit CfNormalizationTreeWriter(TTree &tree);
    CfNormalizationTreeWriter(const CfNormalizationTreeWriter &) = delete;
    CfNormalizationTreeWriter &operator=(const CfNormalizationTreeWriter &) = delete;
    void Set(const CfNormalizationResult &result);

   private:
    CfNormalizationResult value_;
    std::string mode_;
    int version_ = 1;
  };

  class CfNormalizationTreeReader {
   public:
    CfNormalizationTreeReader(TTree &tree, TTreeReader &reader);
    ~CfNormalizationTreeReader();
    CfNormalizationResult Read() const;

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
}  // namespace exp_femto_3d
