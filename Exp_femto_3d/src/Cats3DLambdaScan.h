#pragma once
#include "ProfileLikelihood.h"
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace exp_femto_3d::cats3d_scan {
  /** One actual minimizer invocation; unsuccessful trials remain in the record. */
  struct Attempt {
    std::string seed_origin;
    std::string failure_reason;
    profile_likelihood::MinimizationResult fit;
    int ndf = -1;
  };
  struct Point {
    double coordinate = 0;
    std::vector<Attempt> attempts;
    int winner = -1;
    double delta_objective = std::numeric_limits<double>::quiet_NaN();
    const Attempt *Best() const { return winner < 0 ? nullptr : &attempts.at(winner); }
  };
  using FitFunction = std::function<Attempt(double, const std::vector<double> &)>;
  // Coordinates are evaluated in sorted order, with configured, forward
  // neighbor and reverse neighbor starts. Only valid MIGRAD minima can win.
  std::vector<Point> Run(const std::vector<double> &coordinates,
                         const std::vector<double> &configured_seed, const FitFunction &fit);
  struct GroupPoint {
    double coordinate = 0;
    double objective = std::numeric_limits<double>::quiet_NaN();
    double delta_objective = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
    std::string failure_reason;
  };
  // Pass only differential-phi members. A missing/failed member invalidates the
  // whole coordinate; never drop a member and sum a smaller dataset.
  std::vector<GroupPoint> Sum(const std::vector<double> &coordinates,
                              const std::vector<const std::vector<Point> *> &members);
} // namespace exp_femto_3d::cats3d_scan
