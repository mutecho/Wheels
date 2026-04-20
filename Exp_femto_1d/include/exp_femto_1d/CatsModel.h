#pragma once

#include <string>

#include "exp_femto_1d/Types.h"

class TF1;

namespace exp_femto_1d {

  class PiPiCatsModel {
   public:
    PiPiCatsModel(unsigned num_mom_bins, double k_min_mev, double k_max_mev, bool use_coulomb);
    ~PiPiCatsModel();

    [[nodiscard]] double Evaluate(double kstar_mev, double source_size_fm);
    [[nodiscard]] TF1 *BuildFitFunction(const std::string &name, const FitConfig &config);

   private:
    class Impl;
    Impl *impl_ = nullptr;
  };

}  // namespace exp_femto_1d
