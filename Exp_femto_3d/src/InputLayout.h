#pragma once

#include <array>
#include <memory>
#include <string>

#include "TAxis.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  /** Semantic indices in the producer's seven-dimensional LCMS count histogram.
   * No axis reordering is performed; selections and projections share this layout.
   */
  struct InputAxisLayout {
    std::array<int, 3> q{0, 1, 2};
    int mt = 3;
    int centrality = 4;
    int phi = 6;
    int qn = 5;
  };

  [[nodiscard]] InputAxisLayout ResolveInputAxisLayout(InputFormat format);

  /** Validate dimensions before indexing; Femto additionally requires all SE/ME
   * axis edges to be finite, increasing and equal. Diagnostics retain both paths.
   */
  void ValidateInputSparses(const THnSparseF &se, const THnSparseF &me, const InputConfig &config,
                            const InputAxisLayout &layout, const std::string &se_path,
                            const std::string &me_path);

  /** The existing symmetric label mapping is defined only on raw [0, pi].
   * Reject wider domains before opening build/fit outputs; this never folds counts.
   */
  void ValidateSymmetricPhiDomain(double low, double high, bool mapped);

  /** Effective phi bins used by rebin/selection. fold_pi requires a complete
   * paired 0..2pi axis; the returned axis contains its first half only.
   */
  [[nodiscard]] TAxis ResolvePhiSelectionAxis(const TAxis &source, PhiMode mode);

  /** Project the active selection; fold_pi adds the pi-shifted projection in
   * raw counts and sum(w^2), then restores the selected source phi range.
   */
  [[nodiscard]] std::unique_ptr<TH3D> ProjectInput3D(THnSparseF &source, const InputAxisLayout &layout,
                                                     PhiMode mode);

}  // namespace exp_femto_3d
