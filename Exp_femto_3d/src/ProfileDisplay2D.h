#pragma once

#include <cstddef>
#include <optional>
#include <vector>

namespace exp_femto_3d::profile_display_2d {

  struct Point2D {
    double x = 0.0;
    double y = 0.0;
  };

  struct Segment2D {
    Point2D first;
    Point2D second;
  };

  struct BestCandidate {
    double delta = 0.0;
    double objective = 0.0;
    double x = 0.0;
    double y = 0.0;
    int stage = 0;
    int point_index = 0;
    bool valid = false;
  };

  struct ContourLevel {
    double level = 0.0;
    std::vector<Segment2D> segments;
    std::size_t supported_cell_count = 0;
  };

  /**
   * @brief Build display-cell edges for ordered sampling coordinates.
   *
   * Interior edges are coordinate midpoints. The two outer edges are the
   * resolved scan bounds, so the persisted TH2 domain never extends beyond the
   * numerical scan contract. Sampling coordinates remain authoritative in
   * ProfilePoints and are not replaced by the resulting bin centers.
   */
  [[nodiscard]] std::vector<double> BuildClippedEdges(const std::vector<double>& coordinates,
                                                       double scan_lower,
                                                       double scan_upper);

  /** Select the finite valid minimum using delta, stage, then point index. */
  [[nodiscard]] std::optional<BestCandidate> SelectBestPoint(
      const std::vector<BestCandidate>& candidates);

  /**
   * @brief Generate contour segments on the true sampling-coordinate grid.
   *
   * Values and validity are row-major with x varying fastest. A cell
   * contributes only when all four corners are valid and finite. Saddle cells
   * use their bilinear center value for deterministic connectivity; equality is
   * assigned to the high side. Returned segments are disconnected primitives,
   * never a single polyline spanning missing regions.
   */
  [[nodiscard]] std::vector<ContourLevel> BuildContours(
      const std::vector<double>& x_coordinates,
      const std::vector<double>& y_coordinates,
      const std::vector<double>& values,
      const std::vector<bool>& valid,
      const std::vector<double>& levels);

}  // namespace exp_femto_3d::profile_display_2d
