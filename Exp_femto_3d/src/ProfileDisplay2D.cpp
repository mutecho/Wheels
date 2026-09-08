#include "ProfileDisplay2D.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>

namespace exp_femto_3d::profile_display_2d {

  namespace {

    constexpr double kGeometryTolerance = 1.0e-12;

    bool NearlyEqual(const double left, const double right) {
      return std::abs(left - right)
             <= kGeometryTolerance * (1.0 + std::max(std::abs(left), std::abs(right)));
    }

    bool SamePoint(const Point2D& left, const Point2D& right) {
      return NearlyEqual(left.x, right.x) && NearlyEqual(left.y, right.y);
    }

    bool PointLess(const Point2D& left, const Point2D& right) {
      if (!NearlyEqual(left.x, right.x)) return left.x < right.x;
      return left.y < right.y && !NearlyEqual(left.y, right.y);
    }

    Point2D Interpolate(const Point2D& low,
                        const Point2D& high,
                        const double low_value,
                        const double high_value,
                        const double level) {
      const double denominator = high_value - low_value;
      const double fraction = denominator == 0.0 ? 0.5 : (level - low_value) / denominator;
      const double clipped = std::max(0.0, std::min(1.0, fraction));
      return {low.x + clipped * (high.x - low.x), low.y + clipped * (high.y - low.y)};
    }

    void AddSegment(ContourLevel& output, Point2D first, Point2D second) {
      if (SamePoint(first, second)) return;
      if (PointLess(second, first)) std::swap(first, second);
      const bool duplicate = std::any_of(output.segments.begin(), output.segments.end(),
                                         [&](const Segment2D& existing) {
                                           return SamePoint(existing.first, first)
                                                  && SamePoint(existing.second, second);
                                         });
      if (!duplicate) output.segments.push_back({first, second});
    }

  }  // namespace

  std::vector<double> BuildClippedEdges(const std::vector<double>& coordinates,
                                        const double scan_lower,
                                        const double scan_upper) {
    if (coordinates.empty()) throw std::invalid_argument("2D display axis must contain sampling coordinates");
    if (!std::isfinite(scan_lower) || !std::isfinite(scan_upper) || !(scan_upper > scan_lower)) {
      throw std::invalid_argument("2D display scan bounds must be finite and ordered");
    }
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
      if (!std::isfinite(coordinates[index])) {
        throw std::invalid_argument("2D display sampling coordinates must be finite");
      }
      if (index > 0U && !(coordinates[index] > coordinates[index - 1U])) {
        throw std::invalid_argument("2D display sampling coordinates must be strictly increasing");
      }
    }
    if (coordinates.front() < scan_lower && !NearlyEqual(coordinates.front(), scan_lower)) {
      throw std::invalid_argument("first 2D sampling coordinate lies below the resolved scan bound");
    }
    if (coordinates.back() > scan_upper && !NearlyEqual(coordinates.back(), scan_upper)) {
      throw std::invalid_argument("last 2D sampling coordinate lies above the resolved scan bound");
    }

    std::vector<double> edges(coordinates.size() + 1U);
    edges.front() = scan_lower;
    edges.back() = scan_upper;
    for (std::size_t index = 1; index < coordinates.size(); ++index) {
      edges[index] = 0.5 * (coordinates[index - 1U] + coordinates[index]);
    }
    for (std::size_t index = 1; index < edges.size(); ++index) {
      if (!(edges[index] > edges[index - 1U])) {
        throw std::invalid_argument("2D display cell edges must be strictly increasing");
      }
    }
    return edges;
  }

  std::optional<BestCandidate> SelectBestPoint(const std::vector<BestCandidate>& candidates) {
    std::optional<BestCandidate> best;
    for (const BestCandidate& candidate : candidates) {
      if (!candidate.valid || !std::isfinite(candidate.delta)) continue;
      if (!best || std::tie(candidate.delta, candidate.stage, candidate.point_index)
                       < std::tie(best->delta, best->stage, best->point_index)) {
        best = candidate;
      }
    }
    return best;
  }

  std::vector<ContourLevel> BuildContours(const std::vector<double>& x_coordinates,
                                          const std::vector<double>& y_coordinates,
                                          const std::vector<double>& values,
                                          const std::vector<bool>& valid,
                                          const std::vector<double>& levels) {
    if (x_coordinates.size() < 2U || y_coordinates.size() < 2U) {
      throw std::invalid_argument("2D contours require at least two coordinates on each axis");
    }
    const std::size_t expected = x_coordinates.size() * y_coordinates.size();
    if (values.size() != expected || valid.size() != expected) {
      throw std::invalid_argument("2D contour values and validity must match the coordinate grid");
    }
    for (std::size_t index = 1; index < x_coordinates.size(); ++index) {
      if (!(x_coordinates[index] > x_coordinates[index - 1U])) {
        throw std::invalid_argument("2D contour X coordinates must be strictly increasing");
      }
    }
    for (std::size_t index = 1; index < y_coordinates.size(); ++index) {
      if (!(y_coordinates[index] > y_coordinates[index - 1U])) {
        throw std::invalid_argument("2D contour Y coordinates must be strictly increasing");
      }
    }

    std::vector<ContourLevel> output;
    output.reserve(levels.size());
    const std::size_t nx = x_coordinates.size();
    const auto flat = [nx](const std::size_t ix, const std::size_t iy) { return iy * nx + ix; };
    for (const double level : levels) {
      ContourLevel contour;
      contour.level = level;
      if (!std::isfinite(level)) {
        output.push_back(std::move(contour));
        continue;
      }
      for (std::size_t iy = 0; iy + 1U < y_coordinates.size(); ++iy) {
        for (std::size_t ix = 0; ix + 1U < x_coordinates.size(); ++ix) {
          const std::size_t indices[4] = {
              flat(ix, iy), flat(ix + 1U, iy), flat(ix + 1U, iy + 1U), flat(ix, iy + 1U)};
          bool supported = true;
          for (const std::size_t index : indices) {
            supported = supported && valid[index] && std::isfinite(values[index]);
          }
          if (!supported) continue;
          ++contour.supported_cell_count;

          const Point2D corners[4] = {
              {x_coordinates[ix], y_coordinates[iy]},
              {x_coordinates[ix + 1U], y_coordinates[iy]},
              {x_coordinates[ix + 1U], y_coordinates[iy + 1U]},
              {x_coordinates[ix], y_coordinates[iy + 1U]}};
          const double corner_values[4] = {
              values[indices[0]], values[indices[1]], values[indices[2]], values[indices[3]]};
          const bool high[4] = {corner_values[0] >= level, corner_values[1] >= level,
                                corner_values[2] >= level, corner_values[3] >= level};
          std::vector<int> crossing_edges;
          for (int edge = 0; edge < 4; ++edge) {
            const int next = (edge + 1) % 4;
            if (high[edge] != high[next]) crossing_edges.push_back(edge);
          }
          if (crossing_edges.empty()) continue;

          Point2D intersections[4];
          bool have_intersection[4] = {false, false, false, false};
          for (const int edge : crossing_edges) {
            const int next = (edge + 1) % 4;
            intersections[edge] = Interpolate(corners[edge], corners[next],
                                              corner_values[edge], corner_values[next], level);
            have_intersection[edge] = true;
          }
          if (crossing_edges.size() == 2U) {
            AddSegment(contour, intersections[crossing_edges[0]], intersections[crossing_edges[1]]);
            continue;
          }
          if (crossing_edges.size() != 4U || !have_intersection[0] || !have_intersection[1]
              || !have_intersection[2] || !have_intersection[3]) {
            continue;
          }

          const int pattern = (high[0] ? 1 : 0) | (high[1] ? 2 : 0)
                              | (high[2] ? 4 : 0) | (high[3] ? 8 : 0);
          const double center_value = 0.25 * (corner_values[0] + corner_values[1]
                                              + corner_values[2] + corner_values[3]);
          const bool center_high = center_value >= level;
          if ((pattern == 5 && center_high) || (pattern == 10 && !center_high)) {
            AddSegment(contour, intersections[0], intersections[1]);
            AddSegment(contour, intersections[2], intersections[3]);
          } else {
            AddSegment(contour, intersections[3], intersections[0]);
            AddSegment(contour, intersections[1], intersections[2]);
          }
        }
      }
      output.push_back(std::move(contour));
    }
    return output;
  }

}  // namespace exp_femto_3d::profile_display_2d
