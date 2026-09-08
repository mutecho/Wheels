#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "ProfileDisplay2D.h"

namespace {

  using namespace exp_femto_3d::profile_display_2d;

  void Expect(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
  }

  bool Close(const double left, const double right) {
    return std::abs(left - right) < 1.0e-12;
  }

  bool HasEndpoint(const Segment2D& segment, const double x, const double y) {
    return (Close(segment.first.x, x) && Close(segment.first.y, y))
           || (Close(segment.second.x, x) && Close(segment.second.y, y));
  }

}  // namespace

int main() {
  using namespace exp_femto_3d::profile_display_2d;

  const auto regular_edges = BuildClippedEdges({0.0, 1.0, 2.0}, 0.0, 2.0);
  Expect(regular_edges == std::vector<double>({0.0, 0.5, 1.5, 2.0}),
         "regular display edges must use midpoint interiors and exact scan bounds");
  const auto irregular_edges = BuildClippedEdges({0.0, 1.0, 4.0}, 0.0, 4.0);
  Expect(irregular_edges == std::vector<double>({0.0, 0.5, 2.5, 4.0}),
         "irregular display edges must retain midpoint interiors");

  const std::vector<double> x{0.0, 2.0, 5.0};
  const std::vector<double> y{10.0, 12.0, 16.0};
  std::vector<double> linear;
  for (const double yy : y) {
    for (const double xx : x) linear.push_back(xx + yy);
  }
  const auto linear_contours = BuildContours(x, y, linear, std::vector<bool>(9, true), {14.0, 30.0});
  Expect(linear_contours[0].supported_cell_count == 4U && !linear_contours[0].segments.empty(),
         "linear field must generate supported contour segments");
  Expect(linear_contours[1].supported_cell_count == 4U && linear_contours[1].segments.empty(),
         "uncrossed level must not generate a contour");
  bool uses_true_coordinate = false;
  for (const Segment2D& segment : linear_contours[0].segments) {
    uses_true_coordinate = uses_true_coordinate || HasEndpoint(segment, 2.0, 12.0);
  }
  Expect(uses_true_coordinate, "contours must use true sampling coordinates rather than TH2 bin centers");

  std::vector<bool> hole_valid(9, true);
  hole_valid[4] = false;
  const auto hole = BuildContours(x, y, linear, hole_valid, {14.0});
  Expect(hole[0].supported_cell_count == 0U && hole[0].segments.empty(),
         "an invalid corner must suppress every cell that touches it");

  const std::vector<double> disconnected_x{0.0, 1.0, 2.0, 3.0, 4.0};
  const std::vector<double> disconnected_y{0.0, 1.0};
  const std::vector<double> disconnected_values{0.0, 2.0, 0.0, 0.0, 2.0,
                                                0.0, 2.0, 0.0, 0.0, 2.0};
  const std::vector<bool> disconnected_valid{true, true, false, true, true,
                                             true, true, false, true, true};
  const auto disconnected = BuildContours(disconnected_x, disconnected_y,
                                          disconnected_values, disconnected_valid, {1.0});
  Expect(disconnected[0].supported_cell_count == 2U && disconnected[0].segments.size() == 2U,
         "disconnected valid regions must remain separate contour segments");

  const auto constant = BuildContours({0.0, 1.0}, {0.0, 1.0},
                                      {2.0, 2.0, 2.0, 2.0},
                                      std::vector<bool>(4, true), {2.0});
  Expect(constant[0].supported_cell_count == 1U && constant[0].segments.empty(),
         "an all-equal cell must not generate zero-length contour segments");

  const auto vertex_equal = BuildContours({0.0, 1.0}, {0.0, 1.0},
                                          {1.0, 0.0, 2.0, 2.0},
                                          std::vector<bool>(4, true), {1.0});
  Expect(vertex_equal[0].segments.size() == 1U,
         "a vertex on the level must follow the deterministic high-side rule");
  Expect(!Close(vertex_equal[0].segments[0].first.x, vertex_equal[0].segments[0].second.x)
             || !Close(vertex_equal[0].segments[0].first.y, vertex_equal[0].segments[0].second.y),
         "vertex equality must not leave a zero-length segment");

  const auto saddle_high = BuildContours({0.0, 1.0}, {0.0, 1.0},
                                         {3.0, 0.0, 0.0, 3.0},
                                         std::vector<bool>(4, true), {1.0});
  const auto saddle_low = BuildContours({0.0, 1.0}, {0.0, 1.0},
                                        {2.0, 0.0, 0.0, 2.0},
                                        std::vector<bool>(4, true), {1.5});
  Expect(saddle_high[0].segments.size() == 2U && saddle_low[0].segments.size() == 2U,
         "saddle cells must resolve into two non-duplicated segments");
  Expect(HasEndpoint(saddle_high[0].segments[0], 2.0 / 3.0, 0.0)
             || HasEndpoint(saddle_high[0].segments[1], 2.0 / 3.0, 0.0),
         "saddle interpolation must preserve the analytic edge crossing");

  const std::vector<BestCandidate> candidates{
      {1.0, 11.0, 0.0, 0.0, 1, 8, true},
      {1.0, 10.0, 1.0, 1.0, 0, 9, true},
      {1.0, 9.0, 2.0, 2.0, 0, 4, true},
      {0.5, 8.0, 3.0, 3.0, 1, 20, true},
      {0.1, 7.0, 4.0, 4.0, 0, 1, false},
      {std::numeric_limits<double>::quiet_NaN(), 6.0, 5.0, 5.0, 0, 0, true}};
  const auto best = SelectBestPoint(candidates);
  Expect(best && best->point_index == 20 && best->stage == 1 && Close(best->delta, 0.5),
         "a refined candidate with the smallest finite delta must win");
  const auto tie = SelectBestPoint({candidates[0], candidates[1], candidates[2]});
  Expect(tie && tie->point_index == 4 && tie->stage == 0,
         "equal deltas must use stage and then point index as stable tie-breaks");
  Expect(!SelectBestPoint({candidates[4], candidates[5]}),
         "no best point may be fabricated when every candidate is unavailable");

  std::cout << "profile_display_2d_test passed\n";
  return 0;
}
