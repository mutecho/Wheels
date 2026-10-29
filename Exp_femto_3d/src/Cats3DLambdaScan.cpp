#include "Cats3DLambdaScan.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace exp_femto_3d::cats3d_scan {
  std::vector<Point> Run(const std::vector<double> &coordinates, const std::vector<double> &seed,
                         const FitFunction &fit) {
    std::vector<double> sorted = coordinates;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end() ||
        !std::all_of(sorted.begin(), sorted.end(), [](double x) { return std::isfinite(x); })) {
      throw std::runtime_error("CATS3D scan coordinates must be unique and finite.");
    }
    std::vector<Point> points;
    for (double x : sorted) {
      Point p;
      p.coordinate = x;
      points.push_back(std::move(p));
    }
    auto attempt = [&](std::size_t i, const std::vector<double> &start, const char *origin) {
      auto result = fit(points[i].coordinate, start);
      result.seed_origin = origin;
      auto &point = points[i];
      const bool valid = profile_likelihood::IsValid(result.fit);
      const bool better =
          valid && (!point.Best() || result.fit.objective < point.Best()->fit.objective);
      point.attempts.push_back(std::move(result));
      if (better) {
        point.winner = static_cast<int>(point.attempts.size()) - 1;
      }
    };
    for (std::size_t i = 0; i < points.size(); ++i) {
      attempt(i, seed, "configured");
    }
    for (std::size_t i = 1; i < points.size(); ++i) {
      const auto *previous = points[i - 1].Best();
      if (previous) {
        const auto start = previous->fit.values;
        attempt(i, start, "forward_neighbor");
      }
    }
    for (std::size_t i = points.size(); i > 1; --i) {
      const auto *next = points[i - 1].Best();
      if (next) {
        const auto start = next->fit.values;
        attempt(i - 2, start, "reverse_neighbor");
      }
    }
    double reference = std::numeric_limits<double>::infinity();
    for (const auto &point : points) {
      if (point.Best()) {
        reference = std::min(reference, point.Best()->fit.objective);
      }
    }
    for (auto &point : points) {
      if (point.Best()) {
        point.delta_objective = point.Best()->fit.objective - reference;
      }
    }
    return points;
  }
  std::vector<GroupPoint> Sum(const std::vector<double> &coordinates,
                              const std::vector<const std::vector<Point> *> &members) {
    auto sorted = coordinates;
    std::sort(sorted.begin(), sorted.end());
    std::vector<GroupPoint> out;
    for (std::size_t i = 0; i < sorted.size(); ++i) {
      GroupPoint point;
      point.coordinate = sorted[i];
      point.valid = !members.empty();
      double sum = 0;
      for (const auto *member : members) {
        if (!member || member->size() != sorted.size() ||
            member->at(i).coordinate != point.coordinate) {
          throw std::runtime_error("CATS3D group scan coordinate/member mismatch.");
        }
        const auto *best = member->at(i).Best();
        if (!best) {
          point.valid = false;
        } else {
          sum += best->fit.objective;
        }
      }
      if (point.valid && std::isfinite(sum)) {
        point.objective = sum;
      } else {
        point.valid = false;
        point.failure_reason =
            members.empty() ? "no_differential_phi_members" : "one_or_more_phi_members_invalid";
      }
      out.push_back(std::move(point));
    }
    double reference = std::numeric_limits<double>::infinity();
    for (const auto &point : out) {
      if (point.valid) {
        reference = std::min(reference, point.objective);
      }
    }
    for (auto &point : out) {
      if (point.valid) {
        point.delta_objective = point.objective - reference;
      }
    }
    return out;
  }
} // namespace exp_femto_3d::cats3d_scan
