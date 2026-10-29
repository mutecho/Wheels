#include "Cats3DLambdaScan.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace exp_femto_3d;
static void Require(bool ok, const char *message) {
  if (!ok) {
    throw std::runtime_error(message);
  }
}
int main() {
  try {
    const std::vector<double> coordinates = {.3, .5, .7, .9};
    int calls = 0;
    const auto fit = [&](double x, const std::vector<double> &seed) {
      ++calls;
      cats3d_scan::Attempt a;
      a.fit.values = {2 * x, seed.at(0)};
      a.fit.errors = {.1, .2};
      a.fit.objective = (x - .6) * (x - .6);
      a.fit.migrad_status = 0;
      a.fit.model_domain_valid = x != .5;
      a.fit.objective_valid = true;
      a.fit.parameter_errors_valid = true;
      if (x == .5) {
        a.failure_reason = "synthetic_invalid_domain";
      }
      return a;
    };
    const auto points = cats3d_scan::Run(coordinates, {1}, fit);
    Require(points.size() == 4 && !points[1].Best(),
            "invalid point disappeared or acquired a winner");
    Require(std::isnan(points[1].delta_objective), "invalid delta was fabricated");
    for (const auto &point : points) {
      Require(point.attempts.front().seed_origin == "configured" &&
                  point.attempts.front().fit.values.at(1) == 1,
              "configured seed was not attempted independently at every coordinate");
    }

    Require(points[0].Best()->fit.values[0] == .6 && points[3].Best()->fit.values[0] == 1.8,
            "nuisance was not re-minimized");
    Require(calls > 4 && points[2].attempts.size() > 1, "neighbor retries missing");
    auto sum = cats3d_scan::Sum(coordinates, {&points, &points});
    Require(!sum[1].valid && std::isnan(sum[1].objective),
            "failed group member was dropped from sum");
    Require(std::abs(sum[0].objective - 2 * points[0].Best()->fit.objective) < 1e-14,
            "group objective is not a sum");
    auto empty = cats3d_scan::Sum(coordinates, {});
    Require(!empty[0].valid && empty[0].failure_reason == "no_differential_phi_members",
            "empty group was fabricated");
    bool rejected = false;
    try {
      cats3d_scan::Run({.5, .5}, {1}, fit);
    } catch (const std::exception &) {
      rejected = true;
    }
    Require(rejected, "duplicate coordinates accepted");
    std::cout << "CATS3D scan retries, invalid gaps and complete-member sums passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
