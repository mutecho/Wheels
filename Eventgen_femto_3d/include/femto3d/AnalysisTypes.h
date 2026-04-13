#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace femto3d {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kSqrtHalf = 0.70710678118654752440;
  constexpr double kChargedPionMass = 0.13957039;

  struct ParticleState {
    int pdg = 0;
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double mass = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double t = 0.0;
  };

  struct EventData {
    double centrality = 0.0;
    double event_plane_psi = 0.0;
    bool has_input_event_plane = false;
    std::vector<ParticleState> particles;
  };

  struct PairKinematics {
    double energy1 = 0.0;
    double energy2 = 0.0;
    double pair_px = 0.0;
    double pair_py = 0.0;
    double pair_pz = 0.0;
    double pair_energy = 0.0;
    double pair_pt = 0.0;
    double pair_system_mt = 0.0;
    double femto_kt = 0.0;
    double femto_mt = 0.0;
    double pair_phi = 0.0;
  };

  struct PairSeparationLCMS {
    double rho_out = 0.0;
    double rho_side = 0.0;
    double rho_long = 0.0;
  };

  struct RangeBin {
    double min = 0.0;
    double max = 0.0;
    std::string label;

    [[nodiscard]] bool Contains(const double value) const {
      return value >= min && value < max;
    }

    [[nodiscard]] double Center() const {
      return 0.5 * (min + max);
    }
  };

  struct AxisSpec {
    std::string title;
    double min = 0.0;
    double max = 0.0;
    double bin_width = 1.0;

    [[nodiscard]] int NumBins() const {
      if (max <= min) {
        throw std::invalid_argument("Axis maximum must be greater than minimum.");
      }

      if (bin_width <= 0.0) {
        throw std::invalid_argument("Axis bin width must be positive.");
      }

      return static_cast<int>(std::ceil((max - min) / bin_width));
    }
  };

  struct SliceKey {
    std::size_t cent_index = 0;
    std::size_t mt_index = 0;
    std::size_t phi_index = 0;

    [[nodiscard]] bool operator<(const SliceKey &other) const {
      return std::tie(cent_index, mt_index, phi_index) < std::tie(other.cent_index, other.mt_index, other.phi_index);
    }
  };

  struct ProjectionDefinition {
    std::string short_name;
    std::string title;
    double n_out = 0.0;
    double n_side = 0.0;
    double n_long = 0.0;
  };

  struct ProjectionFitResult {
    std::string projection_name;
    double amplitude = 0.0;
    double amplitude_error = 0.0;
    double r2 = 0.0;
    double r2_error = 0.0;
    bool success = false;
  };

  struct AnalysisStatistics {
    std::size_t events_read = 0;
    std::size_t selected_particles = 0;
    std::size_t events_with_valid_event_plane = 0;
    std::size_t events_rejected_no_event_plane_candidates = 0;
    std::size_t events_rejected_small_q_vector = 0;
    std::size_t events_rejected_missing_input_event_plane = 0;
    std::size_t events_rejected_insufficient_femto_particles = 0;
    std::size_t candidate_pairs = 0;
    std::size_t accepted_pairs = 0;
    std::size_t rejected_close_pairs = 0;
    std::size_t rejected_invalid_pair_kinematics = 0;
    std::size_t r2_summary_points_skipped_invalid_hbt_error = 0;
  };

  inline std::string FormatBinLabel(const std::string &prefix,
                                    const double min,
                                    const double max,
                                    const int precision = 2) {
    std::ostringstream stream;
    stream << prefix << "_" << std::fixed << std::setprecision(precision) << min << "_" << max;
    return stream.str();
  }

  inline std::array<ProjectionDefinition, 6> MakeProjectionDefinitions() {
    return {{
        {"out", "r_{out}", 1.0, 0.0, 0.0},
        {"side", "r_{side}", 0.0, 1.0, 0.0},
        {"long", "r_{long}", 0.0, 0.0, 1.0},
        {"out_side", "(r_{out}+r_{side})/#sqrt{2}", kSqrtHalf, kSqrtHalf, 0.0},
        {"out_long", "(r_{out}+r_{long})/#sqrt{2}", kSqrtHalf, 0.0, kSqrtHalf},
        {"side_long", "(r_{side}+r_{long})/#sqrt{2}", 0.0, kSqrtHalf, kSqrtHalf},
    }};
  }

  inline double ProjectOntoDirection(const PairSeparationLCMS &separation, const ProjectionDefinition &direction) {
    return direction.n_out * separation.rho_out + direction.n_side * separation.rho_side
           + direction.n_long * separation.rho_long;
  }

}  // namespace femto3d
