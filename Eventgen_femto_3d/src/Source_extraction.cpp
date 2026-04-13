#include <cmath>
#include <stdexcept>

#include "femto3d/SourceExtraction.h"

namespace femto3d {

  namespace {

    constexpr double kAxisTolerance = 1.0e-12;

    double ComputeSpatialSeparationSquared(const ParticleState &particle1, const ParticleState &particle2) {
      const double delta_x = particle1.x - particle2.x;
      const double delta_y = particle1.y - particle2.y;
      const double delta_z = particle1.z - particle2.z;
      return delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
    }

    double ResolveFemtoMtReferenceMass(const ParticleState &particle1,
                                       const ParticleState &particle2,
                                       const double configured_reference_mass,
                                       const double mass_tolerance) {
      if (mass_tolerance < 0.0) {
        throw std::invalid_argument("Femto mT mass tolerance must be non-negative.");
      }

      const double average_mass = 0.5 * (particle1.mass + particle2.mass);
      if (std::abs(particle1.mass - particle2.mass) > mass_tolerance) {
        throw std::invalid_argument("Femto mT requires a same-species pair with consistent masses.");
      }

      if (configured_reference_mass < 0.0) {
        return average_mass;
      }

      if (std::abs(average_mass - configured_reference_mass) > mass_tolerance) {
        throw std::invalid_argument("Configured femto mT reference mass is inconsistent with the selected pair.");
      }

      return configured_reference_mass;
    }

  }  // namespace

  double ComputeEnergy(const ParticleState &particle) {
    if (particle.mass < 0.0) {
      throw std::invalid_argument("Particle mass must be non-negative.");
    }

    return std::sqrt(particle.px * particle.px + particle.py * particle.py + particle.pz * particle.pz
                     + particle.mass * particle.mass);
  }

  PairKinematics ComputePairKinematics(const ParticleState &particle1,
                                       const ParticleState &particle2,
                                       const double femto_mt_reference_mass,
                                       const double femto_mt_mass_tolerance) {
    PairKinematics pair;
    pair.energy1 = ComputeEnergy(particle1);
    pair.energy2 = ComputeEnergy(particle2);
    pair.pair_px = particle1.px + particle2.px;
    pair.pair_py = particle1.py + particle2.py;
    pair.pair_pz = particle1.pz + particle2.pz;
    pair.pair_energy = pair.energy1 + pair.energy2;

    if (pair.pair_energy <= 0.0) {
      throw std::invalid_argument("Pair energy must be positive.");
    }

    pair.pair_pt = std::sqrt(pair.pair_px * pair.pair_px + pair.pair_py * pair.pair_py);
    const double pair_mt_squared = pair.pair_energy * pair.pair_energy - pair.pair_pz * pair.pair_pz;
    if (pair_mt_squared <= kAxisTolerance) {
      throw std::invalid_argument("Pair-system transverse mass is not positive, so the LCMS boost is undefined.");
    }

    const double reference_mass =
        ResolveFemtoMtReferenceMass(particle1, particle2, femto_mt_reference_mass, femto_mt_mass_tolerance);
    pair.pair_system_mt = std::sqrt(pair_mt_squared);
    pair.femto_kt = 0.5 * pair.pair_pt;
    pair.femto_mt = std::sqrt(reference_mass * reference_mass + pair.femto_kt * pair.femto_kt);
    pair.pair_phi = std::atan2(pair.pair_py, pair.pair_px);
    return pair;
  }

  double WrapPhiToPi(double angle) {
    while (angle <= -kPi) {
      angle += 2.0 * kPi;
    }

    while (angle > kPi) {
      angle -= 2.0 * kPi;
    }

    return angle;
  }

  double ComputePairPhiMinusPsi(const PairKinematics &pair_kinematics, const double event_plane_psi) {
    return WrapPhiToPi(pair_kinematics.pair_phi - event_plane_psi);
  }

  // Sign convention: all relative quantities are defined as particle1 - particle2.
  PairSeparationLCMS ComputeSinglePairDistanceLCMS(const ParticleState &particle1,
                                                   const ParticleState &particle2,
                                                   const double close_pair_distance_tolerance) {
    if (close_pair_distance_tolerance < 0.0) {
      throw std::invalid_argument("Close pair distance tolerance must be non-negative.");
    }

    const double spatial_separation_squared = ComputeSpatialSeparationSquared(particle1, particle2);
    if (spatial_separation_squared < close_pair_distance_tolerance * close_pair_distance_tolerance) {
      throw std::invalid_argument("Close pair rejection: the spatial separation is too small.");
    }

    const PairKinematics pair = ComputePairKinematics(particle1, particle2);
    if (pair.pair_pt < kAxisTolerance) {
      throw std::invalid_argument("Pair transverse momentum is zero, so the out-side axes are undefined.");
    }

    const double beta_lcms = pair.pair_pz / pair.pair_energy;
    const double gamma_lcms = pair.pair_energy / pair.pair_system_mt;

    const double out_x = pair.pair_px / pair.pair_pt;
    const double out_y = pair.pair_py / pair.pair_pt;
    const double side_x = -out_y;
    const double side_y = out_x;

    const double delta_x = particle1.x - particle2.x;
    const double delta_y = particle1.y - particle2.y;
    const double delta_z = particle1.z - particle2.z;
    const double delta_t = particle1.t - particle2.t;

    PairSeparationLCMS separation;
    separation.rho_out = delta_x * out_x + delta_y * out_y;
    separation.rho_side = delta_x * side_x + delta_y * side_y;
    separation.rho_long = gamma_lcms * (delta_z - beta_lcms * delta_t);
    return separation;
  }

  PairSeparationLCMS ComputeSinglePairDistanceLCMS(double px1,
                                                   double py1,
                                                   double pz1,
                                                   double mass1,
                                                   double x1,
                                                   double y1,
                                                   double z1,
                                                   double t1,
                                                   double px2,
                                                   double py2,
                                                   double pz2,
                                                   double mass2,
                                                   double x2,
                                                   double y2,
                                                   double z2,
                                                   double t2,
                                                   const double close_pair_distance_tolerance) {
    ParticleState particle1;
    particle1.px = px1;
    particle1.py = py1;
    particle1.pz = pz1;
    particle1.mass = mass1;
    particle1.x = x1;
    particle1.y = y1;
    particle1.z = z1;
    particle1.t = t1;

    ParticleState particle2;
    particle2.px = px2;
    particle2.py = py2;
    particle2.pz = pz2;
    particle2.mass = mass2;
    particle2.x = x2;
    particle2.y = y2;
    particle2.z = z2;
    particle2.t = t2;

    return ComputeSinglePairDistanceLCMS(particle1, particle2, close_pair_distance_tolerance);
  }

}  // namespace femto3d
