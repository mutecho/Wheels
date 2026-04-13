#pragma once

#include "femto3d/AnalysisTypes.h"

namespace femto3d {

  [[nodiscard]] double ComputeEnergy(const ParticleState &particle);

  [[nodiscard]] PairKinematics ComputePairKinematics(const ParticleState &particle1,
                                                     const ParticleState &particle2,
                                                     double femto_mt_reference_mass = -1.0,
                                                     double femto_mt_mass_tolerance = 1.0e-3);

  [[nodiscard]] double WrapPhiToPi(double angle);

  [[nodiscard]] double ComputePairPhiMinusPsi(const PairKinematics &pair_kinematics, double event_plane_psi);

  [[nodiscard]] PairSeparationLCMS ComputeSinglePairDistanceLCMS(const ParticleState &particle1,
                                                                 const ParticleState &particle2,
                                                                 double close_pair_distance_tolerance = 1.0e-12);

  [[nodiscard]] PairSeparationLCMS ComputeSinglePairDistanceLCMS(double px1,
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
                                                                 double close_pair_distance_tolerance = 1.0e-12);

}  // namespace femto3d
