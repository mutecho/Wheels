#pragma once

#include "femto3d/AnalysisConfig.h"
#include "femto3d/AnalysisTypes.h"

#include <cstddef>
#include <vector>

namespace femto3d {

enum class EventPlaneFailureReason {
  kNone = 0,
  kDisabled = 1,
  kNoCandidates = 2,
  kInsufficientCandidates = 3,
  kQTooSmall = 4,
  kMissingInputBranch = 5,
};

struct EventPlaneResult {
  bool success = false;
  double psi = 0.0;
  double qx = 0.0;
  double qy = 0.0;
  double q_magnitude = 0.0;
  std::size_t candidate_count = 0;
  EventPlaneFailureReason failure_reason = EventPlaneFailureReason::kNone;
  bool used_input_branch = false;
};

[[nodiscard]] double ComputeParticlePt(const ParticleState& particle);

[[nodiscard]] double ComputeParticlePhi(const ParticleState& particle);

[[nodiscard]] double ComputeParticleEta(const ParticleState& particle);

[[nodiscard]] bool IsEventPlaneCandidate(const ParticleState& particle,
                                         const EventPlaneConfig& config);

[[nodiscard]] EventPlaneResult ReconstructEventPlane(
    const std::vector<ParticleState>& event_plane_candidates,
    const EventPlaneConfig& config);

[[nodiscard]] EventPlaneResult MakeInputEventPlaneResult(double event_plane_psi,
                                                         int harmonic_order);

[[nodiscard]] double WrapAngleToEventPlanePeriod(double angle,
                                                 int harmonic_order);

[[nodiscard]] double WrapPhiMinusPsi2(double angle);

[[nodiscard]] double ComputePairPhiMinusEventPlane(
    const PairKinematics& pair_kinematics,
    double event_plane_psi,
    int harmonic_order);

[[nodiscard]] double ComputePairPhiMinusPsi2(const PairKinematics& pair_kinematics,
                                             double event_plane_psi2);

}  // namespace femto3d
