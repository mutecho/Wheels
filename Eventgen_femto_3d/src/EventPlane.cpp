#include "femto3d/EventPlane.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace femto3d {

namespace {

constexpr double kMomentumTolerance = 1.0e-12;

void ValidateEventPlaneConfig(const EventPlaneConfig& config) {
  if (config.harmonic_order <= 0) {
    throw std::invalid_argument("Event-plane harmonic order must be positive.");
  }

  if (config.eta_max <= config.eta_min) {
    throw std::invalid_argument(
        "Event-plane eta maximum must be greater than minimum.");
  }

  if (config.min_q_magnitude < 0.0) {
    throw std::invalid_argument(
        "Event-plane minimum Q-vector magnitude must be non-negative.");
  }
}

double ComputeEventPlaneWeight(const ParticleState& particle,
                               const EventPlaneConfig& config) {
  switch (config.weight_mode) {
    case EventPlaneWeightMode::kUnit:
      return 1.0;
    case EventPlaneWeightMode::kPt:
      return ComputeParticlePt(particle);
  }

  throw std::invalid_argument("Unsupported event-plane weight mode.");
}

}  // namespace

double ComputeParticlePt(const ParticleState& particle) {
  return std::hypot(particle.px, particle.py);
}

double ComputeParticlePhi(const ParticleState& particle) {
  return std::atan2(particle.py, particle.px);
}

double ComputeParticleEta(const ParticleState& particle) {
  const double pt = ComputeParticlePt(particle);
  if (pt <= kMomentumTolerance) {
    if (std::abs(particle.pz) <= kMomentumTolerance) {
      return 0.0;
    }

    return std::copysign(std::numeric_limits<double>::infinity(), particle.pz);
  }

  return std::asinh(particle.pz / pt);
}

bool IsEventPlaneCandidate(const ParticleState& particle,
                           const EventPlaneConfig& config) {
  ValidateEventPlaneConfig(config);

  const double pt = ComputeParticlePt(particle);
  if (pt <= kMomentumTolerance) {
    return false;
  }

  const double eta = ComputeParticleEta(particle);
  if (!std::isfinite(eta) || eta < config.eta_min || eta >= config.eta_max) {
    return false;
  }

  if (config.allowed_abs_pdg.empty()) {
    return true;
  }

  const int abs_pdg = std::abs(particle.pdg);
  return std::find(config.allowed_abs_pdg.begin(),
                   config.allowed_abs_pdg.end(),
                   abs_pdg) != config.allowed_abs_pdg.end();
}

EventPlaneResult ReconstructEventPlane(
    const std::vector<ParticleState>& event_plane_candidates,
    const EventPlaneConfig& config) {
  ValidateEventPlaneConfig(config);

  EventPlaneResult result;
  result.candidate_count = event_plane_candidates.size();

  if (event_plane_candidates.empty()) {
    result.failure_reason = EventPlaneFailureReason::kNoCandidates;
    return result;
  }

  if (event_plane_candidates.size() < config.min_candidates) {
    result.failure_reason = EventPlaneFailureReason::kInsufficientCandidates;
    return result;
  }

  for (const ParticleState& particle : event_plane_candidates) {
    const double phi = ComputeParticlePhi(particle);
    const double weight = ComputeEventPlaneWeight(particle, config);
    result.qx += weight * std::cos(static_cast<double>(config.harmonic_order) * phi);
    result.qy += weight * std::sin(static_cast<double>(config.harmonic_order) * phi);
  }

  result.q_magnitude = std::hypot(result.qx, result.qy);
  if (result.q_magnitude < config.min_q_magnitude) {
    result.failure_reason = EventPlaneFailureReason::kQTooSmall;
    return result;
  }

  result.psi =
      WrapAngleToEventPlanePeriod(std::atan2(result.qy, result.qx) /
                                      static_cast<double>(config.harmonic_order),
                                  config.harmonic_order);
  result.success = true;
  return result;
}

EventPlaneResult MakeInputEventPlaneResult(const double event_plane_psi,
                                           const int harmonic_order) {
  EventPlaneResult result;
  result.success = true;
  result.psi = WrapAngleToEventPlanePeriod(event_plane_psi, harmonic_order);
  result.used_input_branch = true;
  return result;
}

double WrapAngleToEventPlanePeriod(double angle, const int harmonic_order) {
  if (harmonic_order <= 0) {
    throw std::invalid_argument("Event-plane harmonic order must be positive.");
  }

  const double period = 2.0 * kPi / static_cast<double>(harmonic_order);
  const double lower_bound = -0.5 * period;
  const double upper_bound = 0.5 * period;

  while (angle < lower_bound) {
    angle += period;
  }

  while (angle >= upper_bound) {
    angle -= period;
  }

  return angle;
}

double WrapPhiMinusPsi2(const double angle) {
  return WrapAngleToEventPlanePeriod(angle, 2);
}

double ComputePairPhiMinusEventPlane(const PairKinematics& pair_kinematics,
                                     const double event_plane_psi,
                                     const int harmonic_order) {
  return WrapAngleToEventPlanePeriod(pair_kinematics.pair_phi - event_plane_psi,
                                     harmonic_order);
}

double ComputePairPhiMinusPsi2(const PairKinematics& pair_kinematics,
                               const double event_plane_psi2) {
  return WrapPhiMinusPsi2(pair_kinematics.pair_phi - event_plane_psi2);
}

}  // namespace femto3d
