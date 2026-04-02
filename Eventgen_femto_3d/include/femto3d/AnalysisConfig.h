#pragma once

#include "femto3d/AnalysisTypes.h"
#include "femto3d/ProjectionFit.h"

#include <string>
#include <vector>

namespace femto3d {

enum class EventPlaneWeightMode {
  kUnit = 0,
  kPt = 1,
};

struct InputTreeConfig {
  std::string tree_name = "events";
  std::string centrality_branch = "centrality";
  std::string event_plane_branch = "event_plane_psi";
  std::string pdg_branch = "pdg";
  std::string px_branch = "px";
  std::string py_branch = "py";
  std::string pz_branch = "pz";
  std::string mass_branch = "mass";
  std::string x_branch = "x";
  std::string y_branch = "y";
  std::string z_branch = "z";
  std::string t_branch = "t";
};

struct EventPlaneConfig {
  bool enabled = true;
  int harmonic_order = 2;
  bool use_internal_reconstruction = true;
  bool fallback_to_input_branch = false;
  double eta_min = -2.5;
  double eta_max = -0.5;
  EventPlaneWeightMode weight_mode = EventPlaneWeightMode::kPt;
  std::vector<int> allowed_abs_pdg = {211, 321, 2212};
  std::size_t min_candidates = 2U;
  double min_q_magnitude = 1.0e-6;
};

struct ParticleSelectionConfig {
  int target_pdg = 211;
  double femto_eta_min = -0.5;
  double femto_eta_max = 0.5;
  double femto_pt_min = 0.2;
  double femto_pt_max = 1.0;
  // Reference mass used in femto_mt = sqrt(m_ref^2 + kT^2).
  double femto_mt_reference_mass = kChargedPionMass;
  // Keep target_pdg and m_ref aligned to the same particle species.
  double femto_mt_mass_tolerance = 1.0e-3;
};

struct HistogramConfig {
  AxisSpec rho_out_axis{"#rho_{out}", -20.0, 20.0, 0.5};
  AxisSpec rho_side_axis{"#rho_{side}", -20.0, 20.0, 0.5};
  AxisSpec rho_long_axis{"#rho_{long}", -20.0, 20.0, 0.5};
  AxisSpec projection_axis{"r", -20.0, 20.0, 0.5};
  bool warn_on_overflow = true;
};

struct AnalysisConfig {
  InputTreeConfig input;
  EventPlaneConfig event_plane;
  ParticleSelectionConfig selection;
  HistogramConfig histograms;
  ProjectionFitConfig projection_fit;
  std::vector<RangeBin> centrality_bins;
  // Standard femtoscopic mT slice bins: sqrt(m_ref^2 + kT^2).
  std::vector<RangeBin> mt_bins;
  std::vector<RangeBin> phi_bins;
  double close_pair_distance_tolerance = 1.0e-12;
  std::string top_directory_name = "Femto3D";
  std::string r2_summary_directory_name = "R2Summary";
};

inline std::vector<RangeBin> MakeUniformBins(const double min,
                                             const double max,
                                             const std::size_t count,
                                             const std::string& prefix) {
  if (max <= min) {
    throw std::invalid_argument("Uniform bin maximum must be greater than minimum.");
  }

  if (count == 0U) {
    throw std::invalid_argument("Uniform bin count must be positive.");
  }

  std::vector<RangeBin> bins;
  bins.reserve(count);

  const double width = (max - min) / static_cast<double>(count);
  for (std::size_t index = 0; index < count; ++index) {
    const double low = min + static_cast<double>(index) * width;
    const double high = low + width;
    bins.push_back({low, high, FormatBinLabel(prefix, low, high)});
  }

  return bins;
}

inline AnalysisConfig MakeDefaultAnalysisConfig() {
  AnalysisConfig config;
  config.centrality_bins = {
      {0.0, 10.0, FormatBinLabel("cent", 0.0, 10.0)},
      {10.0, 30.0, FormatBinLabel("cent", 10.0, 30.0)},
      {30.0, 50.0, FormatBinLabel("cent", 30.0, 50.0)},
  };
  config.mt_bins = {
      {0.20, 0.40, FormatBinLabel("femto_mt", 0.20, 0.40)},
      {0.40, 0.60, FormatBinLabel("femto_mt", 0.40, 0.60)},
      {0.60, 0.80, FormatBinLabel("femto_mt", 0.60, 0.80)},
  };
  config.phi_bins = MakeUniformBins(-kPi / 2.0, kPi / 2.0, 12U, "phi");
  return config;
}

}  // namespace femto3d
