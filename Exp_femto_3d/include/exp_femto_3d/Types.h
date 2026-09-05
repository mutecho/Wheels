#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace exp_femto_3d {

  constexpr double kPi = 3.14159265358979323846;

  enum class LogLevel {
    kDebug = 0,
    kInfo = 1,
    kWarn = 2,
    kError = 3,
  };

  enum class FitModel {
    kDiag,
    kFull,
  };

  enum class ProgressMode {
    kAuto,
    kEnabled,
    kDisabled,
  };

  // Rebin modes describe how one output selection is formed from normal sparse-axis bins.
  enum class RebinMode {
    kNative,
    kFactor,
    kRanges,
    kLegacyRanges,
  };

  enum class CoulombMode {
    kNone,
    kGamow,
    kFiniteSource,
  };

  enum class FiniteSourceMode {
    kFixed1D,
    kIterative1D,
  };

  inline bool UsesCoulomb(const CoulombMode mode) {
    return mode != CoulombMode::kNone;
  }

  inline int CoulombModeCode(const CoulombMode mode) {
    switch (mode) {
      case CoulombMode::kNone:
        return 0;
      case CoulombMode::kGamow:
        return 1;
      case CoulombMode::kFiniteSource:
        return 2;
    }
    return 0;
  }

  struct RangeBin {
    double min = 0.0;
    double max = 0.0;
    std::string label;

    [[nodiscard]] bool Contains(const double value) const {
      return value >= min && value < max;
    }
  };

  // Axis rebin settings are resolved against the input ROOT axis at build time.
  struct AxisRebinConfig {
    bool configured = false;
    bool enabled = false;
    RebinMode mode = RebinMode::kNative;
    std::optional<int> factor;
    std::optional<double> min;
    std::optional<double> max;
  };

  struct InputConfig {
    std::string input_root;
    std::string task_name;
    std::string same_event_subtask;
    std::string mixed_event_subtask;
    std::string sparse_object_name;
  };

  struct OutputConfig {
    std::string output_directory;
    std::string cf_root_name = "cf_output.root";
    std::string fit_root_name = "fit_output.root";
    std::string fit_summary_name = "fit_summary.tsv";
    std::string fit_report_directory;
    std::string fit_report_root_name = "fit_report.root";
    // Profile diagnostics are deliberately isolated from the detailed-fit and report files.
    std::string profile_root_name = "profile_likelihood.root";
    LogLevel log_level = LogLevel::kInfo;
  };

  struct BuildCfConfig {
    bool map_pair_phi_to_symmetric_range = false;
    bool write_normalized_se_me_1d_projections = false;
    bool reopen_output_file_per_slice = true;
    bool split_mixed_event_by_phi = false;
    bool split_same_event_by_qn = false;
    bool split_mixed_event_by_qn = false;
    AxisRebinConfig phi_rebin;
    AxisRebinConfig mt_rebin;
    ProgressMode progress = ProgressMode::kAuto;
  };

  // A parameter override changes only explicitly configured fit seeds, bounds, or supported fixed values.
  struct LevyFitParameterOverride {
    std::optional<double> initial;
    std::optional<double> min;
    std::optional<double> max;
    std::optional<double> fixed_value;
  };

  // The typed collection keeps chi-square and PML fit setup on one validated parameter contract.
  struct LevyFitParameterOverrides {
    LevyFitParameterOverride norm;
    LevyFitParameterOverride lambda;
    LevyFitParameterOverride rout2;
    LevyFitParameterOverride rside2;
    LevyFitParameterOverride rlong2;
    LevyFitParameterOverride routside2;
    LevyFitParameterOverride routlong2;
    LevyFitParameterOverride rsidelong2;
    LevyFitParameterOverride alpha;
    LevyFitParameterOverride baseline_q2;
  };

  struct LevyFitOptions {
    CoulombMode coulomb_mode = CoulombMode::kNone;
    FiniteSourceMode finite_source_mode = FiniteSourceMode::kFixed1D;
    bool use_core_halo_lambda = true;
    bool use_q2_baseline = false;
    bool use_pml = false;
    double fit_q_max = 0.15;
    LevyFitParameterOverrides parameters;
  };

  enum class ProfileSliceScope {
    kListed,
    kFitSelection,
  };

  enum class ProfileRetryStrategy {
    kReferenceAndBidirectionalNeighbors,
    kReferenceOnly,
  };

  // Legacy TMinuit uses process-local callback state; safe parallel execution
  // therefore partitions complete Coulomb groups across independent processes.
  enum class ProfileExecutionMode { kAlongsideFit, kProfileOnly };
  enum class ProfileParallelBackend { kSerial, kProcess, kThread };
  enum class ProfileMinimizerBackend { kLegacyTMinuit, kMinuit2 };
  enum class ProfileHesseStrategy { kAllAttempts, kNone };

  struct ProfileCheckpointConfig {
    bool enabled = false;
    bool resume = false;
    std::string run_id;
    std::string directory;
  };

  // One rectangular 1D/2D profile grid. Coordinates use native fit-parameter units.
  struct ProfileScanConfig {
    std::string id;
    std::vector<std::string> parameters;
    std::vector<int> points;
    std::vector<double> min;
    std::vector<double> max;
    bool refine = false;
    std::vector<int> refinement_points;
  };

  // This is opt-in because profiles are intentionally much more expensive than a fit.
  struct ProfileLikelihoodConfig {
    bool enabled = false;
    ProfileSliceScope slice_scope = ProfileSliceScope::kListed;
    std::vector<std::string> slice_ids;
    ProfileRetryStrategy retry_strategy = ProfileRetryStrategy::kReferenceAndBidirectionalNeighbors;
    bool write_likelihood_slice = true;
    ProfileExecutionMode execution_mode = ProfileExecutionMode::kAlongsideFit;
    ProfileParallelBackend parallel_backend = ProfileParallelBackend::kSerial;
    int workers = 1;
    ProfileMinimizerBackend minimizer_backend = ProfileMinimizerBackend::kLegacyTMinuit;
    ProfileHesseStrategy hesse_strategy = ProfileHesseStrategy::kAllAttempts;
    ProfileCheckpointConfig checkpoint;
    std::vector<double> contour_levels{1.0, 2.0, 4.0};
    std::vector<ProfileScanConfig> scans;
  };

  struct FitConfig {
    FitModel model = FitModel::kFull;
    LevyFitOptions options;
    std::optional<bool> map_pair_phi_to_symmetric_range;
    bool reopen_output_file_per_slice = true;
    ProgressMode progress = ProgressMode::kAuto;
    ProfileLikelihoodConfig profile_likelihood;
  };

  struct ApplicationConfig {
    InputConfig input;
    OutputConfig output;
    BuildCfConfig build;
    FitConfig fit;
    std::vector<RangeBin> centrality_bins;
    std::vector<RangeBin> mt_bins;
    std::vector<RangeBin> phi_bins;
    std::vector<RangeBin> qn_bins;
    std::vector<RangeBin> fit_centrality_bins;
    std::vector<RangeBin> fit_mt_bins;
  };

  struct SliceCatalogEntry {
    std::string slice_id;
    std::string group_id;
    std::string slice_directory;
    std::string se_object_path;
    std::string me_object_path;
    std::string cf_object_path;
    std::string projection_x_path;
    std::string projection_y_path;
    std::string projection_z_path;
    int centrality_index = -1;
    int mt_index = -1;
    int qn_index = -1;
    int phi_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    double qn_low = std::numeric_limits<double>::quiet_NaN();
    double qn_high = std::numeric_limits<double>::quiet_NaN();
    std::string qn_label = "qn_all";
    double raw_phi_low = 0.0;
    double raw_phi_high = 0.0;
    double raw_phi_center = 0.0;
    double display_phi_low = 0.0;
    double display_phi_high = 0.0;
    double display_phi_center = 0.0;
    bool build_uses_symmetric_phi_range = false;
    bool split_mixed_event_by_phi = false;
    bool split_mixed_event_by_qn = false;
    bool mt_rebin_enabled = false;
    std::string mt_rebin_mode = "legacy";
    bool phi_rebin_enabled = false;
    std::string phi_rebin_mode = "native";
    bool is_qn_integrated = true;
    bool is_phi_integrated = false;
  };

  struct LevyFitResult {
    std::string fit_model = "diag";
    std::string slice_id;
    std::string group_id;
    std::string slice_directory;
    int centrality_index = -1;
    int mt_index = -1;
    int qn_index = -1;
    int phi_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    double qn_low = std::numeric_limits<double>::quiet_NaN();
    double qn_high = std::numeric_limits<double>::quiet_NaN();
    std::string qn_label = "qn_all";
    double phi = std::numeric_limits<double>::quiet_NaN();
    double raw_phi_low = 0.0;
    double raw_phi_high = 0.0;
    double display_phi_low = 0.0;
    double display_phi_high = 0.0;
    bool fit_uses_symmetric_phi_range = false;
    bool mt_rebin_enabled = false;
    std::string mt_rebin_mode = "legacy";
    bool phi_rebin_enabled = false;
    std::string phi_rebin_mode = "native";
    bool is_qn_integrated = true;
    bool is_phi_integrated = false;
    double norm = std::numeric_limits<double>::quiet_NaN();
    double norm_err = std::numeric_limits<double>::quiet_NaN();
    double lambda = std::numeric_limits<double>::quiet_NaN();
    double lambda_err = std::numeric_limits<double>::quiet_NaN();
    double rout2 = std::numeric_limits<double>::quiet_NaN();
    double rout2_err = std::numeric_limits<double>::quiet_NaN();
    double rside2 = std::numeric_limits<double>::quiet_NaN();
    double rside2_err = std::numeric_limits<double>::quiet_NaN();
    double rlong2 = std::numeric_limits<double>::quiet_NaN();
    double rlong2_err = std::numeric_limits<double>::quiet_NaN();
    double routside2 = std::numeric_limits<double>::quiet_NaN();
    double routside2_err = std::numeric_limits<double>::quiet_NaN();
    double routlong2 = std::numeric_limits<double>::quiet_NaN();
    double routlong2_err = std::numeric_limits<double>::quiet_NaN();
    double rsidelong2 = std::numeric_limits<double>::quiet_NaN();
    double rsidelong2_err = std::numeric_limits<double>::quiet_NaN();
    double alpha = std::numeric_limits<double>::quiet_NaN();
    double alpha_err = std::numeric_limits<double>::quiet_NaN();
    double baseline_q2 = std::numeric_limits<double>::quiet_NaN();
    double baseline_q2_err = std::numeric_limits<double>::quiet_NaN();
    double fit_statistic = std::numeric_limits<double>::quiet_NaN();
    double edm = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int status = -1;
    int minuit_istat = -1;
    bool has_off_diagonal = false;
    bool uses_coulomb = false;
    std::string coulomb_mode = "none";
    std::string finite_source_mode;
    double finite_source_radius_fm = std::numeric_limits<double>::quiet_NaN();
    bool uses_core_halo_lambda = true;
    bool uses_q2_baseline = false;
    bool uses_pml = false;
  };

  struct CoulombKernelCatalogEntry {
    std::string group_id;
    int centrality_index = -1;
    int mt_index = -1;
    int qn_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    double qn_low = std::numeric_limits<double>::quiet_NaN();
    double qn_high = std::numeric_limits<double>::quiet_NaN();
    std::string qn_label = "qn_all";
    bool is_qn_integrated = true;
    std::string finite_source_mode;
    double seed_radius_fm = std::numeric_limits<double>::quiet_NaN();
    double final_radius_fm = std::numeric_limits<double>::quiet_NaN();
    bool cats_enabled = false;
    double kstar_min_mev = std::numeric_limits<double>::quiet_NaN();
    double kstar_max_mev = std::numeric_limits<double>::quiet_NaN();
    int kstar_bin_count = 0;
  };

  struct BuildCfRunStatistics {
    std::size_t requested_groups = 0;
    std::size_t stored_slices = 0;
    std::size_t skipped_zero_mixed_event_groups = 0;
    std::size_t skipped_zero_mixed_event_slices = 0;
    std::size_t skipped_zero_same_event_slices = 0;
    std::size_t mt_input_bins = 0;
    std::size_t mt_output_bins = 0;
    std::size_t phi_input_bins = 0;
    std::size_t phi_output_bins = 0;
    bool mt_rebin_enabled = false;
    bool phi_rebin_enabled = false;
    RebinMode mt_rebin_mode = RebinMode::kNative;
    RebinMode phi_rebin_mode = RebinMode::kNative;
  };

  struct FitRunStatistics {
    std::size_t catalog_slices = 0;
    std::size_t selected_slices = 0;
    std::size_t fitted_slices = 0;
    std::size_t skipped_missing_objects = 0;
    std::size_t skipped_missing_raw_histograms = 0;
    std::size_t profile_selected_slices = 0;
    std::size_t profile_completed_slices = 0;
    std::size_t profile_valid_points = 0;
    std::size_t profile_failed_points = 0;
    std::size_t profile_estimated_attempts = 0;
    std::size_t profile_estimated_slice_evaluations = 0;
    std::size_t profile_estimated_groups = 0;
    std::size_t profile_estimated_coarse_points_per_slice = 0;
    std::size_t profile_estimated_refined_points_per_slice = 0;
    std::size_t profile_configured_workers = 1;
    std::size_t profile_effective_workers = 1;
    std::string profile_output_path;
    bool profile_estimate_only = false;
  };

  inline bool NearlyEqual(const double lhs, const double rhs, const double tolerance = 1.0e-6) {
    return std::abs(lhs - rhs) < tolerance;
  }

  inline bool MatchesRangeBin(const RangeBin &lhs, const RangeBin &rhs, const double tolerance = 1.0e-6) {
    return NearlyEqual(lhs.min, rhs.min, tolerance) && NearlyEqual(lhs.max, rhs.max, tolerance);
  }

  inline bool IsValidRangeBin(const RangeBin &bin) {
    return std::isfinite(bin.min) && std::isfinite(bin.max) && bin.max > bin.min;
  }

  class ConfigError : public std::runtime_error {
   public:
    using std::runtime_error::runtime_error;
  };

}  // namespace exp_femto_3d
