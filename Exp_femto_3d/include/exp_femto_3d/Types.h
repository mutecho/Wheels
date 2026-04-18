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

  struct RangeBin {
    double min = 0.0;
    double max = 0.0;
    std::string label;

    [[nodiscard]] bool Contains(const double value) const {
      return value >= min && value < max;
    }
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
    LogLevel log_level = LogLevel::kInfo;
  };

  struct BuildCfConfig {
    bool map_pair_phi_to_symmetric_range = false;
    bool write_normalized_se_me_1d_projections = false;
    bool reopen_output_file_per_slice = true;
    ProgressMode progress = ProgressMode::kAuto;
  };

  struct LevyFitOptions {
    bool use_coulomb = false;
    bool use_core_halo_lambda = true;
    bool use_q2_baseline = false;
    bool use_pml = false;
    double fit_q_max = 0.15;
  };

  struct FitConfig {
    FitModel model = FitModel::kFull;
    LevyFitOptions options;
    std::optional<bool> map_pair_phi_to_symmetric_range;
    bool reopen_output_file_per_slice = true;
    ProgressMode progress = ProgressMode::kAuto;
  };

  struct ApplicationConfig {
    InputConfig input;
    OutputConfig output;
    BuildCfConfig build;
    FitConfig fit;
    std::vector<RangeBin> centrality_bins;
    std::vector<RangeBin> mt_bins;
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
    int phi_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    double raw_phi_low = 0.0;
    double raw_phi_high = 0.0;
    double raw_phi_center = 0.0;
    double display_phi_low = 0.0;
    double display_phi_high = 0.0;
    double display_phi_center = 0.0;
    bool build_uses_symmetric_phi_range = false;
    bool is_phi_integrated = false;
  };

  struct LevyFitResult {
    std::string fit_model = "diag";
    std::string slice_id;
    std::string group_id;
    std::string slice_directory;
    int centrality_index = -1;
    int mt_index = -1;
    int phi_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    double phi = std::numeric_limits<double>::quiet_NaN();
    bool fit_uses_symmetric_phi_range = false;
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
    bool uses_core_halo_lambda = true;
    bool uses_q2_baseline = false;
    bool uses_pml = false;
  };

  struct BuildCfRunStatistics {
    std::size_t requested_groups = 0;
    std::size_t stored_slices = 0;
    std::size_t skipped_zero_mixed_event_groups = 0;
    std::size_t skipped_zero_same_event_slices = 0;
  };

  struct FitRunStatistics {
    std::size_t catalog_slices = 0;
    std::size_t selected_slices = 0;
    std::size_t fitted_slices = 0;
    std::size_t skipped_missing_objects = 0;
    std::size_t skipped_missing_raw_histograms = 0;
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
