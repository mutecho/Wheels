#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace exp_femto_1d {

  constexpr double kPi = 3.14159265358979323846;

  enum class LogLevel {
    kDebug = 0,
    kInfo = 1,
    kWarn = 2,
    kError = 3,
  };

  enum class ProgressMode {
    kAuto,
    kEnabled,
    kDisabled,
  };

  enum class RegionKind {
    kMinBias = 0,
    kInPlane = 1,
    kOutOfPlane = 2,
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
    double norm_low = 0.5;
    double norm_high = 0.8;
    double kstar_min = 0.0;
    double kstar_max = 0.8;
    bool reopen_output_file_per_slice = true;
    ProgressMode progress = ProgressMode::kAuto;
  };

  struct FitConfig {
    double fit_kstar_max = 0.25;
    bool use_coulomb = true;
    bool reopen_output_file_per_slice = true;
    ProgressMode progress = ProgressMode::kAuto;

    double baseline_p0_init = 1.0;
    double baseline_p0_min = 0.9;
    double baseline_p0_max = 1.1;

    double baseline_p1_init = 0.0;
    double baseline_p1_min = -0.2;
    double baseline_p1_max = 0.2;

    double baseline_p2_init = 1.0e-5;
    double baseline_p2_min = -5.0;
    double baseline_p2_max = 5.0;

    bool baseline_p3_fixed = true;
    double baseline_p3_value = 0.0;

    bool baseline_p4_fixed = true;
    double baseline_p4_value = 0.0;

    double source_size_init = 6.0;
    double source_size_min = 3.0;
    double source_size_max = 16.0;

    unsigned cats_num_mom_bins = 250;
    double cats_kmin_mev = 0.0;
    double cats_kmax_mev = 250.0;
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
    int centrality_index = -1;
    int mt_index = -1;
    int region_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    std::string region_name;
    RegionKind region_kind = RegionKind::kMinBias;
    double ep_low_1 = 0.0;
    double ep_high_1 = 0.0;
    double ep_low_2 = 0.0;
    double ep_high_2 = 0.0;
    bool has_second_interval = false;
    double norm_low = 0.0;
    double norm_high = 0.0;
    double kstar_min = 0.0;
    double kstar_max = 0.0;
  };

  struct PiPiFitResult {
    std::string slice_id;
    std::string group_id;
    std::string slice_directory;
    int centrality_index = -1;
    int mt_index = -1;
    int region_index = -1;
    double cent_low = 0.0;
    double cent_high = 0.0;
    double mt_low = 0.0;
    double mt_high = 0.0;
    std::string region_name;
    double fit_kstar_max = std::numeric_limits<double>::quiet_NaN();
    bool uses_coulomb = true;
    double baseline_p0 = std::numeric_limits<double>::quiet_NaN();
    double baseline_p0_err = std::numeric_limits<double>::quiet_NaN();
    double baseline_p1 = std::numeric_limits<double>::quiet_NaN();
    double baseline_p1_err = std::numeric_limits<double>::quiet_NaN();
    double baseline_p2 = std::numeric_limits<double>::quiet_NaN();
    double baseline_p2_err = std::numeric_limits<double>::quiet_NaN();
    double baseline_p3 = std::numeric_limits<double>::quiet_NaN();
    double baseline_p3_err = std::numeric_limits<double>::quiet_NaN();
    double baseline_p4 = std::numeric_limits<double>::quiet_NaN();
    double baseline_p4_err = std::numeric_limits<double>::quiet_NaN();
    double source_size = std::numeric_limits<double>::quiet_NaN();
    double source_size_err = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    double fit_statistic = std::numeric_limits<double>::quiet_NaN();
    double edm = std::numeric_limits<double>::quiet_NaN();
    int status = -1;
    int minuit_istat = -1;
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
    std::size_t skipped_failed_fits = 0;
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

}  // namespace exp_femto_1d
