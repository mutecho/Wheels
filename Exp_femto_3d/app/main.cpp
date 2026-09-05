#include <exception>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "exp_femto_3d/Config.h"
#include "exp_femto_3d/Logging.h"
#include "exp_femto_3d/Workflow.h"

namespace exp_femto_3d {
  namespace {

    struct CliArgs {
      std::string command;
      std::string config_path;
      std::optional<FitModel> model_override;
      std::optional<std::string> input_cf_root_override;
      bool profile_estimate_only = false;
    };

    void PrintUsage() {
      std::cout << "Usage:\n"
                << "  exp_femto_3d build-cf --config <file.toml>\n"
                << "  exp_femto_3d fit --config <file.toml> [--model full|diag] "
                   "[--input-cf-root <path>] [--profile-estimate-only]\n";
    }

    CliArgs ParseCli(const int argc, char **argv) {
      if (argc < 2) {
        throw std::runtime_error("Missing subcommand.");
      }

      CliArgs args;
      args.command = argv[1];
      for (int index = 2; index < argc; ++index) {
        const std::string token = argv[index];
        if (token == "--config") {
          if (index + 1 >= argc) {
            throw std::runtime_error("Missing value after --config.");
          }
          args.config_path = argv[++index];
          continue;
        }
        if (token == "--model") {
          if (index + 1 >= argc) {
            throw std::runtime_error("Missing value after --model.");
          }
          args.model_override = ParseFitModel(argv[++index]);
          continue;
        }
        if (token == "--input-cf-root") {
          if (index + 1 >= argc) {
            throw std::runtime_error("Missing value after --input-cf-root.");
          }
          args.input_cf_root_override = std::string(argv[++index]);
          continue;
        }
        if (token == "--profile-estimate-only") {
          args.profile_estimate_only = true;
          continue;
        }
        throw std::runtime_error("Unknown argument: " + token);
      }

      if (args.config_path.empty()) {
        throw std::runtime_error("--config is required.");
      }
      if (args.command != "build-cf" && args.command != "fit") {
        throw std::runtime_error("Unknown subcommand: " + args.command);
      }
      if (args.command != "fit" && args.profile_estimate_only) {
        throw std::runtime_error("--profile-estimate-only is valid only for the fit subcommand.");
      }
      return args;
    }

  }  // namespace
}  // namespace exp_femto_3d

int main(int argc, char **argv) {
  using namespace exp_femto_3d;

  try {
    const CliArgs args = ParseCli(argc, argv);
    ApplicationConfig config = LoadApplicationConfig(args.config_path);
    if (const char *worker_slices = std::getenv("EXP_FEMTO_3D_PROFILE_WORKER_SLICES")) {
      config.fit.profile_likelihood.slice_ids.clear();
      std::string encoded(worker_slices);
      std::size_t begin = 0;
      while (begin <= encoded.size()) {
        const std::size_t separator = encoded.find(';', begin);
        const std::string slice = encoded.substr(begin, separator == std::string::npos
                                                            ? std::string::npos : separator - begin);
        if (!slice.empty()) config.fit.profile_likelihood.slice_ids.push_back(slice);
        if (separator == std::string::npos) break;
        begin = separator + 1U;
      }
      const char *worker_output = std::getenv("EXP_FEMTO_3D_PROFILE_WORKER_OUTPUT");
      if (worker_output == nullptr || config.fit.profile_likelihood.slice_ids.empty()) {
        throw std::runtime_error("Incomplete internal profile worker environment.");
      }
      config.output.profile_root_name = worker_output;
      config.fit.profile_likelihood.checkpoint.enabled = false;
      config.fit.profile_likelihood.checkpoint.resume = false;
    }
    const Logger logger(config.output.log_level);

    if (args.command == "build-cf") {
      const BuildCfRunStatistics statistics = RunBuildCf(config, logger);
      std::cout << "build-cf stored_slices=" << statistics.stored_slices
                << " mt_rebin=" << (statistics.mt_rebin_enabled ? "enabled" : "disabled")
                << " mt_rebin_mode=" << ToString(statistics.mt_rebin_mode)
                << " mt_bins=" << statistics.mt_input_bins << "->" << statistics.mt_output_bins
                << " phi_rebin=" << (statistics.phi_rebin_enabled ? "enabled" : "disabled")
                << " phi_rebin_mode=" << ToString(statistics.phi_rebin_mode)
                << " phi_bins=" << statistics.phi_input_bins << "->" << statistics.phi_output_bins
                << " skipped_zero_me_groups=" << statistics.skipped_zero_mixed_event_groups
                << " skipped_zero_mixed_event_slices=" << statistics.skipped_zero_mixed_event_slices
                << " skipped_zero_se_slices=" << statistics.skipped_zero_same_event_slices << "\n";
      return 0;
    }

    const FitRunStatistics statistics =
        RunFit(config, logger, args.model_override, args.input_cf_root_override, args.profile_estimate_only,
               args.config_path);
    if (statistics.profile_estimate_only) {
      std::cout << "profile-estimate slices=" << statistics.profile_selected_slices
                << " groups=" << statistics.profile_estimated_groups
                << " coarse_points_per_slice=" << statistics.profile_estimated_coarse_points_per_slice
                << " refined_points_per_slice_max=" << statistics.profile_estimated_refined_points_per_slice
                << " attempts_max=" << statistics.profile_estimated_attempts
                << " likelihood_slice_evaluations_max=" << statistics.profile_estimated_slice_evaluations
                << " workers=" << statistics.profile_effective_workers
                << "/" << statistics.profile_configured_workers
                << "\n";
      return 0;
    }
    if (config.fit.profile_likelihood.enabled
        && config.fit.profile_likelihood.execution_mode == ProfileExecutionMode::kProfileOnly) {
      std::cout << "profile completed_slices=" << statistics.profile_completed_slices
                << " selected_slices=" << statistics.profile_selected_slices
                << " valid_points=" << statistics.profile_valid_points
                << " failed_points=" << statistics.profile_failed_points
                << " groups=" << statistics.profile_estimated_groups
                << " workers=" << statistics.profile_effective_workers
                << "/" << statistics.profile_configured_workers
                << " output=" << statistics.profile_output_path << "\n";
      return 0;
    }
    std::cout << "fit fitted_slices=" << statistics.fitted_slices << " selected_slices=" << statistics.selected_slices
              << " missing_objects=" << statistics.skipped_missing_objects
              << " missing_raw_histograms=" << statistics.skipped_missing_raw_histograms << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "[error] " << error.what() << "\n\n";
    PrintUsage();
    return 1;
  }
}
