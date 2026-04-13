#include <exception>
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
    };

    void PrintUsage() {
      std::cout << "Usage:\n"
                << "  exp_femto_3d build-cf --config <file.toml>\n"
                << "  exp_femto_3d fit --config <file.toml> [--model full|diag] "
                   "[--input-cf-root <path>]\n";
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
        throw std::runtime_error("Unknown argument: " + token);
      }

      if (args.config_path.empty()) {
        throw std::runtime_error("--config is required.");
      }
      if (args.command != "build-cf" && args.command != "fit") {
        throw std::runtime_error("Unknown subcommand: " + args.command);
      }
      return args;
    }

  }  // namespace
}  // namespace exp_femto_3d

int main(int argc, char **argv) {
  using namespace exp_femto_3d;

  try {
    const CliArgs args = ParseCli(argc, argv);
    const ApplicationConfig config = LoadApplicationConfig(args.config_path);
    const Logger logger(config.output.log_level);

    if (args.command == "build-cf") {
      const BuildCfRunStatistics statistics = RunBuildCf(config, logger);
      std::cout << "build-cf stored_slices=" << statistics.stored_slices
                << " skipped_zero_me_groups=" << statistics.skipped_zero_mixed_event_groups
                << " skipped_zero_se_slices=" << statistics.skipped_zero_same_event_slices << "\n";
      return 0;
    }

    const FitRunStatistics statistics = RunFit(config, logger, args.model_override, args.input_cf_root_override);
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
