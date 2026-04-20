#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "exp_femto_1d/Config.h"
#include "exp_femto_1d/Logging.h"
#include "exp_femto_1d/Workflow.h"

namespace exp_femto_1d {
  namespace {

    struct CliArgs {
      std::string command;
      std::string config_path;
      std::optional<std::string> input_cf_root_override;
    };

    void PrintUsage() {
      std::cout << "Usage:\n"
                << "  exp_femto_1d build-cf --config <file.toml>\n"
                << "  exp_femto_1d fit --config <file.toml> [--input-cf-root <path>]\n";
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
}  // namespace exp_femto_1d

int main(int argc, char **argv) {
  using namespace exp_femto_1d;

  try {
    const CliArgs args = ParseCli(argc, argv);
    const ApplicationConfig config = LoadApplicationConfig(args.config_path);
    const Logger logger(config.output.log_level);

    if (args.command == "build-cf") {
      const BuildCfRunStatistics statistics = RunBuildCf(config, logger);
      std::cout << "build-cf"
                << " stored_slices=" << statistics.stored_slices
                << " skipped_zero_mixed_event_groups=" << statistics.skipped_zero_mixed_event_groups
                << " skipped_zero_same_event_slices=" << statistics.skipped_zero_same_event_slices << "\n";
      return 0;
    }

    const FitRunStatistics statistics = RunFit(config, logger, args.input_cf_root_override);
    std::cout << "fit"
              << " catalog_slices=" << statistics.catalog_slices
              << " selected_slices=" << statistics.selected_slices
              << " fitted_slices=" << statistics.fitted_slices
              << " skipped_missing_objects=" << statistics.skipped_missing_objects
              << " skipped_failed_fits=" << statistics.skipped_failed_fits << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "[error] " << error.what() << "\n\n";
    PrintUsage();
    return 1;
  }
}
