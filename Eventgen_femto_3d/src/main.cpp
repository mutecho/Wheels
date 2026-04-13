#include <exception>
#include <iostream>

#include "femto3d/Config.h"
#include "femto3d/Workflow.h"

int main(int argc, char **argv) {
  using namespace femto3d;

  try {
    const CliOptions cli_options = ParseCliArgs(argc, argv);
    ApplicationConfig config = LoadApplicationConfig(cli_options.config_path);
    ApplyCliOverrides(cli_options, config);

    const AnalysisStatistics statistics = RunAnalysis(config);
    std::cout << "Analysis completed.\n"
              << "  Config file: " << cli_options.config_path << "\n"
              << "  Input file: " << config.input_root_path << "\n"
              << "  Output file: " << config.output_root_path << "\n"
              << "  Input schema: " << ToString(config.input_schema) << "\n"
              << "  Events read: " << statistics.events_read << "\n"
              << "  Selected particles: " << statistics.selected_particles << "\n"
              << "  Events with valid event plane: " << statistics.events_with_valid_event_plane << "\n"
              << "  Rejected no event-plane candidates: " << statistics.events_rejected_no_event_plane_candidates
              << "\n"
              << "  Rejected small Q vector: " << statistics.events_rejected_small_q_vector << "\n"
              << "  Rejected missing input event plane: " << statistics.events_rejected_missing_input_event_plane
              << "\n"
              << "  Rejected insufficient femto particles: " << statistics.events_rejected_insufficient_femto_particles
              << "\n"
              << "  Candidate pairs: " << statistics.candidate_pairs << "\n"
              << "  Accepted pairs: " << statistics.accepted_pairs << "\n"
              << "  Rejected close pairs: " << statistics.rejected_close_pairs << "\n"
              << "  Rejected invalid pair kinematics: " << statistics.rejected_invalid_pair_kinematics << "\n"
              << "  R2 summary points skipped due to invalid HBT errors: "
              << statistics.r2_summary_points_skipped_invalid_hbt_error << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "[error] " << error.what() << "\n\n";
    PrintUsage(std::cerr);
    return 1;
  }
}
