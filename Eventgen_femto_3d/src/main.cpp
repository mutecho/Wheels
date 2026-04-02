#include "femto3d/AnalysisConfig.h"
#include "femto3d/Workflow.h"

#include <exception>
#include <iostream>
#include <string>

int main(int /*argc*/, char** /*argv*/) {
  femto3d::AnalysisConfig config = femto3d::MakeDefaultAnalysisConfig();
  std::string input_root_path = "input.root";
  std::string output_root_path = "output.root";
  std::string tree_name = "events";

  // Future CLI version:
  // if (argc >= 2) {
  //   input_root_path = argv[1];
  // }
  // if (argc >= 3) {
  //   output_root_path = argv[2];
  // }
  // if (argc >= 4) {
  //   tree_name = argv[3];
  // }

  config.input.tree_name = tree_name;

  try {
    const femto3d::AnalysisStatistics statistics =
        femto3d::RunAnalysis(config, input_root_path, output_root_path);
    std::cout << "Analysis completed.\n"
              << "  Input file: " << input_root_path << "\n"
              << "  Output file: " << output_root_path << "\n"
              << "  Tree name: " << tree_name << "\n"
              << "  Events read: " << statistics.events_read << "\n"
              << "  Selected particles: " << statistics.selected_particles << "\n"
              << "  Events with valid event plane: "
              << statistics.events_with_valid_event_plane << "\n"
              << "  Rejected no event-plane candidates: "
              << statistics.events_rejected_no_event_plane_candidates << "\n"
              << "  Rejected small Q vector: "
              << statistics.events_rejected_small_q_vector << "\n"
              << "  Rejected missing input event plane: "
              << statistics.events_rejected_missing_input_event_plane << "\n"
              << "  Rejected insufficient femto particles: "
              << statistics.events_rejected_insufficient_femto_particles << "\n"
              << "  Candidate pairs: " << statistics.candidate_pairs << "\n"
              << "  Accepted pairs: " << statistics.accepted_pairs << "\n"
              << "  Rejected close pairs: " << statistics.rejected_close_pairs
              << "\n"
              << "  Rejected invalid pair kinematics: "
              << statistics.rejected_invalid_pair_kinematics << "\n";
  } catch (const std::exception& error) {
    std::cerr << "Analysis failed: " << error.what() << "\n";
    return 2;
  }

  return 0;
}
