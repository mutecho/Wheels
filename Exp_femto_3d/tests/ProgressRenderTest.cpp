#include "exp_femto_3d/Logging.h"

#include <chrono>
#include <iostream>
#include <string>

namespace {

  bool Contains(const std::string &text, const std::string &expected) {
    return text.find(expected) != std::string::npos;
  }

}  // namespace

int main() {
  using namespace exp_femto_3d;

  // Half-complete snapshots should expose the same operator-facing fields as
  // the live progress line: stage label, percent, activity frame, and ETA.
  const ProgressRenderSnapshot half_done{"fit", 10U, 5U, 1, std::chrono::seconds(20)};
  const std::string half_done_line = FormatProgressLine(half_done);
  if (!Contains(half_done_line, "fit [") || !Contains(half_done_line, "50%")
      || !Contains(half_done_line, " | - | ETA 00:00:20")) {
    std::cerr << "Expected progress line to include label, percent, activity frame, and ETA.\n"
              << half_done_line << "\n";
    return 1;
  }

  // A not-yet-started snapshot cannot estimate remaining time, so it should
  // advertise an unknown ETA rather than pretending the remaining time is zero.
  const ProgressRenderSnapshot not_started{"build-cf", 10U, 0U, 0, std::chrono::seconds(0)};
  const std::string not_started_line = FormatProgressLine(not_started);
  if (!Contains(not_started_line, "0%") || !Contains(not_started_line, "ETA --:--")) {
    std::cerr << "Expected zero-complete progress line to show an unknown ETA.\n"
              << not_started_line << "\n";
    return 2;
  }

  // Completed snapshots should clamp both percent and ETA to the terminal
  // state even if the caller reports more completed steps than planned.
  const ProgressRenderSnapshot complete{"build-cf", 10U, 12U, 3, std::chrono::seconds(30)};
  const std::string complete_line = FormatProgressLine(complete);
  if (!Contains(complete_line, "100%") || !Contains(complete_line, "ETA 00:00")) {
    std::cerr << "Expected completed progress line to show 100% and zero ETA.\n"
              << complete_line << "\n";
    return 3;
  }

  return 0;
}
