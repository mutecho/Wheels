#pragma once

#include "femto3d/AnalysisConfig.h"
#include "femto3d/AnalysisTypes.h"

#include <vector>

namespace femto3d {

// Loads the configured input ROOT file and normalizes every event into
// the EventData boundary consumed by the analysis workflow.
[[nodiscard]] std::vector<EventData> LoadEventData(
    const ApplicationConfig& config);

}  // namespace femto3d
