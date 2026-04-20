#pragma once

#include <optional>
#include <string>
#include <vector>

#include "exp_femto_1d/Logging.h"
#include "exp_femto_1d/Types.h"

namespace exp_femto_1d {

  [[nodiscard]] BuildCfRunStatistics RunBuildCf(const ApplicationConfig &config, const Logger &logger);

  [[nodiscard]] FitRunStatistics RunFit(const ApplicationConfig &config,
                                        const Logger &logger,
                                        std::optional<std::string> input_cf_root_path = std::nullopt);

  [[nodiscard]] std::vector<SliceCatalogEntry> LoadSliceCatalog(const std::string &cf_root_path);

}  // namespace exp_femto_1d
