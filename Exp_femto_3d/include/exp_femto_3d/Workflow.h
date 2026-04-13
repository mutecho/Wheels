#pragma once

#include <optional>
#include <string>
#include <vector>

#include "exp_femto_3d/Logging.h"
#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  [[nodiscard]] BuildCfRunStatistics RunBuildCf(const ApplicationConfig &config, const Logger &logger);

  [[nodiscard]] FitRunStatistics RunFit(const ApplicationConfig &config,
                                        const Logger &logger,
                                        std::optional<FitModel> override_model = std::nullopt,
                                        std::optional<std::string> input_cf_root_path = std::nullopt);

  [[nodiscard]] std::vector<SliceCatalogEntry> LoadSliceCatalog(const std::string &cf_root_path);

}  // namespace exp_femto_3d
