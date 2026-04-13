#pragma once

#include <string>

#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  class Logger {
   public:
    explicit Logger(LogLevel threshold);

    [[nodiscard]] LogLevel threshold() const {
      return threshold_;
    }

    void Debug(const std::string &message) const;
    void Info(const std::string &message) const;
    void Warn(const std::string &message) const;
    void Error(const std::string &message) const;

   private:
    void Log(LogLevel level, const std::string &message) const;

    LogLevel threshold_;
  };

}  // namespace exp_femto_3d
