#include "exp_femto_3d/Logging.h"

#include <iostream>

#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  Logger::Logger(const LogLevel threshold) : threshold_(threshold) {
  }

  void Logger::Debug(const std::string &message) const {
    Log(LogLevel::kDebug, message);
  }

  void Logger::Info(const std::string &message) const {
    Log(LogLevel::kInfo, message);
  }

  void Logger::Warn(const std::string &message) const {
    Log(LogLevel::kWarn, message);
  }

  void Logger::Error(const std::string &message) const {
    Log(LogLevel::kError, message);
  }

  void Logger::Log(const LogLevel level, const std::string &message) const {
    if (static_cast<int>(level) < static_cast<int>(threshold_)) {
      return;
    }

    std::ostream &stream = level == LogLevel::kWarn || level == LogLevel::kError ? std::cerr : std::cout;
    stream << "[" << ToString(level) << "] " << message << "\n";
  }

}  // namespace exp_femto_3d
