#include "exp_femto_3d/Logging.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <unistd.h>

#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  namespace {

    constexpr int kProgressBarWidth = 20;

    bool ShouldEnableProgress(const ProgressMode mode) {
      switch (mode) {
        case ProgressMode::kAuto:
          return ::isatty(STDERR_FILENO) != 0;
        case ProgressMode::kEnabled:
          return true;
        case ProgressMode::kDisabled:
          return false;
      }
      return false;
    }

  }  // namespace

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

  void Logger::BeginProgress(const std::string &label, const std::size_t total_steps, const ProgressMode mode) const {
    AbortProgress();
    progress_state_.label = label;
    progress_state_.total_steps = total_steps;
    progress_state_.last_percent = -1;
    progress_state_.enabled = total_steps > 0 && ShouldEnableProgress(mode);
    progress_state_.drawn = false;
    progress_state_.line_closed = true;
  }

  void Logger::UpdateProgress(const std::size_t completed_steps) const {
    if (!progress_state_.enabled) {
      return;
    }

    const std::size_t clamped_completed_steps = std::min(completed_steps, progress_state_.total_steps);
    const int percent = progress_state_.total_steps > 0
                            ? static_cast<int>((100 * clamped_completed_steps) / progress_state_.total_steps)
                            : 0;
    if (percent == progress_state_.last_percent && progress_state_.drawn && !progress_state_.line_closed) {
      return;
    }

    std::string bar(static_cast<std::size_t>(kProgressBarWidth), '-');
    if (percent >= 100) {
      std::fill(bar.begin(), bar.end(), '=');
    } else {
      const int head_index = std::min((percent * kProgressBarWidth) / 100, kProgressBarWidth - 1);
      for (int index = 0; index < head_index; ++index) {
        bar[static_cast<std::size_t>(index)] = '=';
      }
      bar[static_cast<std::size_t>(head_index)] = '>';
    }

    std::cerr << '\r';
    if (!progress_state_.label.empty()) {
      std::cerr << progress_state_.label << ' ';
    }
    std::cerr << '[' << bar << "] " << percent << '%' << std::flush;
    progress_state_.drawn = true;
    progress_state_.line_closed = false;
    progress_state_.last_percent = percent;
  }

  void Logger::FinishProgress() const {
    if (!progress_state_.enabled) {
      return;
    }

    UpdateProgress(progress_state_.total_steps);
    AbortProgress();
  }

  void Logger::AbortProgress() const {
    CloseProgressLine();
    progress_state_.enabled = false;
  }

  void Logger::CloseProgressLine() const {
    if (!progress_state_.enabled || !progress_state_.drawn || progress_state_.line_closed) {
      return;
    }

    std::cerr << '\n' << std::flush;
    progress_state_.line_closed = true;
  }

  void Logger::Log(const LogLevel level, const std::string &message) const {
    if (static_cast<int>(level) < static_cast<int>(threshold_)) {
      return;
    }

    CloseProgressLine();
    std::ostream &stream = level == LogLevel::kWarn || level == LogLevel::kError ? std::cerr : std::cout;
    stream << "[" << ToString(level) << "] " << message << "\n";
  }

  ProgressReporter::ProgressReporter(const Logger &logger,
                                     const std::string &label,
                                     const std::size_t total_steps,
                                     const ProgressMode mode)
      : logger_(logger) {
    logger_.BeginProgress(label, total_steps, mode);
  }

  ProgressReporter::~ProgressReporter() {
    if (!finished_) {
      logger_.AbortProgress();
    }
  }

  void ProgressReporter::Update(const std::size_t completed_steps) const {
    logger_.UpdateProgress(completed_steps);
  }

  void ProgressReporter::Finish() {
    if (finished_) {
      return;
    }
    logger_.FinishProgress();
    finished_ = true;
  }

}  // namespace exp_femto_3d
