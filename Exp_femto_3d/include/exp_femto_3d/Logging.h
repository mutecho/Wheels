#pragma once

#include <cstddef>
#include <string>

#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  class ProgressReporter;

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
    friend class ProgressReporter;

    struct ProgressState {
      std::string label;
      std::size_t total_steps = 0;
      int last_percent = -1;
      bool enabled = false;
      bool drawn = false;
      bool line_closed = true;
    };

    void BeginProgress(const std::string &label, std::size_t total_steps, ProgressMode mode) const;
    void UpdateProgress(std::size_t completed_steps) const;
    void FinishProgress() const;
    void AbortProgress() const;
    void CloseProgressLine() const;
    void Log(LogLevel level, const std::string &message) const;

    LogLevel threshold_;
    mutable ProgressState progress_state_;
  };

  // Keep workflow progress reporting outside the histogram math so build and fit
  // only expose completed-slice counts to the CLI layer.
  class ProgressReporter {
   public:
    ProgressReporter(const Logger &logger, const std::string &label, std::size_t total_steps, ProgressMode mode);
    ~ProgressReporter();

    void Update(std::size_t completed_steps) const;
    void Finish();

   private:
    const Logger &logger_;
    bool finished_ = false;
  };

}  // namespace exp_femto_3d
