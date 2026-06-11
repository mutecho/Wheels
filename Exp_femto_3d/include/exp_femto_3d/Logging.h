#pragma once

#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <string>
#include <thread>

#include "exp_femto_3d/Types.h"

namespace exp_femto_3d {

  // Value-only render input keeps progress formatting testable without a live
  // terminal or heartbeat thread.
  struct ProgressRenderSnapshot {
    std::string label;
    std::size_t total_steps = 0U;
    std::size_t completed_steps = 0U;
    int activity_frame = 0;
    std::chrono::steady_clock::duration elapsed{};
  };

  [[nodiscard]] std::string FormatProgressLine(const ProgressRenderSnapshot &snapshot);

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
      std::size_t completed_steps = 0;
      int activity_frame = 0;
      int last_percent = -1;
      std::size_t last_rendered_width = 0U;
      std::chrono::steady_clock::time_point start_time{};
      bool enabled = false;
      bool drawn = false;
      bool line_closed = true;
    };

    void BeginProgress(const std::string &label, std::size_t total_steps, ProgressMode mode) const;
    void UpdateProgress(std::size_t completed_steps) const;
    void TickProgress() const;
    void FinishProgress() const;
    void AbortProgress() const;
    void CloseProgressLine() const;
    void CloseProgressLineLocked() const;
    void RenderProgressLineLocked(std::chrono::steady_clock::time_point now) const;
    [[nodiscard]] bool IsProgressEnabled() const;
    void Log(LogLevel level, const std::string &message) const;

    LogLevel threshold_;
    mutable std::mutex output_mutex_;
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
    void StartHeartbeatIfNeeded();
    void RunHeartbeat();
    void StopHeartbeat();

    const Logger &logger_;
    std::mutex heartbeat_mutex_;
    std::condition_variable heartbeat_wake_;
    std::thread heartbeat_thread_;
    bool finished_ = false;
    bool heartbeat_started_ = false;
    bool stop_requested_ = false;
  };

}  // namespace exp_femto_3d
