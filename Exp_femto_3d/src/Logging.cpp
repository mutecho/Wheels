#include "exp_femto_3d/Logging.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "exp_femto_3d/Config.h"

namespace exp_femto_3d {

  namespace {

    constexpr int kProgressBarWidth = 50;
    constexpr std::chrono::seconds kHeartbeatInterval(1);
    constexpr std::array<char, 4> kActivityFrames = {'/', '-', '\\', '|'};

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

    // Convert completed slice counts into the integer percent shown in the terminal line.
    int ProgressPercent(const std::size_t completed_steps, const std::size_t total_steps) {
      if (total_steps == 0U) {
        return 0;
      }

      const std::size_t clamped_completed_steps = std::min(completed_steps, total_steps);
      return static_cast<int>((100ULL * clamped_completed_steps) / total_steps);
    }

    // Use the same moving-head bar style as the Eventgen CLI progress output.
    std::string FormatProgressBar(const int percent) {
      std::string bar(static_cast<std::size_t>(kProgressBarWidth), '-');
      if (percent >= 100) {
        std::fill(bar.begin(), bar.end(), '=');
        return bar;
      }

      const int head_index = std::min((percent * kProgressBarWidth) / 100, kProgressBarWidth - 1);
      for (int index = 0; index < head_index; ++index) {
        bar[static_cast<std::size_t>(index)] = '=';
      }
      bar[static_cast<std::size_t>(head_index)] = '>';
      return bar;
    }

    // Rotate a compact activity marker so long operations show the process is alive.
    char ActivityFrame(const int frame) {
      const int non_negative_frame = std::max(frame, 0);
      const std::size_t frame_index =
          static_cast<std::size_t>(non_negative_frame % static_cast<int>(kActivityFrames.size()));
      return kActivityFrames[frame_index];
    }

    // Keep positive estimates in fixed HH:MM:SS form for stable terminal width.
    std::string FormatPositiveEta(const std::int64_t remaining_seconds) {
      if (remaining_seconds <= 0) {
        return "00:00";
      }

      const std::int64_t hours = remaining_seconds / 3600;
      const std::int64_t minutes = (remaining_seconds % 3600) / 60;
      const std::int64_t seconds = remaining_seconds % 60;

      std::ostringstream output;
      output << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':'
             << std::setw(2) << seconds;
      return output.str();
    }

    // Estimate remaining wall time from the average duration of completed slices.
    std::string FormatEta(const ProgressRenderSnapshot &snapshot) {
      if (snapshot.total_steps == 0U) {
        return "--:--";
      }

      const std::size_t completed_steps = std::min(snapshot.completed_steps, snapshot.total_steps);
      if (completed_steps == 0U) {
        return "--:--";
      }
      if (completed_steps >= snapshot.total_steps) {
        return "00:00";
      }

      const std::int64_t elapsed_millis = std::max<std::int64_t>(
          1, std::chrono::duration_cast<std::chrono::milliseconds>(snapshot.elapsed).count());
      const std::int64_t remaining_steps = static_cast<std::int64_t>(snapshot.total_steps - completed_steps);
      const std::int64_t completed_steps64 = static_cast<std::int64_t>(completed_steps);
      const std::int64_t remaining_millis =
          (elapsed_millis * remaining_steps + completed_steps64 - 1) / completed_steps64;
      return FormatPositiveEta((remaining_millis + 999) / 1000);
    }

  }  // namespace

  // Render a complete single-line progress message from value-only state.
  std::string FormatProgressLine(const ProgressRenderSnapshot &snapshot) {
    const int percent = ProgressPercent(snapshot.completed_steps, snapshot.total_steps);

    std::ostringstream output;
    if (!snapshot.label.empty()) {
      output << snapshot.label << ' ';
    }
    output << '[' << FormatProgressBar(percent) << "] " << percent << "% | "
           << ActivityFrame(snapshot.activity_frame) << " | ETA " << FormatEta(snapshot);
    return output.str();
  }

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

  // Initialize one progress span and leave rendering to UpdateProgress or the heartbeat.
  void Logger::BeginProgress(const std::string &label, const std::size_t total_steps, const ProgressMode mode) const {
    AbortProgress();
    std::lock_guard<std::mutex> lock(output_mutex_);
    progress_state_.label = label;
    progress_state_.total_steps = total_steps;
    progress_state_.completed_steps = 0U;
    progress_state_.activity_frame = 0;
    progress_state_.last_percent = -1;
    progress_state_.last_rendered_width = 0U;
    progress_state_.start_time = std::chrono::steady_clock::now();
    progress_state_.enabled = total_steps > 0 && ShouldEnableProgress(mode);
    progress_state_.drawn = false;
    progress_state_.line_closed = true;
  }

  // Record completed work and redraw only when visible progress changes.
  void Logger::UpdateProgress(const std::size_t completed_steps) const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    if (!progress_state_.enabled) {
      return;
    }

    progress_state_.completed_steps = std::min(completed_steps, progress_state_.total_steps);
    const int percent = ProgressPercent(progress_state_.completed_steps, progress_state_.total_steps);
    if (percent == progress_state_.last_percent && progress_state_.drawn && !progress_state_.line_closed) {
      return;
    }

    RenderProgressLineLocked(std::chrono::steady_clock::now());
  }

  // Refresh the activity frame and ETA even when no new slice has completed.
  void Logger::TickProgress() const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    if (!progress_state_.enabled) {
      return;
    }

    ++progress_state_.activity_frame;
    RenderProgressLineLocked(std::chrono::steady_clock::now());
  }

  // Return to column zero and erase the full terminal row before redrawing.
  // Some PTY frontends use insert-style rendering for a bare carriage return,
  // which otherwise leaves successive heartbeat frames concatenated.
  void Logger::RenderProgressLineLocked(const std::chrono::steady_clock::time_point now) const {
    const ProgressRenderSnapshot snapshot{progress_state_.label,
                                          progress_state_.total_steps,
                                          progress_state_.completed_steps,
                                          progress_state_.activity_frame,
                                          now - progress_state_.start_time};
    const std::string line = FormatProgressLine(snapshot);

    std::cerr << "\r\033[2K" << line;
    if (line.size() < progress_state_.last_rendered_width) {
      std::cerr << std::string(progress_state_.last_rendered_width - line.size(), ' ');
    }
    std::cerr << std::flush;

    progress_state_.drawn = true;
    progress_state_.line_closed = false;
    progress_state_.last_percent = ProgressPercent(progress_state_.completed_steps, progress_state_.total_steps);
    progress_state_.last_rendered_width = line.size();
  }

  // Force the terminal state before closing the live progress line.
  void Logger::FinishProgress() const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    if (!progress_state_.enabled) {
      return;
    }

    progress_state_.completed_steps = progress_state_.total_steps;
    if (!progress_state_.drawn || progress_state_.line_closed || progress_state_.last_percent < 100) {
      RenderProgressLineLocked(std::chrono::steady_clock::now());
    }
    CloseProgressLineLocked();
    progress_state_.enabled = false;
  }

  // Close an unfinished progress span without pretending all work completed.
  void Logger::AbortProgress() const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    CloseProgressLineLocked();
    progress_state_.enabled = false;
  }

  void Logger::CloseProgressLine() const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    CloseProgressLineLocked();
  }

  void Logger::CloseProgressLineLocked() const {
    if (!progress_state_.enabled || !progress_state_.drawn || progress_state_.line_closed) {
      return;
    }

    std::cerr << '\n' << std::flush;
    progress_state_.line_closed = true;
  }

  bool Logger::IsProgressEnabled() const {
    std::lock_guard<std::mutex> lock(output_mutex_);
    return progress_state_.enabled;
  }

  void Logger::Log(const LogLevel level, const std::string &message) const {
    if (static_cast<int>(level) < static_cast<int>(threshold_)) {
      return;
    }

    std::lock_guard<std::mutex> lock(output_mutex_);
    CloseProgressLineLocked();
    std::ostream &stream = level == LogLevel::kWarn || level == LogLevel::kError ? std::cerr : std::cout;
    stream << "[" << ToString(level) << "] " << message << "\n";
  }

  ProgressReporter::ProgressReporter(const Logger &logger,
                                     const std::string &label,
                                     const std::size_t total_steps,
                                     const ProgressMode mode)
      : logger_(logger) {
    logger_.BeginProgress(label, total_steps, mode);
    StartHeartbeatIfNeeded();
  }

  ProgressReporter::~ProgressReporter() {
    StopHeartbeat();
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
    StopHeartbeat();
    logger_.FinishProgress();
    finished_ = true;
  }

  // Start a background heartbeat only for progress modes that actually render.
  void ProgressReporter::StartHeartbeatIfNeeded() {
    if (!logger_.IsProgressEnabled()) {
      return;
    }

    std::lock_guard<std::mutex> lock(heartbeat_mutex_);
    if (heartbeat_started_) {
      return;
    }

    stop_requested_ = false;
    heartbeat_started_ = true;
    heartbeat_thread_ = std::thread(&ProgressReporter::RunHeartbeat, this);
  }

  // Wake once per second so the activity frame and ETA stay live during long ROOT calls.
  void ProgressReporter::RunHeartbeat() {
    std::unique_lock<std::mutex> lock(heartbeat_mutex_);
    while (!stop_requested_) {
      heartbeat_wake_.wait_for(lock, kHeartbeatInterval, [this]() {
        return stop_requested_;
      });
      if (stop_requested_) {
        break;
      }

      lock.unlock();
      logger_.TickProgress();
      lock.lock();
    }
  }

  // Stop the heartbeat before finish/abort so no background thread touches the logger afterward.
  void ProgressReporter::StopHeartbeat() {
    {
      std::lock_guard<std::mutex> lock(heartbeat_mutex_);
      if (!heartbeat_started_) {
        return;
      }
      stop_requested_ = true;
    }
    heartbeat_wake_.notify_all();

    if (heartbeat_thread_.joinable()) {
      heartbeat_thread_.join();
    }

    std::lock_guard<std::mutex> lock(heartbeat_mutex_);
    heartbeat_started_ = false;
  }

}  // namespace exp_femto_3d
