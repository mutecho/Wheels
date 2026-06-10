#include "femto3d/ProgressReporter.h"

#include <unistd.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {

  constexpr int kBarWidth = 50;
  constexpr std::chrono::seconds kHeartbeatInterval(1);
  constexpr std::array<char, 4> kActivityFrames = {'/', '-', '\\', '|'};

  bool ShouldEnableProgress(const femto3d::ProgressMode progress_mode) {
    switch (progress_mode) {
      case femto3d::ProgressMode::kAuto:
        return ::isatty(STDERR_FILENO) != 0;
      case femto3d::ProgressMode::kEnabled:
        return true;
      case femto3d::ProgressMode::kDisabled:
        return false;
    }

    return false;
  }

  int ProgressPercent(const std::size_t completed_events, const std::size_t total_events) {
    if (total_events == 0U) {
      return 0;
    }

    const std::size_t clamped_completed_events = std::min(completed_events, total_events);
    return static_cast<int>((100ULL * clamped_completed_events) / total_events);
  }

  std::string FormatProgressBar(const int percent) {
    std::string bar(static_cast<std::size_t>(kBarWidth), '-');
    if (percent >= 100) {
      std::fill(bar.begin(), bar.end(), '=');
      return bar;
    }

    const int head_index = std::min((percent * kBarWidth) / 100, kBarWidth - 1);
    for (int index = 0; index < head_index; ++index) {
      bar[static_cast<std::size_t>(index)] = '=';
    }
    bar[static_cast<std::size_t>(head_index)] = '>';
    return bar;
  }

  char ActivityFrame(const int frame) {
    const int non_negative_frame = std::max(frame, 0);
    const std::size_t frame_index =
        static_cast<std::size_t>(non_negative_frame % static_cast<int>(kActivityFrames.size()));
    return kActivityFrames[frame_index];
  }

  std::string FormatPositiveEta(const std::int64_t remaining_seconds) {
    if (remaining_seconds <= 0) {
      return "00:00";
    }

    const std::int64_t hours = remaining_seconds / 3600;
    const std::int64_t minutes = (remaining_seconds % 3600) / 60;
    const std::int64_t seconds = remaining_seconds % 60;

    std::ostringstream output;
    output << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2)
           << seconds;
    return output.str();
  }

  std::string FormatEta(const femto3d::ProgressRenderSnapshot &snapshot) {
    if (snapshot.total_events == 0U) {
      return "--:--";
    }

    const std::size_t completed_events = std::min(snapshot.completed_events, snapshot.total_events);
    if (completed_events == 0U) {
      return "--:--";
    }
    if (completed_events >= snapshot.total_events) {
      return "00:00";
    }

    const std::int64_t elapsed_millis = std::max<std::int64_t>(
        1, std::chrono::duration_cast<std::chrono::milliseconds>(snapshot.elapsed).count());
    const std::int64_t remaining_events = static_cast<std::int64_t>(snapshot.total_events - completed_events);
    const std::int64_t completed_events64 = static_cast<std::int64_t>(completed_events);
    const std::int64_t remaining_millis =
        (elapsed_millis * remaining_events + completed_events64 - 1) / completed_events64;
    return FormatPositiveEta((remaining_millis + 999) / 1000);
  }

}  // namespace

namespace femto3d {

  std::string FormatProgressLine(const ProgressRenderSnapshot &snapshot) {
    const int percent = ProgressPercent(snapshot.completed_events, snapshot.total_events);
    std::ostringstream output;
    output << '[' << FormatProgressBar(percent) << "] " << percent << "% | "
           << ActivityFrame(snapshot.activity_frame) << " | ETA " << FormatEta(snapshot);
    return output.str();
  }

  ProgressReporter::ProgressReporter(const ProgressMode progress_mode)
      : output_(&std::cerr),
        start_time_(std::chrono::steady_clock::now()),
        progress_mode_(progress_mode),
        enabled_(ShouldEnableProgress(progress_mode)) {}

  ProgressReporter::~ProgressReporter() {
    StopHeartbeat();

    std::lock_guard<std::mutex> lock(mutex_);
    CloseLineLocked();
  }

  void ProgressReporter::SetTotalEvents(const std::size_t total_events) {
    if (!enabled_) {
      return;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    total_events_ = total_events;
    completed_events_ = std::min(completed_events_, total_events_);
    StartHeartbeatIfNeededLocked();
  }

  void ProgressReporter::UpdateCompletedEvents(const std::size_t completed_events) {
    if (!enabled_ || total_events_ == 0U) {
      return;
    }

    std::lock_guard<std::mutex> lock(mutex_);
    completed_events_ = std::min(completed_events, total_events_);

    const int percent = ProgressPercent(completed_events_, total_events_);
    if (drawn_ && percent == last_rendered_percent_) {
      return;
    }

    RenderLocked(std::chrono::steady_clock::now());
  }

  void ProgressReporter::Finish() {
    if (!enabled_) {
      return;
    }

    StopHeartbeat();

    std::lock_guard<std::mutex> lock(mutex_);
    completed_events_ = total_events_;
    RenderLocked(std::chrono::steady_clock::now());
    CloseLineLocked();
  }

  void ProgressReporter::StartHeartbeatIfNeededLocked() {
    if (!enabled_ || heartbeat_started_ || total_events_ == 0U) {
      return;
    }

    heartbeat_started_ = true;
    heartbeat_thread_ = std::thread(&ProgressReporter::RunHeartbeat, this);
  }

  void ProgressReporter::RunHeartbeat() {
    std::unique_lock<std::mutex> lock(mutex_);
    while (!stop_requested_) {
      heartbeat_wake_.wait_for(lock, kHeartbeatInterval, [this]() {
        return stop_requested_;
      });
      if (stop_requested_) {
        break;
      }

      ++activity_frame_;
      RenderLocked(std::chrono::steady_clock::now());
    }
  }

  void ProgressReporter::StopHeartbeat() {
    if (!enabled_ || !heartbeat_started_) {
      return;
    }

    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_requested_ = true;
    }
    heartbeat_wake_.notify_all();

    if (heartbeat_thread_.joinable()) {
      heartbeat_thread_.join();
    }
  }

  void ProgressReporter::RenderLocked(const std::chrono::steady_clock::time_point now) {
    const ProgressRenderSnapshot snapshot{total_events_, completed_events_, activity_frame_, now - start_time_};
    const std::string line = FormatProgressLine(snapshot);

    (*output_) << '\r' << line;
    if (line.size() < last_rendered_width_) {
      (*output_) << std::string(last_rendered_width_ - line.size(), ' ');
    }
    (*output_) << std::flush;

    drawn_ = true;
    line_closed_ = false;
    last_rendered_percent_ = ProgressPercent(completed_events_, total_events_);
    last_rendered_width_ = line.size();
  }

  void ProgressReporter::CloseLineLocked() {
    if (drawn_ && !line_closed_) {
      (*output_) << '\n' << std::flush;
      line_closed_ = true;
    }
  }

}  // namespace femto3d
