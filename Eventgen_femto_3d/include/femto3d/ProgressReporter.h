#pragma once

#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <iosfwd>
#include <mutex>
#include <string>
#include <thread>

#include "femto3d/Config.h"
#include "femto3d/Workflow.h"

namespace femto3d {

  // Value-only render input keeps progress formatting testable without a live
  // terminal or heartbeat thread.
  struct ProgressRenderSnapshot {
    std::size_t total_events = 0U;
    std::size_t completed_events = 0U;
    int activity_frame = 0;
    std::chrono::steady_clock::duration elapsed{};
  };

  [[nodiscard]] std::string FormatProgressLine(const ProgressRenderSnapshot &snapshot);

  // CLI progress reporter for the analysis event loop. The workflow only sees
  // the AnalysisProgressSink interface, so physics code stays independent of
  // terminal rendering and locking.
  class ProgressReporter final : public AnalysisProgressSink {
   public:
    explicit ProgressReporter(ProgressMode progress_mode);
    ~ProgressReporter() override;

    void SetTotalEvents(std::size_t total_events) override;
    void UpdateCompletedEvents(std::size_t completed_events) override;
    void Finish();

   private:
    void StartHeartbeatIfNeededLocked();
    void RunHeartbeat();
    void StopHeartbeat();
    void RenderLocked(std::chrono::steady_clock::time_point now);
    void CloseLineLocked();

    std::size_t total_events_ = 0U;
    std::size_t completed_events_ = 0U;
    int activity_frame_ = 0;
    int last_rendered_percent_ = -1;
    std::size_t last_rendered_width_ = 0U;
    std::ostream *output_ = nullptr;
    std::chrono::steady_clock::time_point start_time_{};
    std::mutex mutex_;
    std::condition_variable heartbeat_wake_;
    std::thread heartbeat_thread_;
    ProgressMode progress_mode_ = ProgressMode::kAuto;
    bool enabled_ = false;
    bool heartbeat_started_ = false;
    bool drawn_ = false;
    bool line_closed_ = false;
    bool stop_requested_ = false;
  };

}  // namespace femto3d
