# Work Items

## WI-001 强制化 project-state closeout gate

- 状态: open
- 背景:
  - 本轮代码与验证都已完成，但 `project-state` 写回没有在第一次收尾时自动执行。
  - 原因是父线程没有把 closeout gate 作为正常完成前的硬步骤。
- 建议:
  - 在 adopted `project-state/` 仓库中，把“写回 `project-state`”固定放入计划尾项。
  - 结束前显式检查 `state_changed_this_task` 与 `project_state_written_this_task`。

## WI-002 评估可委派的验证子代理边界

- 状态: open
- 背景:
  - 这轮测试设计适合 sidecar delegation，但最终 build/test 执行位于关键路径，主线程本地执行更高效。
- 建议:
  - 区分两类验证：
    - 关键路径 reproduce-fix-rerun：保留父线程执行
    - 已稳定 smoke/regression：可考虑委派给专门验证子代理
  - 在后续任务中把这个分界条件写进启动 plan 或 handoff 模板。
