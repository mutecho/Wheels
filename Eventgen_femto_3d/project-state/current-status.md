# Current Status

## 2026-04-13 project-state 对齐更新与 wrapper 路径核对

- 状态: 已完成本轮台账同步与最新证据补记
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 对照当前 `HEAD` 校对 `project-state` 与代码/构建输出是否一致。
- 修正文档中关于测试数量的旧描述：从早先的 `6/6` 更新为当前构建注册的 `8` 项测试。
- 记录 `run_eventgen_femto_3d.sh` 现在启动的是 `bin/eventgen_femto_3d`，不再指向 `build/` 目录。
- 补记本轮验证证据：本地 `ctest` 与一次真实 `alienv` wrapper smoke。

## 当前代码能力快照

- `Eventgen_femto_3d` 当前仍保持上一轮实现结果：
  - `ApplicationConfig` 已将运行参数和分析参数分离。
  - 入口为 `CLI + TOML`，要求 `--config`，并支持 `--input-root`、`--output-root`、`--input-schema` 覆盖。
  - 输入 schema 支持 `legacy_vector_tree` 与 `blastwave_flat_trees`。
  - `Workflow` 消费归一化后的 `EventData`，不再直接耦合 ROOT 输入树细节。
  - blast-wave 模式读取 `events` + `particles`，并对 `event_id` 连续性、孤儿粒子、缺分支、标量类型做强校验。
  - legacy 单树 `std::vector<>` 输入兼容仍保留。

## 当前验证结论

- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build` 已在当前 `HEAD` 运行：
  - 共注册 `8` 项测试。
  - 其中 `5` 项直接通过。
  - `input_adapter_test`、`legacy_workflow_smoke_test`、`r2_summary_visibility_test` 因未进入 O2/ROOT runtime，被 `run_root_guard.sh` 按设计标记为 `SKIP`，无失败项。
- `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/run_eventgen_femto_3d.sh --config /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /tmp/eventgen_femto_3d_wrapper_smoke_20260413.root` 已在 `alienv` 下运行成功：
  - 成功读取 `5000` 个事件。
  - 成功写出新的 ROOT 输出文件。
  - 说明当前 wrapper 的 `bin/eventgen_femto_3d` 路径切换是有效的。

## 已知剩余事项

- 代码里已经对 blast-wave 标量分支类型做显式校验，但自动化里还没有单独补一个“错误标量类型” fixture 用例。
- `project-state` closeout gate 仍是流程约束而不是工具硬拦截，见 `handoff.md` 与 `work-items.md` 的 `WI-001`。
