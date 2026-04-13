# Current Status

## 2026-04-13 Eventgen 输入适配与 CLI/TOML 配置化

- 状态: 已完成本轮实现与验证收口
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 为 `Eventgen_femto_3d` 增加 `ApplicationConfig`，将运行参数与分析参数分离。
- 新增 `CLI + TOML` 入口，要求 `--config` 必填，并支持 `--input-root`、`--output-root`、`--input-schema` 覆盖。
- 新增 schema-aware 输入层，支持：
  - `legacy_vector_tree`
  - `blastwave_flat_trees`
- `Workflow` 不再直接适配 ROOT 输入树，而是消费归一化的 `EventData`。
- blast-wave 模式读取 `events` + `particles` 两棵树，并对 `event_id` 连续性、孤儿粒子、缺分支、标量类型做强校验。
- 保留 legacy 单树 `std::vector<>` 输入兼容。

## 当前验证结论

- `cmake` 配置通过。
- 全量测试通过：`6/6`。
- wrapper 脚本已做真实参数透传 smoke，能正确把 `--config` 与 CLI 覆盖参数传给二进制。

## 已知剩余事项

- 代码里已经对 blast-wave 标量分支类型做显式校验，但自动化里还没有单独补一个“错误标量类型” fixture 用例。
- 本轮已补回 `project-state`，但也暴露出 closeout 阶段没有被强制执行的流程问题，见 `handoff.md` 与 `work-items.md`。
