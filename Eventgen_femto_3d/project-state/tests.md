# Tests

## 2026-04-13 当前 HEAD 对齐验证记录

- verification_status: `verified`

## 运行命令

- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`
- `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/run_eventgen_femto_3d.sh --config /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /tmp/eventgen_femto_3d_wrapper_smoke_20260413.root`

## 当前构建注册的测试

- `directional_covariance_failure_semantics_test`
- `alpha_hbt_error_semantics_test`
- `fit_summary_directional_error_mask_test`
- `r2_summary_policy_test`
- `config_parse_validation_test`
- `input_adapter_test`
- `legacy_workflow_smoke_test`
- `r2_summary_visibility_test`

## 本轮结果

- `ctest` 共注册 `8` 项测试。
- 直接通过：
  - `directional_covariance_failure_semantics_test`
  - `alpha_hbt_error_semantics_test`
  - `fit_summary_directional_error_mask_test`
  - `r2_summary_policy_test`
  - `config_parse_validation_test`
- 被 `run_root_guard.sh` 按设计标记为 `SKIP`：
  - `input_adapter_test`
  - `legacy_workflow_smoke_test`
  - `r2_summary_visibility_test`
- `alienv` wrapper smoke 成功：
  - 读取 `5000` 个事件
  - 选择 `203546` 个粒子
  - 生成输出 `/tmp/eventgen_femto_3d_wrapper_smoke_20260413.root`
  - 说明 `bin/eventgen_femto_3d` 入口路径与 O2/ROOT 运行时链路可用

## 本轮新增覆盖点

- `CFG-01`: `--config` 缺失、未知参数、schema 合法/非法解析。
- `CFG-02`: CLI 对 `input_root` / `output_root` / `input_schema` 的覆盖优先级。
- `BW-01`: blast-wave `events + particles` 聚合 happy path，含 0 粒子事件。
- `BW-02`: `event_id` 不连续、粒子指向不存在事件的 fail-fast。
- `LEG-01`: legacy 单树向量输入 smoke 与输出契约回归。

## 尚未单独固化的覆盖点

- blast-wave 标量分支类型错误的专门 fixture 还未单独加入自动化。
