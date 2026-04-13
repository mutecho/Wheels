# Tests

## 2026-04-13 验证记录

- verification_status: `verified`

## 运行命令

- `cmake -S /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d -B /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`
- `cmake --build /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`
- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`
- `run_eventgen_femto_3d.sh --config ... --input-root ... --output-root ...`

## 通过的测试

- `directional_covariance_failure_semantics_test`
- `alpha_hbt_error_semantics_test`
- `fit_summary_directional_error_mask_test`
- `config_parse_validation_test`
- `input_adapter_test`
- `legacy_workflow_smoke_test`

## 本轮新增覆盖点

- `CFG-01`: `--config` 缺失、未知参数、schema 合法/非法解析。
- `CFG-02`: CLI 对 `input_root` / `output_root` / `input_schema` 的覆盖优先级。
- `BW-01`: blast-wave `events + particles` 聚合 happy path，含 0 粒子事件。
- `BW-02`: `event_id` 不连续、粒子指向不存在事件的 fail-fast。
- `LEG-01`: legacy 单树向量输入 smoke 与输出契约回归。

## 尚未单独固化的覆盖点

- blast-wave 标量分支类型错误的专门 fixture 还未单独加入自动化。
