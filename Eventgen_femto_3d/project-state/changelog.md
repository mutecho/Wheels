# Changelog

## 2026-06-10

- 将 `projection_fit.alpha_min` 与 `projection_fit.alpha_max` 暴露到 TOML 配置；未设置时默认仍为 `0.2` 与 `2.0`。
- `projection_fit_metadata` 新增写出 `alpha_min` 与 `alpha_max`，方便从结果 ROOT 文件追溯本次 alpha 拟合边界。
- 扩展 `config_parse_validation_test`，覆盖 alpha bounds 默认值、TOML 覆盖以及非法区间 `0 < alpha_min < alpha_max` 校验。
- 在每个 `R2Summary/<centrality>/<femto_mt>/` 下新增 `source_parameters_overview_canvas`，按参考图顺序排列 `alpha`、`Rout2`、`Ros2`、`Rside2`、`Rol2`、`Rlong2`、`Rsl2`。
- 新增 `alpha` 的 per-phi summary 点收集，仅用于 overview canvas；已有单图 `Rout2_vs_phi` 等 graph/canvas 与 `epsf_vs_mt` 输出路径保持不变。
- 扩展 `r2_summary_visibility_test`，递归检查 overview canvas 中的信息面板和 6 个半径图 primitive。
- 在 O2Physics ROOT executor 下完成构建、完整 `ctest`、临时 smoke 输出 ROOT 检查，确认 overview canvas 与 alpha bounds metadata 写出，并导出 overview canvas PNG 做版式确认。

## 2026-06-09

- 新增 CLI 进度条支持，默认 `auto` 模式只在 stderr 是 TTY 时显示，`--progress` 可强制开启，`--no-progress` 可关闭。
- 新增 `ProgressReporter`，事件循环通过 `AnalysisProgressSink` 回调更新进度，分析核心不直接依赖终端渲染。
- 在 `R2Summary/<centrality>/` 下新增 `epsf_vs_mt` 与 `epsf_vs_mt_canvas`，从当前已有 `Rside2_vs_phi` summary 点提取 `#epsilon_{f} = 2R_{s,2}^{2}/R_{s,0}^{2}` 随 femto `m_{T}` 的变化。
- 扩展 `r2_summary_policy_test`，覆盖进度条格式和 `Rside2(phi)` 二阶调制到 `epsf` 的解析提取。
- 将 `r2_summary_visibility_test` 的 marker-size 持久化检查改为 `1e-6` 容差，避免 ROOT 浮点持久化精度导致误报。
- 在 O2Physics ROOT executor 下完成构建、完整 `ctest` 和临时输出 ROOT 检查。

## 2026-06-05

- 将 `config/bw_dense_mix_glauber_15phibin copy.toml` 更新为 `15` 个 phi bin，默认输出路径同步改为 `_15phibin.root`。
- 新增 `scripts/run_dense_mix_glauber_15phibin.sh`，仿照 7-bin 脚本固定 config/input/output 路径并连续运行 b1、b3、b8。
- 真实运行 15-bin dense-mix Glauber 分析并写出 `res/PbPb_b1_dense_mix_glauber_15phibin.root`、`res/PbPb_b3_dense_mix_glauber_15phibin.root`、`res/PbPb_b8_dense_mix_glauber_15phibin.root`。
- 用 ROOT inspector 验证三个 15-bin 输出结构可读，三者均为 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
- 用 ROOT 汇总宏确认三个输出均有 `fit_summary` `30` 行、`unique_phi_bins=15`，HBT xyz 与 alpha 成功行数均为 `30/30`。
- 记录沙箱内 `alienv` `/dev/fd` 权限问题；同一 15-bin 脚本在授权非沙箱环境下成功。

## 2026-06-04

- 新增 `config/bw_dense_mix_glauber_7phibin.toml`，复用 `bw_7phibin_test.toml` 的 7 个 phi bin、2 个 femto-mT bin、input `psi2` event-plane 与拟合容错设置。
- 新增 `scripts/run_dense_mix_glauber_7phibin.sh`，脚本内部固定 config/input/output 路径，调用时不需要额外命令行参数。
- 修正 `scripts/run_eventgen_femto_3d.sh` 的项目根目录解析，使 wrapper 从 `scripts/` 下也能正确找到 `bin/eventgen_femto_3d`。
- 对 `PbPb_b1_dense_mix_glauber.root` 与 `PbPb_b8_dense_mix_glauber.root` 完成真实一键脚本分析。
- 写出 `res/PbPb_b1_dense_mix_glauber_7phibin.root` 与 `res/PbPb_b8_dense_mix_glauber_7phibin.root`，并用 ROOT inspector 验证输出结构可读。
- 记录本轮旧 ROOT runtime 动态库链接问题，并在 `O2Physics/latest-master-o2` 下重新 configure/build。
- 补记本轮 `ctest`、输入 ROOT 检查、输出 ROOT 检查与 `fit_summary` 摘要证据。

## 2026-04-13

- 完成 `Eventgen_femto_3d` 的 `CLI + TOML` 配置化改造。
- 完成 blast-wave `events + particles` 输入适配，并保留 legacy 输入兼容。
- 将输入层从 `Workflow` 中解耦，统一归一化为 `EventData`。
- 新增配置解析、输入适配、legacy smoke 与 ROOT runtime guard 测试。
- 在 O2/ROOT 环境下完成构建、全量测试与 wrapper 脚本 smoke。
- 补写 `project-state`，记录本轮实现结果与流程改进项。
- 补齐 `project-state/guide.md`，完成面向维护者的人类可读概览。
- 调整了可执行二进制文件的位置到 `bin` 目录下。
- 同步 `project-state` 到当前 `HEAD`，将旧的 `6/6` 测试描述更新为当前注册的 `8` 项测试。
- 记录 `run_eventgen_femto_3d.sh` 已切换为调用 `bin/eventgen_femto_3d`。
- 补记一次当前 `HEAD` 下的 `ctest` 与 `alienv` wrapper smoke 证据。
