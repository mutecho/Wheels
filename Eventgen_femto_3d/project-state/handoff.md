# Handoff

## 最新交接

- 交接时间: 2026-06-10 CST
- 当前 owner: 父线程主执行者
- 下一步 owner: 任意继续维护该模块的执行者
- verification_status: `verified`
- project_state_sync_status: `written`

## 已完成事项

- 已将 `projection_fit.alpha_min` 与 `projection_fit.alpha_max` 暴露到 TOML 配置。
  - 未设置时默认仍为 `0.2` 与 `2.0`。
  - 配置校验要求 `0 < alpha_min < alpha_max`。
  - Minuit alpha 参数边界现在来自配置，而不是 `ProjectionFit.cpp` 内部硬编码。
  - `projection_fit_metadata` 会写出 `alpha_min` 与 `alpha_max`，便于结果文件追溯。
  - 临时 smoke 输出 `/private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root` 已确认 metadata 为 `alpha_min=0.20000000000000001`、`alpha_max=2`。
- 已新增 `R2Summary/<centrality>/<femto_mt>/source_parameters_overview_canvas`。
  - canvas 按 `2 x 4` 排列，同一 centrality 与 femto `m_T` 区间下依次显示信息面板、`alpha`、`Rout2`、`Ros2`、`Rside2`、`Rol2`、`Rlong2`、`Rsl2`。
  - `alpha` 来自每个 `cent/femto_mt/phi` slice 的 `fit_products.alpha`，只写入 overview canvas 内部，不新增独立 `alpha_vs_phi` ROOT key。
  - 既有单图 `Rout2_vs_phi` 等 graph/canvas 与 `epsf_vs_mt` 输出路径保持兼容。
- 已用临时 smoke 输出确认:
  - `/private/tmp/eventgen_femto_3d_overview_smoke.root`
  - `R2Summary/cent_0.00-100.00/femto_mt_0.20-0.40/source_parameters_overview_canvas` 已写出。
  - `R2Summary/cent_0.00-100.00/femto_mt_0.40-0.60/source_parameters_overview_canvas` 已写出。
  - `/private/tmp/eventgen_femto_3d_overview_smoke.png` 已导出并目视确认顺序。
- 已在 O2Physics ROOT executor 下完成:
  - CMake build: `STATUS: PRIMARY_OK`
  - `ctest`: `8/8` 通过，`STATUS: PRIMARY_OK`
  - ROOT inspector: `RUNTIME_STATUS: PRIMARY_OK`，`STATUS: OK`
- 已新增 CLI 进度条:
  - 默认 `auto`，stderr 是 TTY 时显示。
  - `--progress` 强制开启。
  - `--no-progress` 关闭。
- 已新增 `R2Summary/<centrality>/epsf_vs_mt` 和 `epsf_vs_mt_canvas`。
  - 当前定义来自已有 `Rside2_vs_phi` summary 点的二阶调制:
    `#epsilon_{f} = 2R_{s,2}^{2}/R_{s,0}^{2}`。
  - 不读取或依赖 Blast_wave 的 `eps2_f` event tree 分支。
- 已用临时 smoke 输出确认:
  - `/private/tmp/eventgen_femto_3d_epsf_smoke.root`
  - `R2Summary/cent_0.00-100.00/epsf_vs_mt` 是 `TGraphErrors`，点数为 `2`。
- 已将 dense-mix Glauber 15-phi-bin 配置落到:
  - [bw_dense_mix_glauber_15phibin copy.toml](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/bw_dense_mix_glauber_15phibin%20copy.toml)
- 新增 dense-mix Glauber 15-phi-bin 一键运行脚本:
  - [run_dense_mix_glauber_15phibin.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh)
- 已用一键脚本完成三个真实输入分析:
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b1_dense_mix_glauber.root`
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b3_dense_mix_glauber.root`
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b8_dense_mix_glauber.root`
- 已生成并检查三个输出 ROOT:
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b1_dense_mix_glauber_15phibin.root`
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b3_dense_mix_glauber_15phibin.root`
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b8_dense_mix_glauber_15phibin.root`
- 分析统计:
  - b1: `events_read=50000`，`selected_particles=4242192`，`accepted_pairs=104060948`。
  - b3: `events_read=50000`，`selected_particles=3893172`，`accepted_pairs=87782661`。
  - b8: `events_read=50000`，`selected_particles=1884790`，`accepted_pairs=20787715`。
  - 三个输出 `fit_summary` 均有 `30` 行，`unique_phi_bins=15`，HBT xyz 与 alpha 成功行数均为 `30/30`。
- 三个输出 ROOT inspector 均返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
- 沙箱内直接运行 15-bin 脚本会触发 `alienv` `/dev/fd` 权限问题；授权非沙箱运行已成功。
- 本轮未修改 C++ 源码、构建系统或测试注册表，因此未重跑 `ctest`。
- 之前实现的 blast-wave 输入适配、legacy 输入兼容、`CLI + TOML` 入口以及相关测试覆盖仍然成立。

## 下一步注意事项

- 若要复现 15-phi-bin 结果，直接运行 `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh`，不需要额外命令行参数。
- 若需要复现上一轮 7-phi-bin 结果，运行 `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh`。
- 如果 O2Physics ROOT 模块更新后再次出现 `dyld` 找不到 ROOT 动态库，先在 O2Physics 环境下重新运行 CMake configure/build，再判断分析问题。
- 若要强制查看分析进度，给主程序或 wrapper 透传 `--progress`；若日志环境不想出现进度条，透传 `--no-progress`。
- 本轮新增 `epsf_vs_mt` 只在某个 centrality 下至少一个 mT bin 有足够有效 `Rside2_vs_phi` 点时写出。
- 本轮新增 `source_parameters_overview_canvas` 只在对应 `cent/femto_mt` 下存在可画 summary 点时写出。

## 流程提醒

- 在任务开始时，只要仓库里已存在 `project-state/`，父线程就应立即把“结束前写回 `project-state`”放进计划，并作为正常完成前的必过步骤。
- 父线程在最终总结前应显式计算三件事：
  - 是否存在 `project-state/`
  - 本轮是否发生了状态变化
  - 本轮是否已经完成 `project-state` 写回
- 若前两项为真而第三项为假，应先写文档，再结束任务，而不是直接给最终答复。
- 多代理场景下，子代理的 handback 需要固定包含：
  - 验证结果
  - durable decisions
  - durable issues
  - next step hint
  父线程据此统一写回，不能指望子代理自己补文档。

## 下一步建议

- 先补一条 blast-wave 错误标量类型的自动化 fixture，用来覆盖 `BW-03`。
- 后续若再改 `scripts/run_eventgen_femto_3d.sh`、`CMakeLists.txt` 或测试注册表，应同步检查 `tests.md` 与 `guide.md` 是否也需要更新。
