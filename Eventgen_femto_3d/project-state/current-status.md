# Current Status

## 2026-06-23 R2 projection fit range 扩大

- 状态: 已完成 R2 拟合输入投影轴范围从 `[-20, 20]` 扩大到 `[-60, 60]`
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 更新默认配置:
  - [AnalysisConfig.h](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/include/femto3d/AnalysisConfig.h)
- 更新运行与示例 TOML:
  - `config/*.toml`
  - `config/examples/*.toml`
- 更新测试:
  - [ConfigParseValidationTest.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/tests/ConfigParseValidationTest.cpp)

## 本轮分析结论

- 当前 R2 拟合没有单独的 fit-window 参数；`ProjectionFit.cpp` 会把一维投影 histogram 中有有效误差的 bin 全部放进 Minuit chi2。
- 一维投影 histogram 的轴来自 `histograms.projection_axis`，此前默认和当前运行配置均为 `[-20, 20]`，因此超出该范围的投影样本只进入 overflow，不参与拟合。
- `Rlong2` 的 Minuit 参数上界是 `1.0e3`，问题不在 R2 参数上界，而在投影输入范围可能截断 long 方向尾部。
- 本轮只扩大 `histograms.projection_axis` 到 `[-60, 60]`，bin width 保持 `0.5`；`rho_out/side/long` 的 3D source 直方图轴暂不扩大，避免 3D bin 数显著增加。

## 本轮验证结论

- O2Physics ROOT executor 首次沙箱运行返回 `STATUS: ESCALATION_REQUIRED`，原因是 `alienv` 入口访问 `/dev/fd` 被 sandbox 拦截。
- 同一构建命令提升权限后返回 `STATUS: PRIMARY_OK`。
- O2Physics ROOT executor 运行 `ctest --test-dir build --output-on-failure` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。

## 2026-06-23 cmake.sh 按需重连

- 状态: 已完成 `scripts/cmake.sh` 自定位、ROOT cache 漂移检测与默认增量构建
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 更新构建脚本:
  - [cmake.sh](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/scripts/cmake.sh)

## 本轮分析结论

- `scripts/cmake.sh` 现在从脚本位置解析 `project_root`，不再硬编码本机绝对源码路径，因此可从任意当前目录调用。
- 默认 build tree 仍为 `${project_root}/build`，可通过 `EVENTGEN_FEMTO_3D_BUILD_DIR` 覆盖。
- 默认并行数仍为 `8`，可通过 `EVENTGEN_FEMTO_3D_BUILD_JOBS` 覆盖。
- 构建默认不再执行 `--clean-first`，而是交给 CMake 增量依赖判断；如需人工全量重建，可显式设置 `EVENTGEN_FEMTO_3D_CLEAN_FIRST=1`。
- 脚本会比较当前 `root-config --prefix` 对应的 `ROOTConfig.cmake` 与 build tree 里缓存的 `ROOT_DIR`；若发现 ROOT runtime 已更新，则清理 ROOT 相关 CMake cache、重新 configure，并 touch 生成的 `link.txt` 让 CMake 触发目标 relink。
- 只要当前 shell 可见 `root-config`，脚本每次都会显式向 CMake 传入当前 `ROOT_DIR`，避免 cache 类型变化后漏检下一次 ROOT 模块升级。
- 本轮附件中的 `TClassTable::Add` / `TCling::LoadPCM` 噪声来自旧 `local7` 链接产物与当前 `local8` wrapper runtime 混用，不是 femto 分析流程自身日志。

## 本轮验证结论

- `bash -n scripts/cmake.sh` 语法检查通过。
- O2Physics ROOT executor 首次运行 `bash scripts/cmake.sh` 检测到 `ROOT/v6-36-10-alice1-local7` -> `ROOT/v6-36-10-alice1-local8` 的 cache 漂移，刷新 ROOT cache 并重新链接输出。
- O2Physics ROOT executor 第二次运行 `bash scripts/cmake.sh` 返回 `STATUS: PRIMARY_OK`，无 `Building CXX object` 或 `Linking CXX executable` 输出，确认无更新时不会重复重编/重连。
- `build/CMakeFiles/eventgen_femto_3d.dir/link.txt` 与 `bin/eventgen_femto_3d` 的 `LC_RPATH` 均指向 `ROOT/v6-36-10-alice1-local8/lib`。
- O2Physics ROOT executor 运行 `ctest --test-dir build --output-on-failure` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- 真实 wrapper smoke 写出 `/private/tmp/eventgen_femto_3d_root_noise_probe.root`，读取 `5000` 个事件、选择 `203546` 个粒子、接受 `3249284` 个 pair；输出中未再出现 `TClassTable::Add` 或 `TCling::LoadPCM` 噪声。

## 2026-06-10 projection-fit alpha bounds 配置化

- 状态: 已完成 `projection_fit.alpha_min` / `projection_fit.alpha_max` 配置化
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 更新配置解析:
  - [Config.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/Config.cpp)
  - [ProjectionFit.h](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/include/femto3d/ProjectionFit.h)
- 更新拟合实现:
  - [ProjectionFit.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/ProjectionFit.cpp)
- 更新输出 metadata:
  - [Workflow.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/Workflow.cpp)
- 更新测试与示例配置:
  - [ConfigParseValidationTest.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/tests/ConfigParseValidationTest.cpp)
  - `config/*.toml`
  - `config/examples/*.toml`

## 本轮分析结论

- `[projection_fit]` 现在可选设置 `alpha_min` 与 `alpha_max`；未设置时仍使用默认 `0.2` 与 `2.0`。
- 配置校验要求两个值有限，且满足 `0 < alpha_min < alpha_max`。
- Minuit alpha 参数边界改为来自 `ProjectionFitConfig`，初值与步长会根据配置区间安全裁剪。
- `projection_fit_metadata` 新增 `alpha_min` 与 `alpha_max` branch，用于结果文件追溯。
- 物理解释上，默认 `alpha_max = 2` 保持 Levy-stable alpha 的物理区间上界；若显式设置超过 `2`，应视为数值诊断扩展，不应直接当作物理 Levy-stable 结果解释。

## 本轮验证结论

- `cmake --build build` 本地构建成功。
- 本地 `ctest --test-dir build --output-on-failure` 通过；ROOT runtime guard 测试在普通 shell 中按设计跳过。
- O2Physics ROOT executor `ctest --test-dir build --output-on-failure` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- 临时 smoke 输出 `/private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root` 成功读取 `5000` 个事件并接受 `3249284` 个 pair。
- ROOT inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`，确认 `projection_fit_metadata` 有 `8` 个 branch。
- ROOT metadata 查询确认 `alpha_min=0.20000000000000001`、`alpha_max=2`。

## 2026-06-10 source-parameter overview canvas

- 状态: 已完成按 `cent/femto_mt` 分区写出 source-parameter overview canvas
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 更新分析输出:
  - [Workflow.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/Workflow.cpp)
- 更新测试:
  - [R2SummaryVisibilityTest.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/tests/R2SummaryVisibilityTest.cpp)

## 本轮分析结论

- 每个有可画 summary 点的 `R2Summary/<centrality>/<femto_mt>/` 现在写出 `source_parameters_overview_canvas`。
- overview canvas 使用 `2 x 4` 排列:
  - pad 1: centrality、femto `m_T` 与相对角说明
  - pad 2: `alpha` vs `phi_pair - Psi_2`
  - pads 3-8: `Rout2`、`Ros2`、`Rside2`、`Rol2`、`Rlong2`、`Rsl2`
- `alpha` 点来自当前 slice 的 `fit_products.alpha`，只作为 overview canvas 内部图元使用；没有新增独立 `alpha_vs_phi` ROOT key。
- 既有 `R2Summary/<centrality>/<femto_mt>/{Rout2,Rside2,Rlong2,Ros2,Rol2,Rsl2}_vs_phi` graph/canvas 与 `R2Summary/<centrality>/epsf_vs_mt` 输出保持兼容。

## 本轮验证结论

- O2Physics ROOT executor 构建返回 `STATUS: PRIMARY_OK`。
- O2Physics ROOT executor `ctest --test-dir build --output-on-failure` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- 临时 smoke 输出 `/private/tmp/eventgen_femto_3d_overview_smoke.root` 成功读取 `5000` 个事件并接受 `3249284` 个 pair。
- ROOT inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`，确认两个 mT 目录均写出 `source_parameters_overview_canvas`。
- 已导出 `/private/tmp/eventgen_femto_3d_overview_smoke.png`，确认 canvas 非空且顺序为信息面板、`alpha`、`Rout/Ros`、`Rside/Rol`、`Rlong/Rsl`。

## 2026-06-09 进度条与 epsf-vs-mT 输出图

- 状态: 已完成 CLI 进度条与 `epsf_vs_mt` 结果图写出
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 新增源码:
  - [ProgressReporter.h](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/include/femto3d/ProgressReporter.h)
  - [ProgressReporter.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/ProgressReporter.cpp)
- 更新入口和配置:
  - [Config.h](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/include/femto3d/Config.h)
  - [Config.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/Config.cpp)
  - [main.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/main.cpp)
- 更新分析输出:
  - [Workflow.h](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/include/femto3d/Workflow.h)
  - [Workflow.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/src/Workflow.cpp)
- 更新测试:
  - [R2SummaryPolicyTest.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/tests/R2SummaryPolicyTest.cpp)
  - [R2SummaryVisibilityTest.cpp](/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/tests/R2SummaryVisibilityTest.cpp)

## 本轮分析结论

- 进度条参考 `Blast_wave` 的实现方式，提供百分比、活动帧和 ETA；默认 `auto`，并新增 `--progress` / `--no-progress` CLI 开关。
- `RunAnalysis` 仍保留旧调用方式；新增可选 `AnalysisProgressSink*`，由入口层传入 `ProgressReporter`。
- 新图只使用当前项目已经写出的 `Rside2_vs_phi` summary 点，不引入新的输入树或外部图像依赖。
- `epsf_vs_mt` 写入位置为 `R2Summary/<centrality>/epsf_vs_mt`，对应 canvas 为 `R2Summary/<centrality>/epsf_vs_mt_canvas`。
- 本轮临时 smoke 输出 `/private/tmp/eventgen_femto_3d_epsf_smoke.root` 中确认 `R2Summary/cent_0.00-100.00/epsf_vs_mt` 为 `TGraphErrors`，点数为 `2`。

## 本轮验证结论

- 沙箱内直接进入 O2Physics ROOT 环境仍会触发 `alienv` `/dev/fd` 权限问题；按既有规则使用授权非沙箱 executor。
- `cmake --build /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/build -j 8` 返回 `STATUS: PRIMARY_OK`。
- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/build` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- 临时 smoke 命令强制 `--progress` 成功读取 `5000` 个事件，输出 `/private/tmp/eventgen_femto_3d_epsf_smoke.root`。
- ROOT inspector 检查临时输出返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`，并确认新增 `epsf_vs_mt` graph/canvas。

## 2026-06-05 dense-mix Glauber 15-phi-bin 分析配置与结果

- 状态: 已完成 15-phi-bin 配置更新、一键脚本分析与输出 ROOT 结构验证
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 更新配置:
  - [bw_dense_mix_glauber_15phibin copy.toml](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/bw_dense_mix_glauber_15phibin%20copy.toml)
- 新增一键运行脚本:
  - [run_dense_mix_glauber_15phibin.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh)
- 脚本内部固定输入:
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b1_dense_mix_glauber.root`
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b3_dense_mix_glauber.root`
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b8_dense_mix_glauber.root`
- 新增输出:
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b1_dense_mix_glauber_15phibin.root`
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b3_dense_mix_glauber_15phibin.root`
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b8_dense_mix_glauber_15phibin.root`

## 本轮分析结论

- 15-bin 配置覆盖 `[-pi/2, pi/2]`，均匀切成 `15` 个 phi bin；当前仍保留 `2` 个 femto-mT bin，因此每个输出 `fit_summary` 有 `30` 行。
- 当前 dense-mix Glauber 输入文件实际各读到 `50000` 个事件，和 2026-06-04 7-bin 台账中的旧 `5000` 事件样本不同。
- `PbPb_b1_dense_mix_glauber_15phibin.root`:
  - `events_read=50000`，`selected_particles=4242192`。
  - `events_with_valid_event_plane=50000`。
  - `candidate_pairs=180284061`，`accepted_pairs=104060948`。
  - `fit_summary` 有 `30` 行，`unique_phi_bins=15`，HBT xyz 与 alpha 成功行数均为 `30/30`。
- `PbPb_b3_dense_mix_glauber_15phibin.root`:
  - `events_read=50000`，`selected_particles=3893172`。
  - `events_with_valid_event_plane=50000`。
  - `candidate_pairs=151957154`，`accepted_pairs=87782661`。
  - `fit_summary` 有 `30` 行，`unique_phi_bins=15`，HBT xyz 与 alpha 成功行数均为 `30/30`。
- `PbPb_b8_dense_mix_glauber_15phibin.root`:
  - `events_read=50000`，`selected_particles=1884790`。
  - `events_with_valid_event_plane=50000`。
  - `candidate_pairs=35957179`，`accepted_pairs=20787715`。
  - `fit_summary` 有 `30` 行，`unique_phi_bins=15`，HBT xyz 与 alpha 成功行数均为 `30/30`。
- 三个输出 ROOT 文件均包含 `Femto3D`、`R2Summary`、`fit_summary`、`analysis_statistics` 与 metadata 树，ROOT inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
- 三个输出结构均有 `TOP_LEVEL_KEYS=8`，对象摘要为 `directories=38`、`histograms=210`、`trees=6`、`graphs=12`、`canvases=42`、`total_keys_seen=488`。
- 本轮运行未出现 close-pair、pair-kinematics、event-plane 或 R2 summary invalid-error 跳过。

## 本轮验证结论

- 沙箱内直接运行 `Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh` 会触发 `alienv` `/dev/fd` 权限问题和 ROOT module loading 级联失败；同一脚本在授权的非沙箱环境下完整跑完 b1、b3、b8。
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh unexpected_arg` 按预期打印 `Usage` 并返回 exit code `2`。
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d --command 'root -l -b -q /private/tmp/summarize_eventgen_femto3d_15phibin.C'` 返回 `STATUS: PRIMARY_OK`。
- 本轮未修改 C++ 源码、构建系统或自动化测试注册表，因此未重跑 `ctest`；配置解析与运行路径已由真实 15-bin wrapper 分析覆盖。

## 2026-06-04 dense-mix Glauber 7-phi-bin 分析配置与结果

- 状态: 已新增并验证本轮 dense-mix Glauber 输入分析配置
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 新增配置:
  - [bw_dense_mix_glauber_7phibin.toml](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/bw_dense_mix_glauber_7phibin.toml)
- 新增一键运行脚本:
  - [run_dense_mix_glauber_7phibin.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh)
- 默认输入:
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b1_dense_mix_glauber.root`
- 脚本内部显式配置第二个输入:
  - `/Users/allenzhou/Research_software/Blast_wave/res/PbPb_b8_dense_mix_glauber.root`
- 新增输出:
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b1_dense_mix_glauber_7phibin.root`
  - `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b8_dense_mix_glauber_7phibin.root`

## 本轮分析结论

- 两个输入 ROOT 文件均可由 O2Physics ROOT 读取，`events` 树各含 `5000` 个事件，`particles` 树分支类型匹配 `blastwave_flat_trees` 输入契约。
- `PbPb_b1_dense_mix_glauber_7phibin.root`:
  - 读取 `5000` 个事件，选择 `389228` 个粒子。
  - `5000` 个事件都有有效 input `psi2` event plane。
  - `candidate_pairs=15188363`，`accepted_pairs=8770036`。
  - `fit_summary` 有 `14` 行，HBT xyz 与 alpha 拟合成功行数均为 `14/14`。
- `PbPb_b8_dense_mix_glauber_7phibin.root`:
  - 读取 `5000` 个事件，选择 `188980` 个粒子。
  - `5000` 个事件都有有效 input `psi2` event plane。
  - `candidate_pairs=3611138`，`accepted_pairs=2086729`。
  - `fit_summary` 有 `14` 行，HBT xyz 与 alpha 拟合成功行数均为 `14/14`。
- 两个输出 ROOT 文件均包含 `Femto3D`、`R2Summary`、`fit_summary`、`analysis_statistics`、event-plane/mT/projection-fit metadata，结构检查返回 `STATUS: OK`。
- 本轮运行未出现 close-pair、pair-kinematics、event-plane 或 R2 summary invalid-error 跳过。

## 本轮验证结论

- 当前 `bin/eventgen_femto_3d` 初次运行时暴露了旧 ROOT runtime 链接问题: `dyld` 找不到旧路径下的 `libMinuit2.6.36.so`。
- 已在 `O2Physics/latest-master-o2` 下重新运行 CMake configure 与 build，两个命令均返回 `STATUS: PRIMARY_OK`。
- 修正 `scripts/run_eventgen_femto_3d.sh` 的项目根目录解析后，`scripts/run_dense_mix_glauber_7phibin.sh` 可在无参数情况下连续跑完 b1 与 b8。
- 一键脚本若收到额外命令行参数，会打印 `Usage` 并以 exit code `2` 退出，避免外部路径覆盖脚本内置配置。
- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`:
  - 共注册 `8` 项测试。
  - `5` 项直接通过。
  - `input_adapter_test`、`legacy_workflow_smoke_test`、`r2_summary_visibility_test` 因当前 shell 未进入 O2/ROOT runtime，被 `run_root_guard.sh` 按设计标记为 `SKIP`。
  - 无失败项。

## 2026-04-13 project-state 对齐更新与 wrapper 路径核对

- 状态: 已完成本轮台账同步与最新证据补记
- verification_status: `verified`
- project_state_sync_status: `written`

## 本轮范围

- 对照当前 `HEAD` 校对 `project-state` 与代码/构建输出是否一致。
- 修正文档中关于测试数量的旧描述：从早先的 `6/6` 更新为当前构建注册的 `8` 项测试。
- 记录 wrapper 现在启动的是 `bin/eventgen_femto_3d`，不再指向 `build/` 目录。
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
- `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_eventgen_femto_3d.sh --config /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /tmp/eventgen_femto_3d_wrapper_smoke_20260413.root` 是当前等价 wrapper smoke 命令：
  - 成功读取 `5000` 个事件。
  - 成功写出新的 ROOT 输出文件。
  - 说明当前 wrapper 的 `bin/eventgen_femto_3d` 路径切换是有效的。

## 已知剩余事项

- 代码里已经对 blast-wave 标量分支类型做显式校验，但自动化里还没有单独补一个“错误标量类型” fixture 用例。
- `project-state` closeout gate 仍是流程约束而不是工具硬拦截，见 `handoff.md` 与 `work-items.md` 的 `WI-001`。
