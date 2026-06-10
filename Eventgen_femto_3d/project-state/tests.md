# Tests

## 2026-06-10 projection-fit alpha bounds 配置化验证

- verification_status: `verified`

## 运行命令

- `cmake --build build`
- `ctest --test-dir build --output-on-failure`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'ctest --test-dir build --output-on-failure'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'bin/eventgen_femto_3d --config /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root --no-progress'`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root --max-keys 260`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'root -l -b -q -e "TFile f(\"/private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root\"); auto t = dynamic_cast<TTree*>(f.Get(\"projection_fit_metadata\")); if (!t) { gSystem->Exit(2); } double amin=0,amax=0; t->SetBranchAddress(\"alpha_min\", &amin); t->SetBranchAddress(\"alpha_max\", &amax); t->GetEntry(0); printf(\"alpha_min=%.17g alpha_max=%.17g\\n\", amin, amax); gSystem->Exit((std::abs(amin-0.2)<1e-12 && std::abs(amax-2.0)<1e-12) ? 0 : 3);"'`

## 本轮结果

- `config_parse_validation_test` 新增覆盖:
  - 未设置 `alpha_min` / `alpha_max` 时默认值为 `0.2` / `2.0`。
  - TOML 可覆盖 alpha bounds。
  - 非法区间不满足 `0 < alpha_min < alpha_max` 时抛出 `ConfigError`。
- 本地构建成功。
- 本地普通 `ctest` 通过；ROOT guard 测试在未进入 O2/ROOT runtime 时按设计跳过。
- O2Physics ROOT executor `ctest` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- 临时 smoke:
  - 读取 `5000` 个事件，选择 `203546` 个粒子，接受 `3249284` 个 pair。
  - 输出 `/private/tmp/eventgen_femto_3d_alpha_bounds_smoke.root`。
- 输出 ROOT 检查:
  - inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
  - `projection_fit_metadata` 从 `6` 个 branch 扩展为 `8` 个 branch。
  - metadata 查询返回 `alpha_min=0.20000000000000001 alpha_max=2`。

## 2026-06-10 source-parameter overview canvas 验证

- verification_status: `verified`

## 运行命令

- `cmake --build build`
- `ctest --test-dir build --output-on-failure`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'cmake --build build'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'ctest --test-dir build --output-on-failure'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'bin/eventgen_femto_3d --config /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /private/tmp/eventgen_femto_3d_overview_smoke.root --no-progress'`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /private/tmp/eventgen_femto_3d_overview_smoke.root --max-keys 260`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'root -l -b -q -e "TFile f(\"/private/tmp/eventgen_femto_3d_overview_smoke.root\"); auto c = dynamic_cast<TCanvas*>(f.Get(\"R2Summary/cent_0.00-100.00/femto_mt_0.20-0.40/source_parameters_overview_canvas\")); if (!c) { gSystem->Exit(2); } c->Print(\"/private/tmp/eventgen_femto_3d_overview_smoke.png\");"'`

## 本轮结果

- 本地构建成功。
- 沙箱内普通 `ctest` 显示 `8/8` 通过，其中 ROOT guard 测试被跳过；该结果不作为 ROOT 输出验证的最终证据。
- O2Physics ROOT executor 构建返回 `STATUS: PRIMARY_OK`。
- O2Physics ROOT executor `ctest` 返回 `STATUS: PRIMARY_OK`，`8/8` 通过。
- `r2_summary_visibility_test` 新增覆盖:
  - `R2Summary/<centrality>/<femto_mt>/source_parameters_overview_canvas` 存在。
  - canvas 内含 `source_parameters_overview_info_box`。
  - canvas 内含 `Rout2_vs_phi`、`Ros2_vs_phi`、`Rside2_vs_phi`、`Rol2_vs_phi`、`Rlong2_vs_phi`、`Rsl2_vs_phi`。
- 临时 smoke:
  - 读取 `5000` 个事件，选择 `203546` 个粒子，接受 `3249284` 个 pair。
  - 输出 `/private/tmp/eventgen_femto_3d_overview_smoke.root`。
- 输出 ROOT 检查:
  - inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
  - `R2Summary/cent_0.00-100.00/femto_mt_0.20-0.40/source_parameters_overview_canvas` 已写出。
  - `R2Summary/cent_0.00-100.00/femto_mt_0.40-0.60/source_parameters_overview_canvas` 已写出。
  - 对象摘要为 `directories=16`、`histograms=56`、`trees=6`、`graphs=13`、`canvases=23`、`total_keys_seen=162`。
- 已导出 `/private/tmp/eventgen_femto_3d_overview_smoke.png` 并确认版式顺序与任务要求一致。

## 2026-06-09 进度条与 epsf-vs-mT 输出验证

- verification_status: `verified`

## 运行命令

- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'cmake --build /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/build -j 8'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command 'ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/build'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d --command '/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/bin/eventgen_femto_3d --config /Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /private/tmp/eventgen_femto_3d_epsf_smoke.root --progress'`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /private/tmp/eventgen_femto_3d_epsf_smoke.root --max-keys 220`

## 本轮结果

- 构建:
  - 沙箱内首次 ROOT/O2Physics executor 返回 `STATUS: ESCALATION_REQUIRED`，原因为 `alienv` `/dev/fd` 权限问题。
  - 授权非沙箱重跑后，构建返回 `STATUS: PRIMARY_OK`。
- `ctest`:
  - 授权非沙箱 O2Physics ROOT 环境下 `8/8` 通过。
  - 新覆盖点包括 `FormatProgressLine` 百分比/ETA 输出，以及 `Rside2(phi)=10+2cos(2phi)` 对应 `epsf=0.2` 的提取。
- 临时 smoke:
  - `--progress` 强制输出进度条，运行成功。
  - 读取 `5000` 个事件，选择 `203546` 个粒子，接受 `3249284` 个 pair。
  - 输出 `/private/tmp/eventgen_femto_3d_epsf_smoke.root`。
- 输出 ROOT 检查:
  - inspector 返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
  - `R2Summary/cent_0.00-100.00/epsf_vs_mt` 为 `TGraphErrors`，点数为 `2`。
  - `R2Summary/cent_0.00-100.00/epsf_vs_mt_canvas` 已写出。

## 2026-06-05 dense-mix Glauber 15-phi-bin 分析验证

- verification_status: `verified`

## 运行命令

- `awk 'BEGIN{n=0} /\[\[bins\.phi\]\]/{n++} END{print n}' '/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/bw_dense_mix_glauber_15phibin copy.toml'`
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh`
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh unexpected_arg`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b1_dense_mix_glauber_15phibin.root --max-keys 120`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b3_dense_mix_glauber_15phibin.root --max-keys 60`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b8_dense_mix_glauber_15phibin.root --max-keys 60`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d --command 'root -l -b -q /private/tmp/summarize_eventgen_femto3d_15phibin.C'`

## 本轮结果

- 配置检查:
  - `bw_dense_mix_glauber_15phibin copy.toml` 中 `[[bins.phi]]` 数量为 `15`。
  - 输出默认路径已从 `_7phibin.root` 改为 `_15phibin.root`。
- 一键脚本:
  - 沙箱内直接运行会触发 `alienv` `/dev/fd` 权限问题和 ROOT module loading 级联失败。
  - 授权非沙箱运行 `Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh` 成功完成 b1、b3、b8 三个输入分析。
  - `unexpected_arg` 调用按预期打印 `Usage` 并返回 exit code `2`。
- b1 分析输出:
  - `events_read=50000`
  - `selected_particles=4242192`
  - `accepted_pairs=104060948`
  - `fit_summary` 行数 `30`，`unique_phi_bins=15`，HBT xyz 成功 `30/30`，alpha 成功 `30/30`
- b3 分析输出:
  - `events_read=50000`
  - `selected_particles=3893172`
  - `accepted_pairs=87782661`
  - `fit_summary` 行数 `30`，`unique_phi_bins=15`，HBT xyz 成功 `30/30`，alpha 成功 `30/30`
- b8 分析输出:
  - `events_read=50000`
  - `selected_particles=1884790`
  - `accepted_pairs=20787715`
  - `fit_summary` 行数 `30`，`unique_phi_bins=15`，HBT xyz 成功 `30/30`，alpha 成功 `30/30`
- 输出 ROOT 检查:
  - 三个输出均返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
  - 三个输出均包含 `Femto3D`、`R2Summary`、`fit_summary`、`analysis_statistics` 和 metadata 树。
  - 三个输出对象摘要一致: `directories=38`、`histograms=210`、`trees=6`、`graphs=12`、`canvases=42`、`total_keys_seen=488`。
- 本轮未修改 C++ 源码、构建系统或测试注册表，未重跑 `ctest`；真实 wrapper 分析已覆盖本轮配置解析与运行路径。

## 2026-06-04 dense-mix Glauber 7-phi-bin 分析验证

- verification_status: `verified`

## 运行命令

- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Blast_wave/res/PbPb_b1_dense_mix_glauber.root --max-keys 120`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Blast_wave/res/PbPb_b8_dense_mix_glauber.root --max-keys 120`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d --command 'cmake -S /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d -B /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build'`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d --command 'cmake --build /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build'`
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh`
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh unexpected_arg`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b8_dense_mix_glauber_7phibin.root --max-keys 120`
- `bash /Users/allenzhou/.codex/skills/cern_root/root-file-inspector/scripts/run_inspect_root_file.sh /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/res/PbPb_b1_dense_mix_glauber_7phibin.root --max-keys 80`
- `bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh --cwd /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d --command 'root -l -b -q /private/tmp/summarize_eventgen_femto3d.C'`
- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`

## 本轮结果

- 输入 ROOT 检查:
  - `PbPb_b1_dense_mix_glauber.root`: `RUNTIME_STATUS: PRIMARY_OK`，`STATUS: OK`，`events` 有 `5000` entries，最高 cycle `particles` 有 `4054118` entries。
  - `PbPb_b8_dense_mix_glauber.root`: `RUNTIME_STATUS: PRIMARY_OK`，`STATUS: OK`，`events` 有 `5000` entries，最高 cycle `particles` 有 `1808258` entries。
- 初次 wrapper 运行暴露旧二进制链接到旧 ROOT runtime 的 `libMinuit2.6.36.so` 路径问题。
- 在 `O2Physics/latest-master-o2` 下重新 configure/build 后，两个命令均为 `STATUS: PRIMARY_OK`。
- 沙箱内直接运行一键脚本会触发 `alienv` 的 `/dev/fd` 权限与 ROOT module loading 问题；同一命令在授权的非沙箱环境下成功。
- `Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh unexpected_arg` 按预期打印 `Usage` 并返回 exit code `2`。
- b1 分析输出:
  - `events_read=5000`
  - `selected_particles=389228`
  - `accepted_pairs=8770036`
  - `fit_summary` 行数 `14`，HBT xyz 成功 `14/14`，alpha 成功 `14/14`
- b8 分析输出:
  - `events_read=5000`
  - `selected_particles=188980`
  - `accepted_pairs=2086729`
  - `fit_summary` 行数 `14`，HBT xyz 成功 `14/14`，alpha 成功 `14/14`
- 输出 ROOT 检查:
  - 两个输出都返回 `RUNTIME_STATUS: PRIMARY_OK` 与 `STATUS: OK`。
  - 两个输出均包含 `Femto3D`、`R2Summary`、`fit_summary`、`analysis_statistics` 和 metadata 树。
- `ctest`:
  - `8` 项测试注册。
  - `5` 项直接通过。
  - `3` 项按 ROOT runtime guard 设计跳过。
  - `0` 项失败。

## 2026-04-13 当前 HEAD 对齐验证记录

- verification_status: `verified`

## 运行命令

- `ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build`
- `/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_eventgen_femto_3d.sh --config /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/blastwave_flat_trees.toml --output-root /tmp/eventgen_femto_3d_wrapper_smoke_20260413.root`

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
