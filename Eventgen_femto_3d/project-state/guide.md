# Eventgen_femto_3d 项目说明

> 说明：本文件基于当前代码结构与本轮实现结果整理，用于帮助后续维护者快速理解项目目的、架构与使用方式。若后续出现更高权威的人写设计文档，应以人写文档和代码行为为准。

## 项目目的

`Eventgen_femto_3d` 的目标是从事件生成器输出的 ROOT 数据中提取三维 femtoscopy 分析所需的源分布与拟合结果，并将结果写回新的 ROOT 输出文件。

本项目当前支持两类输入契约：

- `legacy_vector_tree`
  - 单棵 `events` 树，事件级信息与粒子级 `std::vector<>` 分支放在同一树内。
- `blastwave_flat_trees`
  - `events` 事件树加 `particles` 粒子树，粒子通过 `event_id` 与事件关联。

项目本轮的重点是：

- 让输入契约从分析流程中解耦
- 让运行入口改成 `CLI + TOML`
- 保留 legacy 工作流兼容，避免回归

## 项目架构

当前代码可按下面几层理解：

- 配置层
  - [Config.h](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/include/femto3d/Config.h)
  - [Config.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Config.cpp)
  - 负责 CLI 参数解析、TOML 解析、默认值、配置校验以及 CLI 覆盖 TOML。

- 配置数据模型
  - [AnalysisConfig.h](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/include/femto3d/AnalysisConfig.h)
  - 定义 `ApplicationConfig`、`AnalysisConfig`、输入 schema、event-plane、selection、histograms、projection_fit、bin 配置。

- 输入适配层
  - [InputReader.h](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/include/femto3d/InputReader.h)
  - [InputReader.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/InputReader.cpp)
  - 负责读取 ROOT 输入并归一化为 `EventData`。
  - `blastwave_flat_trees` 模式下会校验：
    - 必需树和分支存在
    - 事件树 `event_id` 连续
    - 粒子 `event_id` 不越界
    - 标量分支类型符合预期

- 分析主流程
  - [Workflow.h](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/include/femto3d/Workflow.h)
  - [Workflow.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Workflow.cpp)
  - 负责：
    - centrality / mT / phi 分 bin
    - event-plane 选择或回退
    - pair building
    - 源分布填充
    - projection fit
    - 输出 ROOT 目录、统计树、`R2Summary` 与 `epsf_vs_mt` 汇总图

- 物理与数值子模块
  - [EventPlane.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/EventPlane.cpp)
  - [Histogramming.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Histogramming.cpp)
  - [ProjectionFit.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/ProjectionFit.cpp)
  - [Source_extraction.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Source_extraction.cpp)

- 入口层
  - [main.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/main.cpp)
  - [run_eventgen_femto_3d.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_eventgen_femto_3d.sh)
  - [run_dense_mix_glauber_7phibin.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh)
  - [run_dense_mix_glauber_15phibin.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh)
  - `main` 只做 CLI glue。
  - `ProgressReporter` 负责 CLI 进度条，分析流程只通过 `AnalysisProgressSink` 上报完成事件数。
  - CMake 当前把主程序输出到 `bin/eventgen_femto_3d`。
  - shell wrapper 负责在 O2/ROOT 环境下透传参数运行 `bin/eventgen_femto_3d`。
  - dense-mix Glauber 一键脚本负责固定 config/input/output 路径并连续运行 b1、b3 与 b8 分析。

- 测试层
  - 纯逻辑/语义测试：
    - `DirectionalCovarianceFailureSemanticsTest.cpp`
    - `AlphaHbtErrorSemanticsTest.cpp`
    - `FitSummaryDirectionalErrorMaskTest.cpp`
    - `R2SummaryPolicyTest.cpp`
    - `ConfigParseValidationTest.cpp`
  - 运行时相关测试：
    - `InputAdapterTest.cpp`
    - `LegacyWorkflowSmokeTest.cpp`
    - `R2SummaryVisibilityTest.cpp`
  - ROOT runtime guard：
    - `RootRuntimeProbe.cpp`
    - `run_root_guard.sh`
  - 测试二进制当前输出到 `bin/tests/`。

## 使用方法

### 1. 配置与构建

推荐在 O2/ROOT 运行时中构建：

```bash
/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/scripts/cmake.sh
```

该脚本会从自身位置解析项目根目录，默认使用
`/Users/allenzhou/Research_software/Code_Base/Eventgen_femto_3d/build`，并将主程序输出到源码树
`bin/eventgen_femto_3d`。可用 `EVENTGEN_FEMTO_3D_BUILD_DIR` 覆盖 build tree，用
`EVENTGEN_FEMTO_3D_BUILD_JOBS` 覆盖并行数。默认走增量构建；需要全量清理时设置
`EVENTGEN_FEMTO_3D_CLEAN_FIRST=1`。

当当前 shell 可见 `root-config` 时，`cmake.sh` 会把当前 ROOT 的 `ROOT_DIR`
显式传给 CMake；若 CMake cache 仍指向旧 ROOT 包，脚本会清理 `ROOT_*` cache 并触发
relink，避免 `bin/eventgen_femto_3d` 链接旧 ROOT 而 wrapper 运行在新 ROOT 时产生
`TClassTable::Add` / `TCling::LoadPCM` 噪声。

构建完成后，主程序默认位于：

```bash
/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/bin/eventgen_femto_3d
```

### 2. 运行入口

命令行接口：

```bash
eventgen_femto_3d --config <file.toml> \
  [--input-root <path>] \
  [--output-root <path>] \
  [--input-schema legacy_vector_tree|blastwave_flat_trees] \
  [--progress|--no-progress]
```

进度条默认是 `auto` 模式：stderr 是 TTY 时显示。批处理日志中可用 `--no-progress` 关闭；需要强制显示时可用 `--progress`。

如果需要在本机 O2/ROOT 环境下直接运行，可使用 wrapper：

```bash
/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_eventgen_femto_3d.sh \
  --config /path/to/config.toml
```

该 wrapper 当前会进入 `alienv setenv O2Physics/latest-master-o2`，然后调用 `bin/eventgen_femto_3d`。

dense-mix Glauber 分析可用一键脚本运行；config、b1/b3/b8 输入 ROOT 与输出 ROOT 路径均在脚本内部固定：

```bash
/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_7phibin.sh

/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/scripts/run_dense_mix_glauber_15phibin.sh
```

15-phi-bin 脚本使用的配置文件名当前包含空格：

```bash
/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/bw_dense_mix_glauber_15phibin copy.toml
```

### 3. 示例配置

示例文件位于：

- [legacy_vector_tree.example.toml](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/examples/legacy_vector_tree.example.toml)
- [blastwave_flat_trees.example.toml](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/config/examples/blastwave_flat_trees.example.toml)

### 4. 测试

```bash
ctest --output-on-failure --test-dir /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build
```

当前自动化重点覆盖：

- 拟合语义与 `R2Summary` 策略
- 进度条格式与 `epsf` 提取
- CLI/TOML 解析
- blast-wave 输入聚合与 fail-fast
- legacy workflow smoke

`projection_fit.accept_hbt_central_value_only_for_summary` 用于控制：

- 当 HBT 半径中心值有效但误差不可用时，`R2Summary` 是否仍保留该点
- 若启用，summary 图会写入中心值且误差条记为 0
- 若关闭，summary 图会跳过该点，并在运行摘要与 `analysis_statistics` 中累计跳过计数

`projection_fit.alpha_min` 与 `projection_fit.alpha_max` 控制 simultaneous Levy
projection fit 中 alpha 参数传给 Minuit 的边界。两个字段都是可选配置；未设置时默认仍为：

```toml
[projection_fit]
alpha_min = 0.2
alpha_max = 2.0
```

配置校验要求 `0 < alpha_min < alpha_max`。默认 `alpha_max = 2.0` 对应
Levy-stable alpha 的物理上界；若为了诊断显式设置 `alpha_max > 2`，输出可以运行，
但该 alpha 不应再按物理 Levy-stable 指数直接解释。输出 ROOT 的
`projection_fit_metadata` 会记录实际使用的 `alpha_min` 与 `alpha_max`。

`epsf_vs_mt` 当前从每个 centrality 下已有的 `Rside2_vs_phi` summary 点提取：

```text
#epsilon_f = 2 R_{s,2}^{2} / R_{s,0}^{2}
```

输出对象位于 `R2Summary/<centrality>/epsf_vs_mt`，对应 canvas 为 `R2Summary/<centrality>/epsf_vs_mt_canvas`。

每个有可画 summary 点的 `R2Summary/<centrality>/<femto_mt>/` 还会写出
`source_parameters_overview_canvas`，用于按图版顺序汇总同一 centrality 与 femto
`m_T` 区间下的 source parameters。当前 pad 顺序为：

1. 区间信息与相对角说明
2. `alpha` vs `phi_pair - Psi_2`
3. `Rout2_vs_phi`
4. `Ros2_vs_phi`
5. `Rside2_vs_phi`
6. `Rol2_vs_phi`
7. `Rlong2_vs_phi`
8. `Rsl2_vs_phi`

`alpha` 点来自每个 `cent/femto_mt/phi` slice 的 `fit_products.alpha`，只作为
overview canvas 内部图元使用；当前没有新增独立 `alpha_vs_phi` ROOT key。

## 维护注意事项

- `participants` 树当前不属于 `Eventgen_femto_3d` 输入契约的一部分。
- 若修改输入 schema，请同步更新：
  - `AnalysisConfig`
  - `Config.cpp`
  - `InputReader.cpp`
  - example TOML
  - 对应测试
- 若仓库继续采用 `project-state/` 作为协调台账，后续状态变化任务结束前应至少同步：
  - `current-status.md`
  - `handoff.md`
  - `changelog.md`
