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
    - 输出 ROOT 目录与统计树

- 物理与数值子模块
  - [EventPlane.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/EventPlane.cpp)
  - [Histogramming.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Histogramming.cpp)
  - [ProjectionFit.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/ProjectionFit.cpp)
  - [Source_extraction.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/Source_extraction.cpp)

- 入口层
  - [main.cpp](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/src/main.cpp)
  - [run_eventgen_femto_3d.sh](/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/run_eventgen_femto_3d.sh)
  - `main` 只做 CLI glue；shell wrapper 负责在 O2/ROOT 环境下透传参数运行二进制。

- 测试层
  - 配置解析：`ConfigParseValidationTest.cpp`
  - 输入适配：`InputAdapterTest.cpp`
  - legacy smoke：`LegacyWorkflowSmokeTest.cpp`
  - ROOT runtime guard：`RootRuntimeProbe.cpp` + `run_root_guard.sh`

## 使用方法

### 1. 配置与构建

推荐在 O2/ROOT 运行时中构建：

```bash
cmake -S /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d \
      -B /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build

cmake --build /Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/build
```

### 2. 运行入口

命令行接口：

```bash
eventgen_femto_3d --config <file.toml> \
  [--input-root <path>] \
  [--output-root <path>] \
  [--input-schema legacy_vector_tree|blastwave_flat_trees]
```

如果需要在本机 O2/ROOT 环境下直接运行，可使用 wrapper：

```bash
/Users/allenzhou/Research_software/Code_base/Eventgen_femto_3d/run_eventgen_femto_3d.sh \
  --config /path/to/config.toml
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

- CLI/TOML 解析
- blast-wave 输入聚合与 fail-fast
- legacy workflow smoke

`projection_fit.accept_hbt_central_value_only_for_summary` 用于控制：

- 当 HBT 半径中心值有效但误差不可用时，`R2Summary` 是否仍保留该点
- 若启用，summary 图会写入中心值且误差条记为 0
- 若关闭，summary 图会跳过该点，并在运行摘要与 `analysis_statistics` 中累计跳过计数

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
