# Exp_femto_1d Legacy 宏工程化重构 + 1D CATS 拟合功能实施计划

## Summary

目标是在 `/Users/allenzhou/Research_software/Code_base/Exp_femto_1d` 内，把 `legacy/get_cf_from_exp.cpp` 重构为一个可维护、可测试、可配置的完整 `C++17 + CMake` 项目，并新增一个基于 CATS 的 1D 实验相关函数拟合流程。新项目的外层工程组织、CLI、TOML 配置、catalog 化输出、测试组织和 `project-state/` 协调方式，整体参考 `/Users/allenzhou/Research_software/Code_base/Exp_femto_3d`；但 1D 项目必须保留自身的输入 ROOT 结构、切片方式和拟合物理模型，不能机械照搬 3D 实现。

本计划已经锁定以下设计决策，执行时不要再改：
- 首版只支持 like-sign `pi-pi` 的 1D CATS 拟合
- 拟合模型固定为 `baseline polynomial × CATS correlation`
- 输出契约采用纯结构化新格式，不保留 legacy 平铺对象名兼容层
- event-plane 切片固定为 `MinBias / InPlane / OutOfPlane`
- CATS 接口只依赖已安装的公开头和 `libCATS.dylib`，不依赖 CATS 源树内未安装的私有 helper

## Environment And Constraints

- 工作目录：`/Users/allenzhou/Research_software/Code_base`
- 目标项目目录：`/Users/allenzhou/Research_software/Code_base/Exp_femto_1d`
- 3D 参考项目：
  - `/Users/allenzhou/Research_software/Code_base/Exp_femto_3d`
- CATS 参考原型：
  - `/Users/allenzhou/Research_software/CATS/work/fit_pipi.cpp`
- CATS 安装位置：
  - include：`/Users/allenzhou/Research_software/CATS/install/include`
  - lib：`/Users/allenzhou/Research_software/CATS/install/lib/libCATS.dylib`
  - helper：`/Users/allenzhou/Research_software/CATS/install/bin/cats-config`
- `cats-config --incdir` 返回：
  - `/Users/allenzhou/Research_software/CATS/install/include`
- `cats-config --libs` 返回：
  - `/Users/allenzhou/Research_software/CATS/install/lib -lCATS`
- `libCATS.dylib` 动态依赖包含：
  - GSL
  - ROOT 运行库
  - libomp
- ROOT 仍应从 ALICE/O2Physics 环境提供，沿用 `Exp_femto_3d` 的环境进入经验
- 仓库约束必须遵守：
  - 功能性代码块前必须有简短注释
  - 只要修改算法功能、配置契约、输出 schema、metadata、测试策略或用户运行方式，必须同步 `project-state/`
- 当前 `Exp_femto_1d/CMakeLists.txt` 明显是误拷贝的 3D 文件，首步必须纠正
- 当前 `Exp_femto_1d` 尚无 `project-state/`，本任务必须补建并采用

## Current-State Facts To Preserve

### Legacy 宏当前行为

`/Users/allenzhou/Research_software/Code_base/Exp_femto_1d/legacy/get_cf_from_exp.cpp` 当前核心流程是：

- 从同一个输入 ROOT 文件中读取：
  - `task_name/same_event_subtask/relPairkstarmTMultMultPercentileQn`
  - `task_name/mixed_event_subtask/relPairkstarmTMultMultPercentileQn`
- `THnSparseF` 轴语义是：
  - axis 0: `k*`
  - axis 1: `mT`
  - axis 2: `centrality`
  - axis 3: `event-plane`
- 对每个 `centrality × mT`：
  - 先构建 mixed-event 的 `Min bias EP`
  - 再构建 same-event 的 `Min bias EP / In_plane / Out_of_plane`
  - 用 `normrange` 对 SE/ME 做归一化
  - 在 `kstarRange` 上重采样/截断，生成 `CF_reranged_*`
- 当前 event-plane 定义固定为：
  - `In_plane = [0, pi/4] U [3pi/4, pi]`
  - `Out_of_plane = [pi/4, 3pi/4]`
  - mixed-event 始终使用 `[0, pi]`
- 当前写出的对象本质上包括：
  - `SE_Minbias_EP`
  - `ME_Minbias_EP`
  - `CF_reranged_Minbias_EP`
  - `CF_reranged_In_plane`
  - `CF_reranged_Out_of_plane`
  - 一个叠图 canvas

### 3D 参考项目当前可复用的工程思想

`Exp_femto_3d` 已经提供了可直接复用的工程模式：

- `include/src/app/tests/config/examples/legacy`
- `bin/` 作为可执行输出目录
- `build-cf` / `fit` 双阶段 CLI
- TOML 驱动
- `meta/SliceCatalog` 和 `meta/FitCatalog`
- 显式目录化 ROOT 输出
- config parse test / catalog roundtrip / workflow smoke / ROOT runtime guard
- `project-state/` 作为长期协调 ledger

执行时应最大化复用这些工程模式，但不要照搬 3D 的物理模型、phi 映射逻辑或 3D histogram 结构。

### CATS 原型当前可复用的物理/实现思想

`fit_pipi.cpp` 当前给出的核心思想是：

- 用一个预先配置好的 `CATS` 对象表示理论相关函数
- 用 `TF1` 桥接 ROOT fit 与 CATS 求值
- 拟合时把 source size 作为 CATS 的动态参数
- 拟合函数形态采用：
  - `baseline(k) × CATS(k; R)`
- 当前 baseline 原型为四阶内截断的多项式：
  - `p0 * (1 + p1*k + p2*k^2 + p3*k^3 + p4*k^4)`
- 当前原型拟合参数设置为：
  - `p0` 浮动 `[0.9, 1.1]`
  - `p1` 浮动 `[-0.2, 0.2]`
  - `p2` 浮动 `[-5, 5]`
  - `p3 = 0`
  - `p4 = 0`
  - `source_size` 浮动 `[3, 16]`
- 当前原型把输入直方图从 `GeV` 轴转换到 `MeV` 轴再与 CATS 连接

这些思想要保留，但要重写为工程化模块，不要直接复制原型程序结构。

## Final Project Layout

执行完成后，`Exp_femto_1d` 应至少包含：

- `CMakeLists.txt`
- `README.md`
- `ROOT_RUNTIME_AGENT_NOTE.md`
- `app/main.cpp`
- `include/exp_femto_1d/Types.h`
- `include/exp_femto_1d/Config.h`
- `include/exp_femto_1d/Logging.h`
- `include/exp_femto_1d/Workflow.h`
- `include/exp_femto_1d/CatsModel.h`
- `src/Config.cpp`
- `src/Logging.cpp`
- `src/Workflow.cpp`
- `src/CatsModel.cpp`
- `config/examples/exp_femto_1d.example.toml`
- `config/examples/minimal.toml`
- `config/pbpb_build_and_fit.toml`
- `tests/ConfigParseValidationTest.cpp`
- `tests/SliceCatalogRoundTripTest.cpp`
- `tests/WorkflowSmokeTest.cpp`
- `tests/CatsFitSmokeTest.cpp`
- `tests/RootRuntimeProbe.cpp`
- `tests/run_root_guard.sh`
- `legacy/README.md`
- `legacy/get_cf_from_exp.cpp`
- `project-state/guide.md`
- `project-state/current-status.md`
- `project-state/decisions.md`
- `project-state/tests.md`
- `project-state/issues.md`
- `project-state/work-items.md`
- `project-state/handoff.md`
- `project-state/changelog.md`

如果实现中发现某个辅助类还需要额外 `src/BuildCf.cpp` 或 `src/FitWorkflow.cpp` 之类拆分，可以加，但对外接口和职责边界必须保持本计划定义。

## Public CLI Contract

可执行文件固定为：

```bash
./bin/exp_femto_1d
```

支持两个子命令：

```bash
./bin/exp_femto_1d build-cf --config <config.toml>
./bin/exp_femto_1d fit --config <config.toml>
./bin/exp_femto_1d fit --config <config.toml> --input-cf-root <path>
```

不支持 `--model`。v1 只有一个固定模型：`cats-pipi-poly`。

CLI 行为固定为：

- `build-cf`
  - 读取实验 `THnSparseF`
  - 构建并写出结构化 1D CF ROOT 文件
  - stdout 打印摘要：
    - `stored_slices`
    - `skipped_zero_mixed_event_groups`
    - `skipped_zero_same_event_slices`
- `fit`
  - 读取 `build-cf` 输出的 `SliceCatalog`
  - 对选中的 slice 执行 CATS 拟合
  - 写出 `FitCatalog`、per-slice objects 和 summary
  - stdout 打印摘要：
    - `catalog_slices`
    - `selected_slices`
    - `fitted_slices`
    - `skipped_missing_objects`
    - `skipped_failed_fits`

出错时：
- stderr 打印 `[error] <message>`
- 随后打印 usage
- 返回非零状态

## Configuration Schema

### Top-Level Tables

TOML 结构固定为：

- `[input]`
- `[output]`
- `[build]`
- `[fit]`
- `[[bins.centrality]]`
- `[[bins.mt]]`
- `[[fit_selection.centrality]]`
- `[[fit_selection.mt]]`

### `[input]`

字段：

- `input_root`
  - 实验输入 ROOT 文件路径
  - 必填
- `task_name`
  - 例如 `femto-dream-pair-task-track-track`
  - 必填
- `same_event_subtask`
  - 例如 `SameEvent_EP`
  - 必填
- `mixed_event_subtask`
  - 例如 `MixedEvent_EP`
  - 必填
- `sparse_object_name`
  - 默认推荐：`relPairkstarmTMultMultPercentileQn`
  - 必填

### `[output]`

字段：

- `output_directory`
  - 必填
- `cf_root_name`
  - 默认 `cf_output.root`
- `fit_root_name`
  - 默认 `fit_output.root`
- `fit_summary_name`
  - 默认 `fit_summary.tsv`
- `log_level`
  - `debug|info|warn|error`
  - 默认 `info`

扩展名行为固定：
- `cf_root_name` 缺 `.root` 时自动补齐
- `fit_root_name` 缺 `.root` 时自动补齐
- `fit_summary_name` 缺 `.tsv` 时自动补齐

### `[build]`

字段：

- `norm_low`
  - 默认 `0.5`
- `norm_high`
  - 默认 `0.8`
- `kstar_min`
  - 默认 `0.0`
- `kstar_max`
  - 默认 `0.8`
- `reopen_output_file_per_slice`
  - 默认 `true`
- `progress`
  - `true|false|"auto"|"enabled"|"disabled"`
  - 默认 `"auto"`

行为固定：
- `norm_low < norm_high`
- `kstar_min < kstar_max`
- `fit_selection` 若省略，默认等于 build bins
- 不支持在配置中改 event-plane 区间

### `[fit]`

字段：

- `fit_kstar_max`
  - 默认 `0.25`
  - 单位 `GeV/c`
- `use_coulomb`
  - 默认 `true`
- `reopen_output_file_per_slice`
  - 默认 `true`
- `progress`
  - 默认 `"auto"`

参数初始化与约束字段固定加入：

- `baseline_p0_init`
  - 默认 `1.0`
- `baseline_p0_min`
  - 默认 `0.9`
- `baseline_p0_max`
  - 默认 `1.1`

- `baseline_p1_init`
  - 默认 `0.0`
- `baseline_p1_min`
  - 默认 `-0.2`
- `baseline_p1_max`
  - 默认 `0.2`

- `baseline_p2_init`
  - 默认 `1.0e-5`
- `baseline_p2_min`
  - 默认 `-5.0`
- `baseline_p2_max`
  - 默认 `5.0`

- `baseline_p3_fixed`
  - 默认 `true`
- `baseline_p3_value`
  - 默认 `0.0`

- `baseline_p4_fixed`
  - 默认 `true`
- `baseline_p4_value`
  - 默认 `0.0`

- `source_size_init`
  - 默认 `6.0`
- `source_size_min`
  - 默认 `3.0`
- `source_size_max`
  - 默认 `16.0`

额外固定字段：

- `cats_num_mom_bins`
  - 默认 `250`
- `cats_kmin_mev`
  - 默认 `0`
- `cats_kmax_mev`
  - 默认 `250`

v1 不开放 source type 切换，不开放 unlike-sign，不开放 baseline 模型切换。

### `[[bins.centrality]]` and `[[bins.mt]]`

每个 bin 支持：

- `min`
- `max`
- `label` 可选

约束固定：
- `max > min`
- 不允许 exact duplicate
- 允许 overlap，和 3D 项目保持一致
- `fit_selection` 必须精确匹配 build bins 中某个 bin，不允许局部截断

## Public Types And Interfaces

### `Types.h`

需要定义并稳定以下类型：

- `enum class LogLevel`
- `enum class ProgressMode`
- `enum class RegionKind`
  - `kMinBias`
  - `kInPlane`
  - `kOutOfPlane`
- `struct RangeBin`
- `struct InputConfig`
- `struct OutputConfig`
- `struct BuildCfConfig`
- `struct FitConfig`
- `struct ApplicationConfig`
- `struct SliceCatalogEntry`
- `struct PiPiFitResult`
- `struct BuildCfRunStatistics`
- `struct FitRunStatistics`
- `class ConfigError`

### `SliceCatalogEntry`

字段至少包括：

- `slice_id`
- `group_id`
- `slice_directory`
- `se_object_path`
- `me_object_path`
- `cf_object_path`
- `centrality_index`
- `mt_index`
- `region_index`
- `cent_low`
- `cent_high`
- `mt_low`
- `mt_high`
- `region_name`
- `region_kind`
- `ep_low_1`
- `ep_high_1`
- `ep_low_2`
- `ep_high_2`
- `has_second_interval`
- `norm_low`
- `norm_high`
- `kstar_min`
- `kstar_max`

### `PiPiFitResult`

字段至少包括：

- `slice_id`
- `group_id`
- `slice_directory`
- `centrality_index`
- `mt_index`
- `region_index`
- `cent_low`
- `cent_high`
- `mt_low`
- `mt_high`
- `region_name`
- `fit_kstar_max`
- `uses_coulomb`
- `baseline_p0`
- `baseline_p0_err`
- `baseline_p1`
- `baseline_p1_err`
- `baseline_p2`
- `baseline_p2_err`
- `baseline_p3`
- `baseline_p3_err`
- `baseline_p4`
- `baseline_p4_err`
- `source_size`
- `source_size_err`
- `chi2`
- `ndf`
- `fit_statistic`
- `edm`
- `status`
- `minuit_istat`

### `Workflow.h`

公开接口固定为：

- `ApplicationConfig LoadApplicationConfig(const std::string& path);`
- `void ValidateApplicationConfig(ApplicationConfig& config);`
- `BuildCfRunStatistics RunBuildCf(const ApplicationConfig& config, const Logger& logger);`
- `FitRunStatistics RunFit(const ApplicationConfig& config, const Logger& logger, std::optional<std::string> input_cf_root_path = std::nullopt);`
- `std::vector<SliceCatalogEntry> LoadSliceCatalog(const std::string& cf_root_path);`

### `CatsModel.h`

公开接口固定为窄模型接口，不暴露 CATS 细节到 workflow 外：

- `class PiPiCatsModel`
- 构造时接受：
  - `num_mom_bins`
  - `k_min_mev`
  - `k_max_mev`
  - `use_coulomb`
- 提供：
  - `double Evaluate(double kstar_mev, double source_size_fm);`
  - `TF1* BuildFitFunction(const std::string& name, const FitConfig& config);`
  - 供拟合桥接的静态/私有回调
- 该类内部维护 CATS 对象和 source 更新逻辑

## Build-CF Workflow Detail

### Stage Inputs

对每个配置文件，`build-cf` 固定读取：

- ROOT file：`input.input_root`
- SE sparse path：
  - `input.task_name + "/" + input.same_event_subtask + "/" + input.sparse_object_name`
- ME sparse path：
  - `input.task_name + "/" + input.mixed_event_subtask + "/" + input.sparse_object_name`

### Region Definitions

固定 region catalog：

- `MinBias`
  - `[0, pi]`
  - 单区间
- `InPlane`
  - `[0, pi/4]`
  - `[3pi/4, pi]`
  - 双区间
- `OutOfPlane`
  - `[pi/4, 3pi/4]`
  - 单区间

规则固定：
- ME 只构建 `MinBias`
- SE 构建三类 region
- 最终 slice 为：
  - `MinBias`
  - `InPlane`
  - `OutOfPlane`

### Projection Logic

对每个 `centrality × mT` group：

1. 读取 ME sparse
2. 在 centrality 和 mT 上设 range
3. 在 EP 上设 `MinBias`
4. 投影出 `ME_raw1d`
5. 若 `ME_raw1d` 归一化区间积分为 0，整个 group 跳过，并累计 `skipped_zero_mixed_event_groups`
6. 读取 SE sparse
7. 先构建 `SE_raw1d` 的 `MinBias` 版本
8. 再构建 `InPlane` 和 `OutOfPlane`
9. 若某个 SE slice 归一化区间积分为 0，则仅该 slice 跳过，并累计 `skipped_zero_same_event_slices`

### CF Construction Logic

对每个有效 slice：

1. clone `SE_raw1d` 与 group-level `ME_raw1d`
2. 开启 `Sumw2`
3. 在 `[norm_low, norm_high]` 上计算
   - `factor = integral(ME) / integral(SE)`
4. 用 constant TF1 或等效缩放对 ME 做 renormalize
5. 计算 `CF = SE / ME_norm`
6. 在 `[kstar_min, kstar_max]` 上重建输出 histogram，确保 bin 内容与误差被复制到新范围
7. 输出命名固定为：
   - `SE_raw1d`
   - `ME_raw1d`
   - `CF1D`

### Build Output ROOT Contract

ROOT 文件内结构固定为：

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw1d`
- `slices/<slice_id>/ME_raw1d`
- `slices/<slice_id>/CF1D`

不写 legacy 风格 `cent=..._CF_reranged_...` 平铺对象。

`SliceCatalog` 必须能让 `fit` 在不猜 histogram 名的情况下遍历全部 slice。

## Fit Workflow Detail

### Fit Input

`fit` 读取：

- 默认：
  - `output.output_directory/output.cf_root_name`
- 若传入 `--input-cf-root`：
  - 优先使用 CLI 覆盖

然后读取：

- `meta/SliceCatalog`
- 每个选中 slice 的 `cf_object_path`

### Fit Selection

- 只对匹配 `fit_selection.centrality × fit_selection.mt` 的 group 执行
- 对每个匹配 group，三个 region slice 全部拟合
- 不提供按 region 过滤的配置

### Unit Handling

拟合接口的单位规则固定：

- ROOT 存储和用户配置中的 `k*` 全部使用 `GeV/c`
- `PiPiCatsModel` 内部在调用 CATS 前把 `GeV/c` 转为 `MeV/c`
- `fit_kstar_max` 配置也用 `GeV/c`
- `cats_kmin_mev/cats_kmax_mev` 明确使用 `MeV/c`

### PiPiCatsModel Design

内部实现固定为：

- CATS object：
  - like-sign pion pair
  - `pdg = 211, 211`
  - `Q1Q2 = +1`
  - reduced mass 按原型处理
- source：
  - analytic Gaussian source
- notifications：
  - error level
- 每次求值时：
  - 更新 source size
  - 执行 `KillTheCat()`
  - 返回 `EvalCorrFun(kstar_mev)`

拟合桥接函数固定为：

- 参数数组顺序：
  - `0 = p0`
  - `1 = p1`
  - `2 = p2`
  - `3 = p3`
  - `4 = p4`
  - `5 = source_size`
- 返回：
  - `baseline(k) * cats_cf(k_mev, source_size)`

### ROOT Fit Procedure

对每个选中 slice：

1. 读取 `CF1D`
2. clone 为 `DataCF`
3. 构建 `TF1 FitFunction`
4. 设定参数初值、limits 和 fixed 状态
5. 在 `[0, fit.fit_kstar_max]` 上执行 fit
6. 抽取所有参数与误差
7. 构建：
   - `Baseline`
   - `PureFemto`
8. 构建 `FitCanvas`，叠加：
   - data
   - fit function
   - baseline
   - pure femto
9. 写入 fit ROOT 文件
10. 把结果写入 `PiPiFitResult`

### Failure Handling

- 缺少 `CF1D`：
  - 跳过该 slice
  - `skipped_missing_objects += 1`
- fit 执行失败但对象存在：
  - 仍写一条 `FitCatalog` 记录
  - 参数允许为 NaN
  - `status/minuit_istat` 保留真实失败值
  - `skipped_failed_fits += 1`
- 不因单 slice 失败中断整个 fit 阶段

## Fit Output ROOT And TSV Contract

### ROOT Contract

ROOT 文件结构固定为：

- `meta/FitCatalog`
- `fits/<slice_id>/DataCF`
- `fits/<slice_id>/FitFunction`
- `fits/<slice_id>/Baseline`
- `fits/<slice_id>/PureFemto`
- `fits/<slice_id>/FitCanvas`
- `summary/by_region/<group_id>/SourceSize`
- `summary/by_region/<group_id>/Norm`
- `summary/by_region/<group_id>/P1`
- `summary/by_region/<group_id>/P2`

summary 约定：

- 每个 group 生成 4 个 summary histogram
- 横轴 3 个 bin，label 固定：
  - `MinBias`
  - `InPlane`
  - `OutOfPlane`
- `P3/P4` 默认固定为 0，不单独出 summary 图
- 如果后续想加，可在 v2 做，不在本轮范围内

### `FitCatalog` Contract

`meta/FitCatalog` 必须包含前文 `PiPiFitResult` 全部字段。

### `fit_summary.tsv` Contract

header 固定至少包含：

- `sliceId`
- `groupId`
- `centLow`
- `centHigh`
- `mTLow`
- `mTHigh`
- `regionName`
- `fitKstarMax`
- `usesCoulomb`
- `p0`
- `p0Err`
- `p1`
- `p1Err`
- `p2`
- `p2Err`
- `p3`
- `p3Err`
- `p4`
- `p4Err`
- `sourceSize`
- `sourceSizeErr`
- `chi2`
- `ndf`
- `fitStatistic`
- `edm`
- `status`
- `minuitIstat`
- `covarianceQuality`

`covarianceQuality` 复用 3D 项目的字符串映射风格。

## File-By-File Implementation Plan

### `CMakeLists.txt`

- 改正 project 名称为 `ExpFemto1D`
- 启用 testing
- C++17
- `compile_commands.json`
- runtime 输出到 `Exp_femto_1d/bin`
- `find_package(ROOT REQUIRED ...)`
- `find_package(tomlplusplus CONFIG REQUIRED ...)`
- 新增 CATS 发现逻辑：
  - 优先 `CATS_ROOT`
  - 回退 `/Users/allenzhou/Research_software/CATS/install`
  - 读取 include/lib
- 必要时显式链接：
  - `libCATS.dylib`
  - GSL
  - GSL CBLAS
- 添加 target：
  - `exp_femto_1d_core`
  - `exp_femto_1d`
  - `config_parse_validation_test`
  - `root_runtime_probe`
  - `slice_catalog_roundtrip_test`
  - `workflow_smoke_test`
  - `cats_fit_smoke_test`
- ROOT 相关测试继续走 `run_root_guard.sh`
- config-only test 直接进 `ctest`

### `app/main.cpp`

- 复用 3D 项目 CLI 风格
- 解析 `build-cf` / `fit`
- 支持 `--config`
- `fit` 支持 `--input-cf-root`
- stdout 输出简短统计
- error path 打印 usage

### `src/Config.cpp`

- 复用 3D 解析模式
- 增加 build/fit 新字段解析与默认值
- 加入字段级验证：
  - `fit_kstar_max > 0`
  - `source_size_min < source_size_init < source_size_max` 不强制 strict 包含，但要检查 init 在区间内
  - `baseline_p0_min < baseline_p0_max`
  - 其余类似
- `fit_selection` 缺省时回退 build bins
- example config 应由测试加载校验

### `src/Workflow.cpp`

实现这些私有 helper：

- path resolve
- directory ensure
- ROOT open/create
- `BuildSparseObjectPath`
- `BuildGroupId`
- `BuildSliceId`
- `GetOrCreateDirectoryPath`
- `WriteSliceCatalogTree`
- `ReadSliceCatalogTree`
- `WriteFitCatalogTree`
- `WriteFitResultsSummaryTsv`
- `WriteRegionSummaryHistograms`
- `BuildMinBiasProjection`
- `BuildEpSliceProjection`
- `BuildCf1D`
- `LoadStoredHistogram1D`
- `MatchSelectedBin`

保持 `RunBuildCf` / `RunFit` 为主入口，不把主逻辑塞进 `main.cpp`。

### `src/CatsModel.cpp`

- 实现 `PiPiCatsModel`
- 实现 Gaussian source helper
- 实现 CATS object 初始化
- 实现 ROOT `TF1` bridge
- 实现 `BaselinePolynomial`
- 实现 `BuildPureFemtoFunction`
- 实现单位转换

### `tests/ConfigParseValidationTest.cpp`

至少覆盖：

- 最小合法 config
- 扩展名补全
- progress mode 解析
- 非法 `fit_kstar_max`
- 非法 parameter limits
- duplicate bins
- invalid `fit_selection`
- example config 文件可解析

### `tests/SliceCatalogRoundTripTest.cpp`

用 toy `THnSparseF`：

- 创建 4 轴 sparse：`k* / mT / cent / EP`
- 填充足够统计
- 执行 `RunBuildCf`
- 检查：
  - 生成 3 个 slice
  - `SliceCatalog` 存在
  - `slice_directory` 在 `slices/`
  - region metadata 正确
  - `CF1D` path 正确
  - `norm/kstar` metadata 正确

### `tests/WorkflowSmokeTest.cpp`

- 构造 toy input
- 执行 `RunBuildCf`
- 检查统计：
  - 组数
  - slice 数
- 打开 ROOT 输出检查：
  - `SE_raw1d`
  - `ME_raw1d`
  - `CF1D`
- 至少验证一个 `CF1D` 非空且 bin 数合理

### `tests/CatsFitSmokeTest.cpp`

固定策略：

- 使用 `PiPiCatsModel` 生成一条理论曲线
- 再乘上已知 baseline
- 填充成带误差 toy histogram
- 用 `fit` 模型回归
- 验证：
  - `source_size` 回收到容差内
  - `p0` 回收到容差内
  - `status/minuit_istat` 为成功状态
- 该测试不依赖实验 sparse，只依赖 ROOT + CATS

### `README.md`

必须写清：

- 项目目标
- 与 legacy 的关系
- 目录结构
- 构建方法
- 运行方法
- 配置 schema 简介
- output contract
- CATS 环境要求
- ROOT/O2 运行注意事项

### `ROOT_RUNTIME_AGENT_NOTE.md`

复用 3D 项目精神，写清：

- ROOT/CATS 依赖完整 O2/ROOT 环境
- sandboxed `alienv` 失败不应直接判定为代码问题
- 优先在非沙箱 O2Physics 环境下做 authoritative 测试

## Example Configs To Create

### `config/examples/minimal.toml`

最小可运行 demo，包含：

- `[input]`
- `[output]`
- `[build]`
- `[fit]`
- 1 个 centrality bin
- 1 个 mT bin

### `config/examples/exp_femto_1d.example.toml`

完整注释版示例，包含全部可配置字段和默认语义说明。

### `config/pbpb_build_and_fit.toml`

基于 legacy 宏当前 PbPb 使用习惯给出一个真实示例，建议默认：

- bins 与 legacy 宏一致：
  - cent:
    - `0-10`
    - `10-30`
    - `30-50`
    - `50-80`
    - `80-100`
  - mT:
    - `0.244131-0.331059`
    - `0.331059-0.423792`
    - `0.423792-0.51923`
    - `0.51923-0.713863`
    - `0.244131-0.713863`
- `same_event_subtask = "SameEvent_EP"`
- `mixed_event_subtask = "MixedEvent_EP"`
- `sparse_object_name = "relPairkstarmTMultMultPercentileQn"`

不要把 legacy 宏里那些 dataset 切换硬编码继续带进新项目。

## Project-State Writeback Plan

本任务执行结束前必须建立并同步：

### `project-state/guide.md`

写清：

- 项目目的
- 两阶段 workflow
- 配置方式
- 输出契约
- CATS v1 范围限制
- 验证方式

### `project-state/current-status.md`

记录：

- 当前完成了工程化重构
- 当前 fit v1 仅支持 like-sign `pi-pi`
- 当前真实数据等价性是否已完成

### `project-state/decisions.md`

至少记录：

- DEC-001：采用 `project-state/`
- DEC-002：v1 仅支持 like-sign `pi-pi`
- DEC-003：v1 模型固定为 `baseline × CATS`
- DEC-004：不依赖未安装的 `CommonAnaFunctions`
- DEC-005：不保留 legacy 平铺输出名

### `project-state/tests.md`

记录：

- config parse 测试
- slice catalog roundtrip
- workflow smoke
- cats fit smoke
- 非沙箱 O2/ROOT 环境测试结果
- 真实数据回归结果

### `project-state/issues.md`

至少保留一个 issue：

- 真实数据 legacy 数值等价验证未闭环时，标为 open

### `project-state/work-items.md`

至少保留：

- real-data equivalence validation
- inspect fit stability on real slices
- optional future support for unlike-sign or source-model selection

### `project-state/handoff.md`

记录：
- 最新完成内容
- 当前残留风险
- 推荐下一个动作
- 推荐命令

### `project-state/changelog.md`

按日期记录关键变化。

## Validation And Acceptance Criteria

### Build Acceptance

满足以下条件才算 `build-cf` 完成：

- `ctest` 中 config test 通过
- toy sparse 的 roundtrip test 通过
- workflow smoke test 通过
- 输出 ROOT 文件内存在：
  - `meta/SliceCatalog`
  - `slices/<slice_id>/SE_raw1d`
  - `slices/<slice_id>/ME_raw1d`
  - `slices/<slice_id>/CF1D`
- `SliceCatalog` 不依赖对象名反解 metadata

### Fit Acceptance

满足以下条件才算 `fit` 完成：

- `cats_fit_smoke_test` 通过
- `meta/FitCatalog` 存在
- per-slice fit objects 存在
- `fit_summary.tsv` 存在且 header 完整
- summary histograms 存在
- 至少一个真实 slice 的 fit canvas 人工可读且 baseline/pure femto 分离合理

### Full Task Acceptance

满足以下条件才算整个任务完成：

- `Exp_femto_1d` 成为独立 CMake 项目
- `legacy/get_cf_from_exp.cpp` 保留为参考实现，但不进入主构建
- README 可指导首次构建和运行
- 所有新增或变更的算法/配置/schema/test/user-run 语义已经同步到 `project-state/`
- 至少一次在非沙箱 ROOT/O2 环境中完成 `ctest --output-on-failure`
- 完成一轮真实数据对照：
  - 新 `build-cf` 对比 legacy 宏，至少在一个已知 good dataset 上确认 `SE/ME/CF` 一致到可接受浮点误差
  - 新 `fit` 在至少一个 `MinBias` 和一个 EP-differential slice 上给出物理上合理的拟合结果

## Recommended Execution Order

1. 修正 `CMakeLists.txt`，建立基础 target 和目录结构
2. 从 `Exp_femto_3d` 迁移 `Types/Config/Logging/main` 的工程壳
3. 实现 1D `SliceCatalog` 契约
4. 实现 `RunBuildCf`
5. 补齐 config test 与 slice catalog roundtrip test
6. 实现 `PiPiCatsModel`
7. 实现 `RunFit`
8. 补齐 `cats_fit_smoke_test`
9. 写 README 和 `ROOT_RUNTIME_AGENT_NOTE.md`
10. 建立并同步 `project-state/`
11. 在非沙箱 ROOT/O2 环境中运行 `ctest`
12. 做真实数据 legacy 等价回归并更新 `project-state/tests.md` 与 `issues.md`

## Explicit Non-Goals For This Pass

- unlike-sign `pi-pi`
- 多粒子体系抽象壳
- 可切换 source family
- legacy 输出对象命名兼容
- GUI、Python 绑定、notebook 支持
- 复杂全局拟合或 simultaneous fit
- systematics 合并与误差传播扩展

## Assumptions And Defaults

- `Exp_femto_1d` 的用户仍主要在 ALICE/O2Physics 环境中运行
- v1 的 baseline 自由度按 `fit_pipi.cpp` 原型固定，不在本轮做模型比较
- region summary 采用分类轴图，而不是 phi 连续图
- 如果真实数据 fit 表现表明 `p2` 范围过宽或不稳，允许在实现末尾微调默认 limits，但必须同步：
  - example config
  - README
  - `project-state/decisions.md`
  - `project-state/changelog.md`
- 若发现某些真实输入文件并非 `relPairkstarmTMultMultPercentileQn`，只通过 TOML 配置改路径，不再在代码里加 dataset-specific hardcode
