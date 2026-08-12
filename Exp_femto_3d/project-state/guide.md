# Exp_femto_3d 项目说明

## 项目目的

`Exp_femto_3d` 是对原始 ROOT macro `3d_cf_from_exp.cpp` 的工程化重构版本。

这个项目的目标不是改变原有物理分析流程，而是把已经验证过的实验分析逻辑整理成一个常规的
`CMake + C++17` 项目，便于：

- 提高代码可读性和结构化程度
- 规范输入输出与参数传入方式
- 降低运行时无效输出
- 支持后续 agent 或人工对 build / fit 两阶段分别验证

## 当前分析流程

项目保留原始宏的两阶段工作流：

1. `build-cf`

   从实验 ROOT 文件中的 `THnSparseF` 读取 same-event 和 mixed-event 统计，
   按中心度、`mT`、`phi_pair - Psi_EP` 分区切片，生成：

   - `SE_raw3d`
   - `ME_raw3d`
   - `CF3D`
   - 1D projections

2. `fit`

   读取上一步生成的 CF ROOT 文件，基于显式 `SliceCatalog` 遍历切片，执行
   `diag` 或 `full` 的 3D Levy fit，并输出：

   - per-slice fit 结果
   - `FitCatalog`
   - `R2_vs_phi` 等 summary 图
   - `fit_summary.tsv`

## 代码结构

- `include/exp_femto_3d/`

  公共类型、配置接口、日志接口、workflow 接口。

- `src/`

  核心实现，包括：

  - TOML 配置解析
  - 日志输出
  - CF 构建
  - Slice catalog 读写
  - Levy fit 与 summary 输出

- `app/`

  CLI 主入口。

- `config/examples/`

  示例 TOML 配置。

- `tests/`

  配置解析测试、catalog roundtrip 测试、workflow smoke 测试、ROOT runtime probe。

- `legacy/`

  原始 ROOT macro 参考实现，不进入新项目主构建。

## 运行入口

项目主程序为：

```bash
./bin/exp_femto_3d
```

支持两个子命令：

```bash
./bin/exp_femto_3d build-cf --config <config.toml>
./bin/exp_femto_3d fit --config <config.toml>
```

可选 fit 覆盖参数：

```bash
--model full|diag
--input-cf-root <path>
```

## 配置方式

项目通过 TOML 读取配置，主要分为：

- `[input]`
- `[output]`
- `[build]`
- `[fit]`
- `[[bins.centrality]]`
- `[[bins.mt]]`
- `[[fit_selection.centrality]]`
- `[[fit_selection.mt]]`
- 可选 `[fit.parameters.<name>]`

示例配置可参考：

- `config/examples/exp_femto_3d.example.toml`
- `config/pbpb_build_and_fit.toml`

关键配置语义包括：

- `[fit.parameters.<name>]`

  控制主 3D Levy fit 参数的初值和上下限。支持的参数名为：

  - `norm`
  - `lambda`
  - `rout2`
  - `rside2`
  - `rlong2`
  - `routside2`
  - `routlong2`
  - `rsidelong2`
  - `alpha`
  - `baseline_q2`

  每个子表可写 `initial`、`min`、`max`；省略字段保持旧默认值。
  `lambda` 和 `alpha` 额外支持 `fixed_value`。`lambda` 覆盖要求
  `use_core_halo_lambda = true`，`baseline_q2` 覆盖要求
  `use_q2_baseline = true`。

- `[build].progress` / `[fit].progress`

  控制进度条输出，支持：

  - `true`
  - `false`
  - `"auto"`

  其中 `"auto"` 表示只有在 `stderr` 连接到 TTY 时才显示进度条。
  当前进度条会显示 stage label、百分比、活动帧和 ETA；启用时还有 1 秒心跳，
  长时间 ROOT 操作没有完成新 slice 时也会继续刷新运行提示。

- `[build.rebin.mt]` / `[build.rebin.phi]`

  控制 build 阶段是否对 mT / phi 轴做可配置 rebin。`enabled = false`
  表示按输入 sparse 的原始 normal bins 切片；`enabled = true` 时必须
  选择 `factor` 或 `[[bins.mt]]` / `[[bins.phi]]` 的显式范围模式之一。
  factor 模式只合并连续 normal bins，不会对最终 `TH3D` 调用 ROOT 的
  `Rebin`；显式范围模式要求区间与输入轴边界对齐。phi 侧还会把
  `build.rebin.phi` 是否启用写进 `SliceCatalog`，方便一眼看出当前
  build 是否在用 rebin。

  - `enabled = false`：不得再写 `factor`、`min/max` 或该轴的
    `[[bins.*]]`；build 使用输入轴的原生 normal bins。
  - `enabled = true`：必须恰好选择 `factor` 或 `[[bins.*]]` ranges。
  - factor 的可选 `min/max` 必须成对出现并对齐输入轴边界；窗口内 bin 数
    必须能被 factor 整除。
  - ranges 可以重叠以表达积分区间，但不能完全重复；每条边界都必须与输入轴
    边界对齐。
  - 对称 phi 映射启用时，native、factor 和 ranges 的任何非积分区间都不能
    跨过 `pi/2` 接缝。
  - 兼容旧配置：缺少新表时 phi 使用 native bins，已有 `[[bins.mt]]` 继续按
    legacy ranges 解释。

- `[fit].map_pair_phi_to_symmetric_range`

  这是一个可选覆盖项。

  - 不写时：`fit` 默认跟随输入 CF 文件在 `SliceCatalog` 中记录的 build 阶段
    phi 映射语义
  - 显式写为 `true` 或 `false` 时：`fit` 会基于 `raw_phi_*` 重新解释切片坐标，
    从而在不重建 CF 文件的情况下切换 `R2_vs_phi` 等 summary 的横轴语义

- `[build].split_mixed_event_by_phi`

  这是 build-cf 阶段的 ME 分母切分开关，默认值为 `false`。

  - `false`：保留原有行为，每个中心度/`mT` 组只投影一个覆盖全部 phi bin
    的 mixed-event 分母，供 phi 积分切片和逐 phi SE 切片共同使用
  - `true`：逐 phi SE 切片使用相同 phi bin 的 mixed-event 分母；phi 积分切片
    仍然使用完整 phi 范围
  - 该开关不改变 `raw_phi_*` / `display_phi_*` 的坐标映射语义；
    `map_pair_phi_to_symmetric_range` 只控制展示和 fit summary 的 phi 坐标解释

- `[build].split_mixed_event_by_qn`

  这是 build-cf 阶段的 ME qn 分母切分开关，默认值为 `false`。

  - `false`：保持旧行为，`split_same_event_by_qn = true` 写出的
    qn-specific SAME 切片仍使用 qn-integrated ME 分母
  - `true`：要求 `split_same_event_by_qn = true`，并让 ME 分母跟随当前
    qn-specific SAME 切片的 qn 区间；qn-all 切片仍使用完整 qn 轴
  - 该开关只改变 build-cf 写出的 raw ME / CF 内容，不改变 qn label、
    `group_id` 或 fit selection 规则

## 输出约定

### build-cf 输出

CF ROOT 文件中使用显式目录结构，而不是依赖 histogram 名称反解析 metadata：

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw3d`
- `slices/<slice_id>/ME_raw3d`
- `slices/<slice_id>/CF3D`
- `slices/<slice_id>/CF3D_ProjX/Y/Z`

此外，`SliceCatalog` 现在不仅保存 `raw_phi_*` 与 `display_phi_*`，还会保存
本次 build 是否使用对称区间映射的文件级 metadata
`build_uses_symmetric_phi_range`，以及本次 build 是否让 ME 分母跟随 SE phi
切片的 `split_mixed_event_by_phi`、是否让 ME 分母跟随 qn 切片的
`split_mixed_event_by_qn`。每个条目还保存
`mt_rebin_enabled/mt_rebin_mode` 与
`phi_rebin_enabled/phi_rebin_mode`。旧 catalog 缺少这些分支时按 legacy mT
ranges、native phi、integrated phi ME 和 integrated qn ME 读取。

### fit 输出

- `meta/FitCatalog`
- `fits/<slice_id>/...`
- `summary/R2_vs_phi/...`
- `fit_summary.tsv`

`FitCatalog` 会额外记录 `fit_uses_symmetric_phi_range`，明确当前 fit summary
到底使用的是原始 `[0, pi]` 坐标，还是对称 `[-pi/2, pi/2]` 坐标；它和 TSV
也会传递 build rebin mode 与最终 raw/display phi 范围。native/factor 模式下
省略 `[[fit_selection.mt]]` 表示选择 catalog 中的全部 mT groups，显式区间则
在读取 catalog 后校验，校验失败不会先清空已有 fit ROOT 输出。

## 构建与验证要求

### 构建要求

项目依赖：

- ALICE/O2 环境中的 ROOT
- Homebrew 安装的 `toml++`

### ROOT 调用要求

ROOT 相关命令必须在完整进入的 O2Physics 环境中执行。推荐方式：

```bash
alienv setenv O2Physics/latest-master-o2 -c sh -lc '
  cd /Users/allenzhou/Research_software/Code_base/Exp_femto_3d/build &&
  ctest --output-on-failure
'
```

### 已确认的环境注意事项

对本项目而言，agent 在沙箱里执行 `alienv` 时，可能出现：

- `/dev/fd/... Operation not permitted`

这会导致 ROOT/O2 环境只完成部分初始化，进一步诱发：

- PCM / dictionary 报错
- `THnSparseF::Projection(...)` 崩溃

这类报错在当前项目中应优先判断为环境进入失败，而不是代码逻辑错误。

详细诊断请见：

- `ROOT_RUNTIME_AGENT_NOTE.md`

## 当前状态

目前已经完成：

- 从 ROOT macro 到 CMake 项目的主流程迁移
- CLI + TOML 配置化
- SliceCatalog / FitCatalog 显式 metadata 设计
- build-side phi 映射状态持久化到 `SliceCatalog`
- fit 侧可跟随输入 CF metadata，也可显式覆盖 phi 映射语义
- build 侧可通过 `build.split_mixed_event_by_phi` 选择保留积分 ME 分母或让
  ME 分母跟随 SE phi 切片
- build 侧支持 mT/phi 的 native、factor、ranges/legacy-ranges sparse 轴分组，
  并在启动日志与 CLI 摘要打印开关、mode 和输入/输出 bin 数
- progress mode 配置支持 `true` / `false` / `"auto"`，进度条显示活动帧和 ETA
- 本地 unit / integration smoke tests；`2026-08-12` O2Physics ROOT executor
  完整 `ctest --test-dir build --output-on-failure` 六项全通过
- ROOT 调用失败模式诊断文档化
- `2026-04-19` 非沙箱 O2Physics `ctest` 三项全通过

当前仍建议继续完成的工作：

- 用真实实验数据重新做一次与 legacy macro 的数值等价验证
- 在真实数据回归里同时覆盖“跟随输入映射”和“fit 显式覆盖映射”两种语义
- 在真实数据回归里同时覆盖
  `build.split_mixed_event_by_phi = false` 和 `true` 两种 ME 分母语义
- 在真实数据上比较 native、factor 与显式 ranges rebin 的积分守恒和 fit 输出
- 继续把沙箱里的 `alienv` 失败视为环境进入问题，除非非沙箱复现出相同症状
