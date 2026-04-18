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

示例配置可参考：

- `config/examples/exp_femto_3d.example.toml`
- `config/examples/pbpb_build_and_fit.toml`

## 输出约定

### build-cf 输出

CF ROOT 文件中使用显式目录结构，而不是依赖 histogram 名称反解析 metadata：

- `meta/SliceCatalog`
- `slices/<slice_id>/SE_raw3d`
- `slices/<slice_id>/ME_raw3d`
- `slices/<slice_id>/CF3D`
- `slices/<slice_id>/CF3D_ProjX/Y/Z`

### fit 输出

- `meta/FitCatalog`
- `fits/<slice_id>/...`
- `summary/R2_vs_phi/...`
- `fit_summary.tsv`

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
- 本地 unit / integration smoke tests
- ROOT 调用失败模式诊断文档化

当前仍建议继续完成的工作：

- 用真实实验数据重新做一次与 legacy macro 的数值等价验证
- 决定 ROOT runtime guard 在 agent 测试里是否长期保留
