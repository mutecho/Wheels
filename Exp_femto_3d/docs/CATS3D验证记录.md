# CATS3D 接入验证记录（2026-10-19）

本记录区分接入测试、真实数据上的数值收敛和物理精度验证。所有 ROOT 构建、CTest、拟合、读取和绘图都通过 O2Physics guarded executor；环境为 `O2Physics/latest-master-o2`、ROOT 6.36.10。外部 CATS 的安装库来自本地 Debug 构建，未修改或重编译外部项目。

## 通用核与 16 GiB 内存预算更新

随后按用户要求移除了库二进制、自定义合同和内容指纹的缓存加载门槛。保留 CATS3D 原生网格与物理元数据检查、q 记录完整性、向量长度和有限非负值检查。新增 `kernel_cache_file`（只读、任意文件名的原生 ByQ 输入）和 `max_peak_gib=8`；当前实现详见 [核缓存与内存说明](CATS3D拟合与lambda扫描.md#核缓存与验证边界)。以下初版测试和真实拟合记录保留为当时的观测，不将旧的缓存指纹约束视为当前要求。

对原文件 `/Users/allenzhou/Research_software/CATS/CATS_3D/OutputFiles/dataKernel_qData_r31_coulomb_byQ.root`，使用其自身 20³ q 网格、31³ 空间网格只读加载，并计算一个高斯源关联函数值。该探针没有更改 cuttest 数据分箱，也没有生成全尺寸新核：

- 文件 754,996,901 字节；核数据 953,312,000 字节（0.88784 GiB）。
- 加载和一次源卷积共 6.67 秒；最大 RSS **2,259,468,288 字节（约 2.104 GiB）**，0 swaps。
- 该尺寸的峰值预算估算为 2.776 GiB；实际加载通过，不需要任何 ExpFemto 私有指纹。
- 当前 cuttest 的 30³ 拟合 q 网格对应核数据 **2.996 GiB**；磁盘 **约 2.37 GiB** 为旧核压缩比外推，峰值 **约 6.99 GiB** 为两份核加工作区估算，均不是全尺寸实测。

探针源码和日志位于 `build/cats3d_memory_probe/`。全部操作通过 O2Physics guarded executor，`PRIMARY_OK`。

本次缓存与内存改动的自动验证：

- 开启 CATS3D 的完整 CTest **13/13 通过**（364.05 秒），隔离的无 CATS / 无 CATS3D 构建 **13/13 通过**（49.27 秒）。日志为 `build/cats3d_existing_kernel_audit/portable_ctest.log` 和 `portable_no_cats_ctest.log`。
- 最后补充了显式核文件的写出前加载检查，并让 CATS3D 流程直接持有、及时释放 ROOT 读取的 CF/SE/ME，避免逐切片累计原始直方图。重新构建后，模型与完整工作流测试 **2/2 通过**（302.62 秒）；无 CATS 的配置/工作流测试 **1/1 通过**（2.03 秒）。对应 `final_memory_ctest.log`、`final_no_cats_test.log`，最终构建没有编译警告。
- 模型测试覆盖任意原生缓存文件名、无私有指纹输入、忽略旧私有标签、缓存开关数值一致、q 网格和物理失配拒绝、NaN 内容拒绝，以及分配前单核和峰值预算限制。缺少指定核文件时拒绝请求，不创建替代文件。
- 用真实 CF 和独立输出路径执行 CLI 保护检查：原有 20³ q 核与当前 30³ q 网格不匹配、指定文件缺失、`max_peak_gib=6.5` 低于约 6.99 GiB 的估算，这三种情况都在核分配和输出初始化前报错；正式 ROOT、TSV、报表和扫描文件四份哨兵内容均逐字节保持不变。证据为 `build/cats3d_existing_kernel_audit/protection/result.log`。

以上没有生成 30³×31³ 的全尺寸新核，也没有执行全规模 cuttest 拟合或扫描。8 GiB 是分配前的估算预算，不是操作系统强制的 RSS 上限；16 GiB 电脑应串行运行并为系统及其他应用留出内存。同一进程的各切片及扫描点复用核，扫描 ROOT 不嵌入整份核。

## 自动测试

- 最终主构建 CTest：**13/13 通过**，359.82 秒；其中 CATS3D 工作流测试 299.60 秒。包含独立配置种子与双向邻点重试、缓存写入检查和输出目录创建的最终实现。另有 CLI χ² 检查，4/4 切片有有效最小值和准确协方差。
- 隔离的无 CATS / 无 CATS3D 构建：**13/13 通过**，45.18 秒。可执行文件放在独立构建目录，未覆盖主构建。CATS3D 请求明确拒绝，旧 Bowler 工作流继续通过。
- 模型测试：各向异性及交叉项的无库仑解析闭合、与本地 CATS3D 直接库仑调用闭合、GeV/MeV 和 q/2 约定、λ=0/1、baseline、输入覆盖/插值/纯度换算、缓存开关及损坏/失配拒绝。
- 工作流测试：原始整数 SE/ME toy 的 diag/full PML 和 χ² 参数恢复与自由度；组覆盖和表误差保存；常数/倍率扫描、所有尝试和失败点；积分切片不进入组总和；`scan_only` 前后正式 ROOT、TSV 和报表逐字节一致；输出路径冲突在写入前拒绝。

## 真实输入和设置

输入：`/Users/allenzhou/ALICE/alidata/femtoep_res/OO/EP_dependence_CF_3d_25ae_pass2_sel8_6phi.root`。
选择一个组：centrality 0–10%，mT 0.4–0.6 GeV/c²，qn integrated；6 个 φ 切片和 1 个 φ 积分切片。原始 q 轴为 40³、每轴 [−0.2, 0.2] GeV/c；不重分箱。

λ=0.7 是演示固定输入，**不是 OO 标定结果**。高斯 α=2，库仑及 QS 开启，baseline 关闭，max_pw=16、subdivisions=1。扫描坐标为最终 λ=[0.65, 0.70]。本次范围受控，不能代表完整生产配置。

| 运行 | 几何 / q 盒 | 源网格 | MIGRAD 上限 | 有效最小值 / 准确误差 | 墙钟 / 最大 RSS |
|---|---|---|---|---|---|
| coarse | full，±0.034 GeV/c，6³ bins | 16³，±15 fm | 500 | 0/7，0/7 | 30.91 s / 384.9 MiB |
| coarse_refined | 同上 | 同上 | 5000 | 0/7，0/7 | 26.77 s / 385.8 MiB |
| diag_coarse | diag，±0.064 GeV/c，12³ bins | 16³，±15 fm | 5000 | 7/7，7/7 | 70.85 s / 537.1 MiB |
| diag_fine（不扫描） | 同上 | 32³，±15 fm | 5000 | 7/7，7/7 | 356.14 s / 850.2 MiB |

前两轮保留了全部不收敛记录及扫描尝试，公开目标函数和不可用误差为 NaN。提高调用上限并未改变这些失败；没有把它们当成成功的 full 真实拟合。diag 同时改变了几何和 q 范围，因此这些运行不能用于判断失败由哪一个因素单独导致。

## 空间网格检查

在每个粗网格有效拟合的固定参数处，用两个生产适配器重新计算最终模型；7 个切片的最大绝对差为 **1.00509×10⁻⁵**。两套网格分别重拟合后，三个 R² 元素的最大位移为细网格条件 HESSE 误差的 **0.0136487σ**（含积分切片）。逐切片结果见 `build/cats3d_validation/grid_comparison.csv`。

此检查只改变空间间距（16³→32³），保持 ±15 fm、q 分箱、λ、max_pw 和其他设置相同；没有验证空间截断范围或 partial-wave 收敛。两轮都使用外部固定 λ，不传播 λ 标定误差。有效最小值和准确协方差只说明数值状态，不能单独证明真实数据模型拟合质量。这组数据的 PML 目标函数/ndf 为积分切片约 1.419、各 φ 切片约 1.032–1.087；均未命中参数边界。

12³ q bins 对应的波函数核为粗网格 27 MiB、细网格 216 MiB，启动峰值估算分别约 81.06 MiB、648.50 MiB；进程 RSS 还包括 ROOT、输入直方图、库和其他对象。细网格计时包含首次波函数核生成；不能直接用它与带短扫描的粗网格总时长推算单次拟合加速比。

## 扫描隔离与图像核对

最终代码使用配置种子和双向相邻点种子，对同组执行独立 `scan_only`（`diag_coarse/configured_scan.toml`）：7 个切片×2 个坐标的 **14/14 点有效且误差有效**，全部 28 次实际尝试均保留。组总曲线逐点核对为恰好 6 个非积分 φ 成员的目标函数之和，积分切片单独保存。λ=0.65、0.70 的组 Δ目标函数分别为 0、约 3.6339；这是两点诊断，不能选作正式 λ 或置信区间。

扫描前后正式 `fit.root`、`summary.tsv`、`report.root` 的 SHA-256 全部一致。原始 CF 和外部 `libCATS` 的大小及 SHA-256 也再次核对一致。

分别用独立 ROOT 进程渲染并查看了模型投影、组扫描和总览图。视觉核对后，仅补充 CATS3D 总览标识、修正投影图例/标题并关闭 CF 投影的默认计数统计框；同组随后再次 7/7 拟合成功（`label_probe/`）。没有更改模型或扫描数值。图像为 `real_projection.png`、`real_lambda_scan.png`、`real_overview.png`。

结构化证据、输入指纹、精确配置和命令见 [CATS3D_validation_record.json](CATS3D_validation_record.json)。

## 可重复运行和产物

验证配置、日志、ROOT、TSV、缓存和图像均在 `build/cats3d_validation/` 下，未写入正式分析输出。示例运行：

```sh
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command '/usr/bin/time -l bin/exp_femto_3d fit --config build/cats3d_validation/diag_coarse/config.toml > build/cats3d_validation/diag_coarse/run.log 2>&1'
```

主测试为 `ctest --test-dir build --output-on-failure`；隔离构建配置为：

```sh
cmake -S . -B build/no_cats_validation \
  -DEXP_FEMTO_3D_ENABLE_CATS=OFF -DEXP_FEMTO_3D_ENABLE_CATS3D=OFF \
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/build/no_cats_validation/bin
cmake --build build/no_cats_validation --parallel 8
ctest --test-dir build/no_cats_validation --output-on-failure
```

上述构建与测试命令也须放在 guarded executor 中执行。使用目录 `build/no_cats_validation` 是为了保留主构建。

## 已知边界

默认 31³、±50 fm 的源网格尚不能视为生产精度已验证。此记录中的真实数据仅覆盖一个组、较小 q 区域和演示 λ；未执行全数据生产拟合。空间范围、间距和 partial-wave 都需针对实际半径检查。

本次沿用当前 CATS3D 核心；其库仑相移 TODO、LCMS/PRF 坐标近似和外部 λ 的物理标定仍是独立验证事项。接入闭合只说明适配器与该核心一致。
