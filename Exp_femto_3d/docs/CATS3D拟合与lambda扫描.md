# CATS3D 高斯拟合与外部 λ 扫描

`fit.correlation_backend = "cats3d"` 使用本机 CATS3D 的三维源与波函数卷积。旧配置未指定此字段时，继续使用 `bowler`。两种后端共享输入 CF、切片选择和拟合统计量；CATS3D 不调用外部 `FitDataCATS3D` 的 FCN、背景预拟合或归一化窗口。

## 模型与单位

\[
C_{\rm core}(\mathbf q;M)=\int d^3r\,S(\mathbf r;M)|\psi_{\mathbf q}(\mathbf r)|^2,
\qquad
C_{\rm model}=N\,[1+\lambda_{\rm main}(C_{\rm core}-1)]B(\mathbf q).
\]

本地 CATS3D 在有限空间网格上对源权重归一化；该数值归一化不能替代网格截断和离散化收敛检查。QS 已在 CATS3D 波函数中，不再外乘 `1+exp(-E)` 或一维 Coulomb 表。

高斯源使用

\[
M=\begin{pmatrix}R_o^2&R_{os}^2&R_{ol}^2\\R_{os}^2&R_s^2&R_{sl}^2\\R_{ol}^2&R_{sl}^2&R_l^2\end{pmatrix},\qquad
S=\frac{\exp[-\mathbf r^TM^{-1}\mathbf r/4]}{(4\pi)^{3/2}\sqrt{\det M}}.
\]

矩阵元素单位为 fm²，交叉项是有符号矩阵元素。`diag` 令三个交叉项为零；`full` 拟合六个矩阵元素。要求严格正定，不能用奇异的半正定矩阵定义该源密度。α 自动固定为 2；显式冲突的 α 设置会报错。无库仑极限为 `1+exp[-qᵀMq/(ℏc)²]`。

输入三个 q 轴沿用当前分析的 GeV/c 单位，适配层乘 1000 后传给 CATS3D。库内部已包含 q/2，不再次除 2。仅使用普通 bins 中满足 `abs(q_i) <= fit_q_max` 的原有 bin 中心，不重分箱；当前核心只能准确表达均匀 q 轴，非均匀轴会报错。模型与曲线只覆盖实际计算区域。

λ 表查询量沿用本地约定 `k_table[MeV/c] = 500*sqrt(q_out²+q_side²+q_long²)`。这只是明确沿用的查询约定；输入为 LCMS q 时，不把此量宣称为经过严格 PRF 变换的 k*。

`with_coulomb=true` 使用本地相同电荷 ππ 库仑波函数；false 使用平面波。两者均启用 identical-pion QS，质量与本地参考一致（139.57 MeV）。当前没有其他粒子或强相互作用配置。

## 配置和运行

参考 `config/examples/cats3d_gaussian.toml` 和 OO 专用配置 `config/oo_build_and_fit_6bins_cats3d.toml`。后者保留原有 sel8 CF 输入和切片选择，给新输出加 `_cats3d` 后缀。**示例 λ=0.7 不是生产标定。**

```sh
bash scripts/run_exp_femto_3d_cats3d.sh
# 或使用已有入口，显式复用已构建 CF：
bash scripts/run_exp_femto_3d.sh --stage fit --config config/examples/cats3d_gaussian.toml
```

两个 CATS3D 入口脚本的 `config_path`、`log_directory` 均可直接编辑；运行时以实际变量和命令行覆盖为准。`run_exp_femto_3d_cats3d.sh` 与 `run_oo_cats_6bins_shared_lambda_cuttest_scan.sh` 共用 `run_cats3d_with_progress.sh`：继续通过 guarded executor 进入 O2Physics，在其内部使用 macOS 伪终端启用已有的自动进度条，并将输出实时转发至终端。

- 启动时打印实际命令和两份日志路径。`logs/cats3d/cats3d_<时间>_<PID>.log` 保存运行输出和进度；同名前缀的 `.runtime.log` 保存环境检查、实际调用和 executor 最终状态。同名覆盖通过时间与 PID 避免。
- 现有进度条显示切片完成比例、活动状态和 ETA；每个 λ 内部的 Minuit 迭代或源卷积进度目前没有独立计数。现有 `log_level` 和模型/统计量设置保持不变。
- `--dry-run` 必须放在入口脚本参数的第一位，只显示将使用的命令及日志路径，不创建日志、不启动 ROOT；其他拟合参数继续转发。
- 启动脚本的修改只影响之后启动的任务。已经运行的任务保留其原输出方式和环境，不会因脚本更新而自动重启或获得新的进度条。

最小新增配置：

```toml
[fit]
correlation_backend = "cats3d"
model = "full"
use_pml = true
use_q2_baseline = true
fit_q_max = 0.15

[fit.cats3d]
with_coulomb = true
grid_bins = 31
grid_min_fm = -50.0
grid_max_fm = 50.0
subdivisions = 1
max_pw = 16
pml_max_calls = 100000 # PML MIGRAD call budget; unsuccessful points remain recorded.
max_kernel_gib = 4.0
max_peak_gib = 8.0
use_disk_cache = true
# kernel_cache_file = "/path/to/matching_native_kernel.root"

[fit.cats3d.lambda]
kind = "constant"
meaning = "final_lambda"
value = 0.7 # 示例
```

CATS3D 仅接受独立切片流程。`shared_phi`、`combined_profile`、旧 `profile_likelihood` 和 `fit.parameters.lambda` 会报错。也不能指定 Bowler 专用的 `fit.coulomb_mode`、`fit.finite_source_mode`、`fit.use_coulomb` 或 `fit.use_core_halo_lambda`；库仑开关由 `fit.cats3d.with_coulomb` 唯一控制。

半径、normalization、baseline 的初值及界限沿用 `fit.parameters`。λ 与 α 不计入自由参数。PML 继续使用项目现有原始 SE/ME 计数及 raw-to-normalized 比例；χ² 继续使用 CF 的已存误差。CATS3D 的 legacy-Minuit PML 调用设置 `SET EPS=1e-8`，使数值导数避开网格累加及 near-zero count-log 目标函数的舍入噪声；不改变 PML 公式或协方差质量判据。CATS3D 的 χ² 使用 ROOT 默认拟合，不启用额外的 `M`/IMPROVE 后搜索；该后搜索有独立的状态码（见 [ROOT Fit 状态定义](https://root.cern/doc/master/classTH1.html)）。无需引入外部 CATS 拟合器的第二套统计量。

## λ 输入含义和组绑定

`kind="constant"` 需要 `value`；`kind="histogram"` 需要 `file`、`histogram`、`axis_unit="GeV"|"MeV"`。输入含义必须显式指定：

- `final_lambda`：直接使用，禁止再给 `weak_decay_fraction`。
- `pair_purity`：输入是对纯度，必须显式给 `weak_decay_fraction`；最终 `lambda_main = input*(1-f_weak)^2`。这里不接收单粒子纯度，也不内置 0.12。

表在 bin 中心间线性插值；从轴端点到首末中心使用最近中心值。超出表轴范围则在运行前报覆盖不足，不能默默夹到端值。普通 bins 的值必须有限并在 [0,1] 内，误差必须有限非负；不使用 underflow/overflow。

```toml
[fit.cats3d.lambda]
kind = "histogram"
meaning = "pair_purity"
file = "/path/to/purity.root"
histogram = "hPairPurity"
axis_unit = "GeV"
weak_decay_fraction = 0.12 # 仅为显式配置示例

[[fit.cats3d.lambda_groups]]
group_id = "从 meta/SliceCatalog 读取的精确 group_id"
kind = "histogram"
meaning = "final_lambda"
file = "/path/to/lambda.root"
histogram = "hLambdaForThisGroup"
axis_unit = "GeV"
```

组覆盖优先于全局默认值；没有默认值时，所选每组都必须有精确绑定。重复、空或 catalog 中不存在的 group_id 会报错。同组的各 φ 和 φ 积分切片共用输入。相对输入文件及缓存目录相对于 TOML 所在目录解析。

原始表及其误差会保存在输出中。HESSE 协方差和参数误差**条件于外部固定 λ 和数值核**，不自动传播表的误差。常数 λ 的拟合误差为 0 表示“固定参数”；并不意味着外部标定误差为零。表输入的旧标量 `lambda`、`lambda_err` 栏为 NaN，不填平均值。

## 独立诊断扫描

```toml
[fit.cats3d.lambda_scan]
enabled = true
execution_mode = "alongside_fit" # 或 scan_only
root_name = "cats3d_lambda_scan.root"
constant_values = [0.6, 0.7, 0.8]
scale_values = [0.9, 1.0, 1.1]
```

- 常数坐标就是最终 λ，在扫描中替代原始常数或表，包括替代原先的纯度换算结果。
- 倍率坐标作用于已完成纯度换算的 λ 表：`lambda_scan(k)=s*lambda_main_input(k)`。所选每组都必须使用 histogram 输入。
- 两个列表是两条独立扫描，不是笛卡尔积。坐标必须有限且不重复，内部升序执行。
- 每点固定 λ，重新最小化其余自由参数；使用配置种子、正向相邻点、反向相邻点重试。只在成功的有效最小值中选最小目标函数。
- 越界 λ、失败最小化及坏协方差均保留。有限但超出 [0,1] 的试探 λ 被记为无效点，不夹值或替换为惩罚高度；未执行 Minuit 的点也有失败原因。λ=0 时半径可能不可识别，不能把其协方差当有效误差。
- 每组总目标函数是所有非积分 φ 成员的条件最小目标函数之和；φ 积分切片另存。成员任一失败则整个组点无效。组扫描要求使用独立的 φ-binned ME，防止把重复使用的 ME 当作独立联合似然。
- Δ目标函数仅相对于同一条曲线的有效扫描点最小值；没有自动区间估计、置信度标签或正式 λ 选优。曲线在无效点断开。

`scan_only` 私下计算所需参考点并写扫描文件，不创建或改写正式 fit ROOT、TSV 和 report ROOT。扫描文件必须与所有输入及正式输出路径不同；λ 扫描默认关闭。

## ROOT 输出

正式文件继续使用 `meta/FitCatalog`、`fits/<slice_id>` 和原有 summary 组织，TSV 原有列顺序不变，在末尾追加后端与 λ 来源字段。

每个有效切片新增：

- `CATS3DCoreCF3D`、`ModelCF3D`、`LambdaMap3D`：同一计算 q 网格上的核心、最终模型和实际 λ。
- `LambdaInput`：ROOT 表原始副本（保留原轴、bin 内容及误差）。
- `LambdaInputMetadata`、`CATS3DContract`、`PhysicalValidationBoundary`：输入含义、校正、核网格/物理设置及物理边界。`CATS3DKernelFile` 保存实际缓存来源。
- `Fit_ProjX/Y/Z` 与 `Canvas_ProjX/Y/Z`：从最终模型的实际 bins 做投影；不使用回退到 Bowler 的 TF3。CATS3D 输出不创建 `LevyFit3D`。

失败的正式拟合仍写 catalog 行与 `FitFailure`，不制造有效模型。组内 φ 拟合或误差不完整时，summary 写明状态，不对不完整成员做 harmonic/epsilon 拟合。报表保存 `cats3d_lambda_inputs/<group_id>` 及详细模型文件位置。

扫描文件：

- `inputs/<slice_id>`：名义外部 λ 来源和 map。
- `lambda_scan/{constant|scale}/<group_id>/<slice_id>/{ScanPoints,AttemptPoints}`：逐点结果及所有实际尝试。参数顺序见同一切片输入目录的 `ParameterContract`；λ/α 为固定位置，表模式下 `ExternalLambdaCarrierFixed` 不是表的平均值。
- `lambda_scan/{constant|scale}/<group_id>/GroupPoints`：完整成员列表、组目标函数及有效标志。
- `DeltaObjectiveSegment_*`、`DeltaObjectiveCanvas`：带缺口的诊断曲线。逐点表还标记扫描端点和位于端点的扫描最小值。

正式文件与扫描文件用 `CATS3DRunStatus` 标记 incomplete/complete；扫描输入目录保存参数初值、界限、固定标志、几何与统计量定义。正式切片的 `FitPerformance` 保存名义拟合耗时、源卷积次数、核字节数和峰值内存估算。

### 从已完成扫描生成逐 λ 报表

扫描主流程只自动写 Δ目标函数曲线和逐点参数，不自动调用正式拟合报表。可以在扫描结束后运行独立后处理：

```bash
bash scripts/run_cats3d_scan_report.sh \
  --scan-root /path/to/lambda_scan_cats3d.root \
  --output-directory /path/to/new_cats3d_report_directory
```

脚本通过 O2Physics guarded executor 构建并运行 `bin/cats3d_scan_report`，再生成可离线打开的 `index.html`。默认输入为 OO cuttest 扫描；省略输出目录时使用 `res/cats3d_scan_report_<时间戳>`。`--dry-run` 仅预览命令。输出目录必须是新目录，以保留已有报表。

每个扫描坐标生成独立 PDF 和 ROOT 文件，包含组内非积分 φ 的参数分布，以及每个切片的 out/side/long 数据与模型投影、参数、误差、边界和 Minuit 状态。ROOT 内的 `slices/<slice_id>` 保留该点的原始 `ScanPoints`、全部 `AttemptPoints`、投影及 `Canvas_FitDiagnostics`；`groups/<group_id>/ParameterVsPhi` 保存参数总览。`scan_points.csv` 与 `manifest.json` 提供逐点数据和图像索引。

PDF 使用同一批 ROOT PNG 的无损快照，以保持页面尺寸和图像布局一致；需要矢量画布时使用配套 ROOT 文件。网页/PDF 组装只使用 Python 标准库，在 ROOT 环境外执行，并清除子进程的 `PYTHONHOME`/`PYTHONPATH`，避免 O2 Python 模块与普通 Python 可执行文件冲突。

后处理不重新最小化切片模型，而是从扫描文件读取参数、几何、baseline、原始 CF 路径和核设置，使用相同的 `Cats3DModel` 适配器重新卷积。核采用记录路径的 load-only 模式；缺失或不匹配时明确停止，不自行生成。每次只持有一个核，沿用 4 GiB 单核 / 8 GiB 峰值估算限制。表输入从扫描文件中保存的 `LambdaInput` 复制，避免原始标定表后来变化。

有效解绘制橙色实线；未收敛但模型域仍有效的保存试探点仅绘制灰色虚线并标明诊断用途。缺失误差保持不可用，不补零，也不绘制误差带。数据和模型平均同样的横向普通 bins（±0.04 GeV/c）；模型曲线止于已计算区域。φ 积分切片单独显示，不混入组内参数图。该报表不自动选择正式 λ，也不改变扫描有效标志或正式拟合输出。

今后可在扫描命令成功后串联此脚本；`alongside_fit` 的正式报表仍只对应名义 λ，不能代替逐扫描坐标的报表。目前没有“开启 TOML 开关便自动生成逐 λ 投影报表”的主流程配置项。

该工具还会为每个坐标生成传统的 `fit_report_cats3d_*.root`，与 `res/fit_report_*.root` 使用相同的 `meta/FitCatalog`、`source_parameters/cent_*/mt_*/source_parameters_overview_canvas`、`summary/R2_vs_phi` 和 `eps_vs_mt` 组织。源参数总览沿用 α、三个对角项和三个交叉项的 2×4 布局，明确标注 CATS3D、高斯源和固定 λ。φ 积分切片保留在目录和参数表中，不参与 φ 谐波拟合。

只有所有所选非积分 φ 成员的极小值和条件误差均有效、且至少有三个成员的组才进行原报表的谐波与偏心率拟合。不完整组绘制有有效条件误差的点，另存 `conditional_error_status_canvas`，不对剩余点做谐波/偏心率拟合。全部失败解和缺失误差仍保存于 `FitCatalog` 及逐切片扫描记录。表输入的 λ 标量栏保持 NaN，并保留原输入表与坐标的倍率含义。

只生成这种传统 ROOT 汇总可以运行：

```bash
bash scripts/run_cats3d_scan_report.sh \
  --scan-root /path/to/lambda_scan_cats3d.root \
  --output-directory /path/to/new_fit_report_directory \
  --summary-only
```

此模式不加载 CF、kernel，不计算源卷积或重新最小化切片模型，也不生成 PDF/HTML；仍会执行符合条件的源参数谐波拟合。常规报表命令仍会重建模型投影，并附带上述传统 ROOT 汇总。

## 核缓存与验证边界

缓存采用本地 CATS3D 原生 `CATS3DFastKernel6DByQ` 格式，不再要求 `ExpFemto3DContract`、内容摘要或 `libCATS` 二进制指纹。自动文件名直接表示 q 网格、空间网格及波函数设置；`CATS3DContract` 只是拟合结果的来源说明，不是加载门槛。

默认缓存目录仍是当前 CMake build 下的 `cats3d_kernels`，可用 `kernel_cache_directory` 改位置。显式指定原生文件：

```toml
# 在已有 [fit.cats3d] 段内设置
use_disk_cache = true
kernel_cache_file = "/path/to/matching_native_kernel.root"
```

`kernel_cache_file` 优先于自动缓存目录，按只读输入处理；文件不存在或不匹配就报错，不生成或覆盖它。原 CATS_3D 生成的 ByQ 文件、以前带私有标记的缓存都可指定任意文件名读取。读取时保留原生元数据核对（q 轴、r 中心/间距、质量、电荷、QS、波函数、max_pw、subdivisions）、q 记录完整性、向量长度及有限非负值检查；这些是在检查科学模型是否一致，不绑定生成文件的电脑或库二进制。跨库版本的物理等价性需要另行核验。新自动缓存仍以临时文件原子提交。

已有旧核 `dataKernel_qData_r31_coulomb_byQ.root` 的 q 网格是 20³、20 MeV/c 间隔，当前 cuttest 数据是 40³、10 MeV/c 间隔。这一差别在去除私有指纹后仍会明确拒绝，不自动改变数据分箱或插值波函数。

每次只持有一份当前 q 网格的 CATS3D 实例，跨切片/组/扫描复用波函数核；源矩阵变化重新卷积。只有六个源矩阵元素完全相同时才复用关联函数，不采用半径容差缓存。normalization、baseline 和 λ 的变化不重建源卷积。

默认源网格 31³、±50 fm 和 max_pw=16 是兼容参考值，不是通用精度保证。40³ 的 cuttest CF 在 ±0.15 GeV/c 拟合盒内有 30³ q bins；其浮点核为 **3,217,428,000 字节（2.996 GiB）**。按已检查旧核的压缩比粗估，磁盘约 2.37 GiB；新核尚未生成，该值不是实测。

原生 ByQ 加载器先读入一份 float 向量，再复制进 CATS3D，因此加载瞬时可能同时持有两份完整核。当前峰值估算为 `2*kernel_bytes + 16*空间采样点数 + 1 GiB`，最后 1 GiB 为 ROOT、输入对象及其他工作区预留。当前网格估算约 **6.99 GiB**；卷积时只保留一份完整核。以前的 `3*kernel_bytes` 估算已被上述经加载代码核对的两份副本估算替代。

`max_kernel_gib=4` 限制单份核，新增 `max_peak_gib=8` 限制上述估算峰值；任意一项超限都在分配波函数核前拒绝。它不是操作系统的 RSS 硬上限，估算也不包含其他应用的内存占用。输入 CF、SE、ME 在预检和逐切片处理时由局部对象持有并及时释放，避免随所选切片数积累副本。本机 16 GiB 使用单进程串行扫描，同一核跨切片和 λ 点复用；不要为每个 λ 启动一个并行进程。全尺寸峰值尚未实测，不以估算承诺可用物理内存；需要更严格余量时可降低 `max_peak_gib` 使启动预检明确拒绝。

本次保持给定 CATS3D 波函数实现，包括当前未启用的库仑相移项 TODO。适配器与该实现闭合，不代表此波函数、LCMS/PRF 近似或任意生产网格已独立完成物理验证。正式分析需检查空间范围、间距和 partial-wave 收敛，并提供真实 λ 标定。

数值和工作流测试通过同一 O2Physics executor 运行：

```sh
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d \
  --command 'ctest --test-dir build --output-on-failure'
```

本次自动测试、真实组短扫描和空间网格检查见 [CATS3D 验证记录](CATS3D验证记录.md)。

扫描结果 ROOT 保存参数、诊断表和 λ map，不嵌入整份波函数核；增加 λ 点数不会按每点再存一份约 3 GiB 核。
