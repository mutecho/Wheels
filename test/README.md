# OO same / mixed 分布比较

A 为参考文件，B 为本地结果：

- A: /Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO/AnalysisResults_3d_25ae_pass2_loosecuttest.root
- B: /Users/allenzhou/ALICE/alidata/hyperloop_setup/femtoep_oo/run/res_20261020/AnalysisResults.root

默认比较两者的 **femto-dream-pair-task-track-track**，读取
**SameEvent_3Dqn/relPair3dRmTMultPercentileQnPairphi** 和
**MixedEvent_3Dqn/relPair3dRmTMultPercentileQnPairphi**。
两个输入的普通 SameEvent/MixedEvent 分布为空，不能代表实际 3D 配对结果。
A 另外含有 _04thres task；默认不与主 task 混合。

## 运行

~~~bash
bash /Users/allenzhou/Research_software/Code_Base/test/run_compare.sh
~~~

入口通过共享 O2Physics ROOT executor 加载环境，使用 CMake 编译、运行数值自检，
再生成结果。可执行文件位于 bin/compare_same_mixed。
若执行器提示 ESCALATION_REQUIRED，应在允许提权的环境中重跑同一命令。
已经进入正确 O2Physics 环境时也可运行 bash build_and_run.sh。

默认结果目录为本目录下的 comparison_25ae_vs_local_20261020/。
重复运行会更新该目录中的同名生成文件。

比较 A 的另一个 task：

~~~bash
bash /Users/allenzhou/Research_software/Code_Base/test/run_compare.sh \
  --task-a femto-dream-pair-task-track-track_04thres \
  --label-a "A: 25ae 04thres" \
  --output comparison_04thres_vs_local
~~~

也支持 --file-a、--file-b、--task-b、--label-b；
--help 查看参数，--self-test 仅运行数值自检。

## 图和数据

- **comparison.pdf**：完整多页图册。
- **png/00_overview.png**：六个主要变量的总览。
- **png/00_pair_yields.png**：原始 mixed/same 总数比及其中心度、mT 分布，保留总体产额差异。
- **png/01_*.png**：逐变量详图；上方分别为 same / mixed 的 A、B 叠加分布，下方为 B/A 和双比值。
- **png/10_*.png**：小动量差立方体、各中心度、各 mT、中心度与 mT 的组合分组图。
- **png/20_*.png**：动量分量二维投影、mT–中心度、方位角–中心度、方位角–mT 的二维分布和 B/A。
- **png/30_*.png、png/31_*.png**：close-pair rejection 前后、PV 前后的 Δη / Δφ QA 投影和二维图。轴标题沿用输入文件。
- **comparison.root**：全部原始投影、归一化密度、比值 TGraphErrors 及完整 bin 有效性 mask；也保存所有二维原始投影及 B/A。
- **metrics.csv**：每个分组、每个变量的原始 pair 总数和形状距离 TV。
- **bin_values.csv**：原始 bin 计数、same/mixed B/A、双比值、近似误差和有效性标记。
- **pair_yields.csv**：未做单位面积归一化的 mixed/same 产额比。
- **inventory.csv、flows.csv、axes.csv**：输入计数、各轴欠流/溢流权重、完整 bin edges。
- **inputs.txt**：输入路径、选用 task、各文件实际包含的 pair task。
- **plot_index.csv**：PDF 页码与 PNG 名称对应表。
- **summary.md**：结果数值和比较口径。

## 比较口径

1. 三个动量差轴的原始范围为 [-0.2, 0.2) GeV/c。主比较要求七个轴全部位于普通 bin；
   不将其他轴的欠流或溢流混入某个投影。原始计数与排除比例保存在 inventory / flows。
   不同轴的 flow 类别可能重叠，不可直接相加。QA 使用各自二维普通 bin 的范围，未施加 7D 选择。
2. 保留所有原始 bin edges，不做平滑、插值或重新分 bin。mT 与中心度分组采用完整原始 bin。
   小动量差图要求三个分量都在 [-0.05, 0.05) GeV/c；只有边界与原始 bin 对齐时才生成。
3. 每个变量、每个分组、每种 sample 分别归一化至单位面积；
   纵轴为 n_i/(N Δx_i)，可直接比较变宽 mT / 中心度 bin 的密度。
   原始 pair 数保留在详图、ROOT 和 CSV 中。
4. R_SE = (SE_B/N_SE_B)/(SE_A/N_SE_A)，mixed 同理。
   双比值 D = R_ME/R_SE；偏离 1 表示 mixed 的形状变化与 same 不同。
   这是形状诊断，没有进行物理相关函数所用的高动量区间归一化。
5. 分母为零的 bin 不画比值点，CSV 留空并标记无效，ROOT 保存 mask。
   二维图用灰色斜线标记无定义 bin。二维比值色阶固定为 0.8–1.2 并标注超出范围的 bin 数；
   ROOT 中保存完整、未经裁剪的比值。
6. 保留输入已有 Sumw2；没有 Sumw2 时采用 ROOT 对无权计数的方差约定。
   单位面积误差包含 Cov(n_i,N)=Var(n_i)。两个输入可能重叠，
   same/mixed 共享 tracks，且 pairs 可能被重复使用；这些协方差无法从现有直方图恢复。
   图中误差仅描述计数精度，不能直接解释为差异的显著性。
7. TV=0.5*Σ|p_A-p_B|，范围 [0,1]，描述归一化形状差异；
   它包含统计涨落且依赖 binning。qn 仅有一个 bin，内部形状没有信息。

自检覆盖单位面积归一化的协方差、变宽 bin 密度、零分母 mask、
完整 bin edges 检查、七维 flow 排除、小动量差选择和权重 Sumw2 保留。
真实运行还检查每个七维投影的计数守恒。
validate_comparison.C 另外使用 ROOT 原生 Projection，对四个输入 sample 的
两个一维投影与一个二维投影逐 bin 交叉核对。通过共享执行器运行：

~~~bash
bash /Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh \
  --cwd /Users/allenzhou/Research_software/Code_Base/test \
  --command 'root -l -b -q validate_comparison.C'
~~~
