# OO：PbPb λ(k)、共同 mT 分箱与 box 归一化

新配置：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/config/oo_cats3d_6phi_pbpb_lambda_box_commonmt.toml`。

输入为 loosecut05thres OO，中心度 0–10、10–30、30–50%，φ 原生 18 bins 合并为 6 bins。mT 使用 0.2–0.4、0.4–0.6、0.6–1.1。CATS3D full 高斯源、α=2、PML、q² baseline，fit_q_max=0.15 GeV/c；λ 扫描关闭，波函数核缓存沿用原路径。

## 纯度表

准备宏：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/scripts/prepare_oo_pbpb_lambda_tables.C`。输入目录为 `/Users/allenzhou/Research_software/CATS/CATS_3D/purity_PbPb`。

| 合并 mT | PbPb 原始 mT |
|---|---|
| 0.2–0.4 | 0.2–0.3 + 0.3–0.4 |
| 0.4–0.6 | 0.4–0.5 + 0.5–0.6 |
| 0.6–1.1 | 0.6–0.7 + 0.7–0.9 + 0.9–1.1 |

每个 k bin 使用对应 `hWeightKstar` 加权：`P=Σ(WᵢPᵢ)/ΣWᵢ`，不额外乘 mT 区间宽度。采用 PbPb 文件内的权重；这是借用 PbPb 输入的 OO 诊断尝试，未进行 OO 权重重估。

纯度误差采用固定权重、输入误差独立的近似：`σ²=Σ(Wᵢ²σᵢ²)/(ΣWᵢ)²`。原文件未提供输入协方差；权重误差单独保存，不计入该纯度误差。外部纯度表误差不自动传播到拟合的条件误差。

合并表保留原 k 轴边界，仅移除总权重为零的尾部。内部或低 k 零权重、缺失对象、轴不一致及覆盖不足均报错。正常 bin 参与计算，underflow/overflow 不参与。保留所有原始对象到 `sources/cent*/`，来源路径和 MD5 随文件保存；九张纯度表及总权重在 `tables/`，合并规则在 `meta/`。

存储值仍是 pair purity。TOML 对九组逐一指定 `meaning="pair_purity"`、`axis_unit="GeV"`、`weak_decay_fraction=0.12`，使最终 λ(k)=P(k)×0.7744。准备宏不预先乘这一因子。

## 新输出路径

运行标识：`oo_cats3d_pbpb_lambda_box_commonmt_20261023`。

- CF：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/cf/oo_cats3d_pbpb_lambda_box_commonmt_20261023/cf_oo_cats3d_pbpb_lambda_box_commonmt_20261023.root`
- 拟合 ROOT：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/data/oo_cats3d_pbpb_lambda_box_commonmt_20261023/fit_oo_cats3d_pbpb_lambda_box_commonmt_20261023.root`
- 汇总 TSV：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/data/oo_cats3d_pbpb_lambda_box_commonmt_20261023/fit_summary_oo_cats3d_pbpb_lambda_box_commonmt_20261023.tsv`
- 报表 ROOT：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/res/oo_cats3d_pbpb_lambda_box_commonmt_20261023/fit_report_oo_cats3d_pbpb_lambda_box_commonmt_20261023.root`
- 合并纯度表：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/data/oo_cats3d_pbpb_lambda_box_commonmt_20261023/pbpb_pair_purity_oo_cats3d_pbpb_lambda_box_commonmt_20261023.root`

box_shell 的 out/side/long 范围均为 [0.15,0.20] GeV/c，区域是外盒减内盒。正式启动使用 `--stage all`，先从原始 OO 数据重建 CF。fit-only 会读取输入 CF 自带的归一化元数据；因此不能使用旧 CF 来代替这一步。

## 准备与启动命令

先设置以下变量（可在任意工作目录执行）：

```bash
project="/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d"
tag="oo_cats3d_pbpb_lambda_box_commonmt_20261023"
config="$project/config/oo_cats3d_6phi_pbpb_lambda_box_commonmt.toml"
source_dir="/Users/allenzhou/Research_software/CATS/CATS_3D/purity_PbPb"
purity_root="$project/data/$tag/pbpb_pair_purity_$tag.root"
executor="/Users/allenzhou/.codex/skills/cern_root/o2physics-root/scripts/run_root_command.sh"

mkdir -p "$project/cf/$tag" "$project/data/$tag" "$project/res/$tag"
```

合并纯度表已准备好。以下是可复用的生成命令；宏拒绝覆盖现有文件，重新生成时使用新的输出路径，并同步修改 TOML 中九组的 `file`：

```bash
bash "$executor" --cwd "$project" --command \
  "root -l -b -q '$project/scripts/prepare_oo_pbpb_lambda_tables.C(\"$source_dir\",\"$purity_root\")'"
```

独立启动脚本：`/Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/scripts/run_oo_cats3d_pbpb_lambda_box_commonmt.sh`。可在任意工作目录直接运行，默认 `--stage all`，使用本配置及专属日志目录，正式启动时自动建立 CF/data/res 目录。

启动预览（不创建目录或日志、不构建 CF、不拟合）：

```bash
bash /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/scripts/run_oo_cats3d_pbpb_lambda_box_commonmt.sh --dry-run
```

正式运行（经 guarded executor，显示实时进度，先构建新 CF 再拟合）：

```bash
bash /Users/allenzhou/Research_software/Code_Base/Exp_femto_3d/scripts/run_oo_cats3d_pbpb_lambda_box_commonmt.sh
```

需要复用本次已构建的新 CF 时，显式加 `--stage fit`。`--model diag`、`--config FILE` 等参数会转发至原入口；自定义配置若使用其他绝对 CF 目录，应先建立该目录。启动脚本不重新生成或覆盖已准备好的纯度表。

可重复的表验证：

```bash
bash "$executor" --cwd "$project" --command \
  "root -l -b -q '$project/tests/PrepareOoPbPbLambdaTablesTest.C(\"$source_dir\",\"$purity_root\")'"
```

本次交付执行了纯度表准备、独立回读验证和启动预览；CF 重建及正式拟合留待上述正式命令执行。
