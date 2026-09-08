# 二维 profile/slice 诊断显示合同

状态：2026-09-07 已实施。本文记录问题 2 的当前实现边界。

## 数值语义

Profile 在固定两个目标参数后重新最小化其余 free nuisance；
fixed-nuisance slice 保持 nominal nuisance。两者共用同一 slice 内的
nominal/profile 全局诊断 reference，绘图阶段不再各自减最小值。
`ProfilePoints` 和 `AttemptPoints` 保留精确坐标、objective、状态和失败原因，
是数值真源。

## 显示语义

- TH2 外边界等于 resolved scan 上下限，内部边界是相邻采样坐标中点。
- coarse 有效 delta 写原值，不可用 bin 保持 NaN；underflow/overflow 不使用。
- profile 阈值画布为 likelihood/status 双面板，另存全有限值范围画布。
- slice 使用自身有效性 mask，不复用 profile 收敛状态。
- invalid 显示为不透明灰色，不使用 0 或 penalty 伪装 likelihood。
- contour 仅在四角均有效且有限的 coarse 单元中以真实坐标插值，
  saddle 用双线性中心值消歧，不跨越缺失区。
- best 从 coarse/refined 全部有效点按 delta、stage、point index 稳定选择。
- 1/2/4 等阈值只是诊断线，不作置信区域解释。

## ROOT 对象与兼容性

Profile 固定写出 `DeltaNeg2LogL2D`、`PointStatus2D`、`Canvas_2D`、
`Canvas_2D_FullRange`；有有效点时写 `BestProfileGridPoint`。
`write_likelihood_slice=true` 时额外写 `SliceDeltaNeg2LogL2D`、
`Canvas_Slice2D`、`Canvas_Slice2D_FullRange`，并在有效时写
`BestSliceGridPoint`。`NominalPoint` 保留真实坐标，越出 scan 时只由画布裁剪。

Checkpoint 使用 `profile-contract-v4|display-contract-v3`；旧 v2 chunk 必须拒绝。
本轮不改变 PML、Levy/PSD/Coulomb 模型、Minuit seed/retry、网格密度或
refinement 策略，不执行真实 OO profile 作业。
