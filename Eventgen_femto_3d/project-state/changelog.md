# Changelog

## 2026-04-13

- 完成 `Eventgen_femto_3d` 的 `CLI + TOML` 配置化改造。
- 完成 blast-wave `events + particles` 输入适配，并保留 legacy 输入兼容。
- 将输入层从 `Workflow` 中解耦，统一归一化为 `EventData`。
- 新增配置解析、输入适配、legacy smoke 与 ROOT runtime guard 测试。
- 在 O2/ROOT 环境下完成构建、全量测试与 wrapper 脚本 smoke。
- 补写 `project-state`，记录本轮实现结果与流程改进项。
- 补齐 `project-state/guide.md`，完成面向维护者的人类可读概览。
- 调整了可执行二进制文件的位置到 `bin` 目录下。
