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
- 同步 `project-state` 到当前 `HEAD`，将旧的 `6/6` 测试描述更新为当前注册的 `8` 项测试。
- 记录 `run_eventgen_femto_3d.sh` 已切换为调用 `bin/eventgen_femto_3d`。
- 补记一次当前 `HEAD` 下的 `ctest` 与 `alienv` wrapper smoke 证据。
