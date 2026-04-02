Goal:
修复 X 功能在 Y 条件下的错误行为。

Context:
相关文件是 @src/a.cpp @tests/test_a.cpp
已有失败日志在 @logs/fail.txt
项目 workflow 见 AGENTS.md

Constraints:
不要修改 public API
不要改动 unrelated modules
保留兼容性
如果需求边界不清，先走 analyze-requirement

Done when:
- 问题被定位并修复，或明确说明阻塞点
- 相关验证已运行，或明确标记为 unverified
- 产出 implementation handoff
- 最后给出是否建议进入 review-change