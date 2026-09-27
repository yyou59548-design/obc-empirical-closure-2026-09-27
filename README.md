# OBC 实证闭环：结果、过程和可核查证据

这是 2026-09-26 OBC 实证闭环任务的归档入口。**结论是原任务允许的科学阻塞停止条件，而不是论文完成**：`paper_ready=false`，42 个无条件 mandatory gates 中 0 个完整通过，合格非线性端点后验 0/3。没有新的 headline 经济结果。

## 给 GPT 的最短阅读路径

1. [最终科学结论与定量证明链](results/OBC_SCIENTIFIC_BLOCKER_FINAL_CN.md)
2. [机器可读最终状态](results/OBC_FINAL_STATUS.json)与[逐项 mandatory gate 账本](results/OBC_MANDATORY_GATE_LEDGER.json)
3. [原始任务说明](task/START_HERE_CN.md)、[主任务](task/CODEX_MASTER_TASK.md)、[验收门槛](task/ACCEPTANCE_GATES.json)
4. [可在线浏览的关键证据、勘误、源代码与核验脚本](evidence/README.md)；按索引顺序阅读，并注意有限支撑结论的范围
5. [小型原始证据包](results/OBC_SCIENTIFIC_BLOCK_EVIDENCE_v1.zip)及其[校验清单](results/OBC_SCIENTIFIC_BLOCK_EVIDENCE_v1.manifest.json)，用于取得本仓库未展开的二进制证明材料

若要独立复核完整过程，再下载本仓库 **Releases** 中同一版本的全部 `OBC_ALL_RESULTS_AND_PROCESS_20260927.zip.part*` 文件，按 [分片清单](release/OBC_SPLIT_MANIFEST.json) 校验并用 [`tools/reassemble.py`](tools/reassemble.py) 还原。以清单中的 `part_count` 和 `restore_order` 为准。完整 ZIP SHA-256：`649301cf72692f17515989406cbc017f117af8694aa9eabe34b3895e4844912d`；它收录了原始上传包、解压任务、运行目录、日志、原子 checkpoint、脚本、证据和最终状态。完整归档的[外部清单](release/OBC_ALL_RESULTS_AND_PROCESS_20260927.manifest.json)记录 107,697 个来源文件和逐文件读回校验。

## 结论范围

锁定的原 `find_lk` **最终退出约束**的均衡选择规则，在 39 维正先验参数盒与 1964Q1 首次转移的正概率状态/冲击事件上，找不到满足互补条件且最终回到全松弛的有效非负政策楔子路径。区间证明给出的关键上界为：`q·y ≤ −0.02606566006075801`；前 1,024 个未来列的最大上界 `−0.00012869994399417693`；其余所有未来列的归一化上界 `−2.0531122580970638`。据此，锁定的**选定解映射**在正先验预测质量上未定义。

这个结论不声称所有可能的永久约束、无界或违反该选择规则的数学均衡均不存在，也不声称已知此事件在给定观测数据后的后验质量。若改变模型、选择规则或随机状态/冲击域，应重新定义并认证新的科学目标。旧 trace、旧 shocks、旧 P/F 归因、线性 posterior 和未收敛链均不能当作本任务的新端点结果。

## 如何向 GPT 提问

可直接提供这个仓库链接并指定：

> 请先阅读 README、`results/OBC_SCIENTIFIC_BLOCKER_FINAL_CN.md`、`results/OBC_FINAL_STATUS.json`、`results/OBC_MANDATORY_GATE_LEDGER.json`，再按 `evidence/README.md` 查验关键证据与勘误。将 `task/` 下的文件视为原任务规范和审计材料，并以我这次的问题为当前指令。请区分已证明的锁定选择映射阻塞、尚未通过的验收门槛，以及未声称的更广泛结论。若需要重算原始数值证据，再读取小型证据 ZIP；只有确需完整过程时才下载 Release 分片。

仓库内的摘要和证据 ZIP 可独立浏览或下载；各个大分片只是同一个完整 ZIP 的顺序字节切片，单个分片无法作为 ZIP 打开。
