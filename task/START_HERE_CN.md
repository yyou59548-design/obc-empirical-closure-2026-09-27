# Codex：OBC 端点学习实证闭环总任务（2026-09-26）

## 0. 任务目标

你的任务不是“继续试几步 MCMC”，也不是写论文。你的任务是：**在不伪造、不借用旧结果、不降低验收标准的前提下，把本项目成文前的全部实证任务推进到底，并交付一个可审计的 paper-ready 结果包。**

当前真实状态为：`paper_ready=false`，三个主端点的合格 nonlinear posterior 为 **0/3**。最新状态见 `REFERENCE/CURRENT_STATUS_CN.md`。不要相信文件名里的 `final`、`certified`、`production` 等历史措辞；只以本包的任务规范和重新执行的证据为准。

## 1. 开始前必须做

1. 阅读：
   - `CODEX_MASTER_TASK.md`
   - `ACCEPTANCE_GATES.json`
   - `DO_NOT_REUSE.md`
   - `OUTPUT_REQUIREMENTS.md`
   - `REFERENCE/CURRENT_STATUS_CN.md`
2. 运行：
   ```bash
   bash SCRIPTS/bootstrap_and_verify.sh
   ```
3. 将两个 INPUT zip 解压到独立目录。不要覆盖原 ZIP；原输入只读。
4. 首先复核输入哈希、环境、依赖和已有 checkpoint 可读性。
5. 所有新运行放在新的工作目录，例如 `WORK_20260926/`；旧 trace 不得拼入新正式 trace，除非本规范明确允许“仅作为位置初始化”。

## 2. 最重要的科学边界

### 2.1 绝对禁止

- 不得把旧 v5 trace 当正式 posterior。
- 不得把历史 `F.T + pseudodensity` 链当新的 scientific posterior。
- 不得使用“分别生成三套 Sobol 再按同一行配对”的旧 QMC 构造。
- 不得把旧 6408/6420 的 shocks/residuals 套到新参数上。
- 不得把线性 Kalman posterior 改名为 nonlinear OBC posterior。
- 不得用超时自动拒绝 proposal。
- 不得用未收敛 source chain 做 PSIS 后宣布 target posterior 完成。
- 不得把失败 coalition/path 填成 0。
- 不得把旧 P/F attribution 改名填进新三端点结果表。
- 不得因为 risk-premium 排名“好看”而改变 sampler、N、seed、branch range、draw selection 或 stopping rule。

完整清单见 `DO_NOT_REUSE.md`。

### 2.2 当前允许作为 scientific target 候选的结构

保留同一个 6420 OBC-DSGE 模型、prior、数据版本、`filter_R`、1964Q1 样本起点和同一经济 equilibrium-selection semantics。数值实现必须采用：

- 完整 8 维 Gaussian forecast score；
- 稳定 QR / square-root 线性代数，不显式求病态协方差逆；
- 与实际无约束线性 transition `F` 一致的 stationary initial covariance（历史 `F.T` 版本只作 robustness/replication control）；
- **单一联合** QMC/random-number design：同一个点的一行同时包含初始状态、全部结构创新、全部观测扰动；禁止三个独立 Sobol block 按行配对；
- endpoint 只切数据/随机流前缀；追加未来数据不能改变旧期原始随机数；
- OBC 求解不能以“超时=拒绝”改变 target；求解范围若扩大，必须保留相同经济选择规则并版本化。

最终 main `N` 不预设为 4096 或 8192；必须按 Phase A 认证后锁定。低保真 `N` 可以用于 proposal/screening，但最终 posterior 必须对锁定的 high-fidelity target 做 exact MH/DA correction，或者经过合法且通过 gate 的 importance correction。

## 3. 必须按顺序完成的阶段

### Phase A — Final target certification

目的：先认证“要采样的数值 posterior 是什么”，再长跑。

必须完成：

1. 联合随机数/QMC 构造回归测试：
   - 已知独立 Gaussian 积分；
   - 关键坐标联合投影；
   - powers-of-two prefix nesting；
   - endpoint prefix identity；
   - 同 seed 重建逐位一致。
2. score / filter 数值测试：
   - full Gaussian 单位等变性；
   - QR 与高精度参考；
   - 线性极限与 exact Kalman；
   - grouped transition 与逐粒子实现逐路径/flag 等价。
3. OBC branch/support audit：
   - 在后验相关区域而不是只在 archive mode 检查 fallback；
   - 若扩大 `(l_max,k_max)`，比较 density/economic-function 变化；
   - 正式 retained draws 中不得静默接受 solver failure。
4. ensemble-size / scramble sensitivity：
   - 候选 `N` 至少比较 `N` 与 `2N`；必要时继续增加；
   - 使用 posterior-relevant validation states，不能只用旧未收敛链的少数固定点；
   - 同时报 relative log density 与 headline economic function stability。
5. 将最终选择写入 `OUTPUT/target_contract.json`，包括模型/数据/prior/filter_R/initialization/score/QMC design/N/scramble/branch rule/source commits/hash。

**Phase A 未通过时，不允许启动“正式 production”。** 可以跑 pilot 来获得 posterior-relevant validation states，但必须标记 `VALIDATION_ONLY`。

### Phase B — 2007Q4 nonlinear posterior

这是第一个硬门槛。

- estimation endpoint `a=2007Q4`；
- 使用 Phase A 锁定的 high-fidelity nonlinear target；
- 可以使用 surrogate、低 N、linear posterior、DE bank、transport proposal 等作为 proposal/screening；
- 任何 delayed-acceptance 最终接受都必须包含 exact high-fidelity correction；
- adaptation 只允许在 warmup；production 开始前冻结；warmup 全丢弃；
- 至少 4 条独立 cold chains；可以采用 replica exchange/tempering；
- checkpoint 至少每 5–10 sweeps/固定批次保存一次，并保存 RNG state；
- 每次 resume 先重算当前状态 exact density；
- 不能因为工具/作业超时把 proposal 记为 reject；作业异常必须从最后原子 checkpoint 重跑。

正式 gate 见 `ACCEPTANCE_GATES.json`。若混合失败，可以重构 proposal，但所有 tuning trace 作废后从 frozen kernel 重新计 production。

### Phase C — 2009Q4 and 2019Q4 nonlinear posteriors

2007Q4 通过后，复制同一科学 target 和预先写好的 tuning policy到：

- `a=2009Q4`
- `a=2019Q4`

proposal 参数可以在各 endpoint 的 warmup 中自适应，但 adaptation policy 必须相同且不看经济结果；production 仍从冻结后重新计数。

三个端点都必须单独通过 gate。不得因某 endpoint 难混而只报告另外两个。

### Phase D — High-fidelity correction / sensitivity

如果正式生产直接以 high-fidelity exact target 校正，则记录 DA second-stage diagnostics；无需再把重要性重权作为主路径。

如果采用 source posterior + importance/PSIS bridge：

- source posterior 必须先独立收敛；
- 在 retained draws 上计算 exact high-fidelity log density；
- 报 Pareto-k、raw/PSIS weight ESS、MCMC autocorrelation-aware precision；
- 对所有 headline functions 报重权 MCSE；
- 若 gate 不过，必须切换 exact delayed acceptance / SMC / 其他严格 target-preserving 方法，不能继续“修权重”直到好看。

### Phase E — Conditional smoothing and factual replay

主比较定义：

- 参数估计端点 `a ∈ {2007Q4, 2009Q4, 2019Q4}`；
- 主 smoothing information endpoint 固定 `b=2019Q4`；
- 危机创新窗口：`2007Q4–2009Q1`（6 个季度）；
- 短期结果窗口：`2008Q1–2009Q2`（6 个季度）；
- persistence 结果延伸至 `2016Q4`。

对于**每一个 posterior draw**，重新在同一 theta 下提取 inherited state 和 shock history。禁止旧 residual/new theta 交叉配对。

必须保存：

- theta；
- inherited state；
- smoothed shocks；
- smoothing RNG/selector metadata；
- factual replay diagnostics；
- OBC solver flags；
- path/continuation validity。

如果 smoothing 本身是随机的，必须用预设独立 inner draws 或证明所采用 deterministic selector 的解释边界；不得把单一路径偷偷解释成已对 smoothing uncertainty 积分。

### Phase F — New endpoint attribution

对于每个合格 posterior/smoothing record：

1. 从实际 `mod.shocks` 读取玩家列表，禁止硬编码 7/8；预期 8，但以实际为准。
2. 短期窗口运行全部 `2^K` coalition；若 K=8，则 256。
3. OBC 与 linear benchmark 使用相同 theta、inherited state、shock path、future innovations，唯一差别为约束求解。
4. 保存全部 coalition values 与 solver validity。
5. 重新计算并交叉验证：
   - Shapley（边际贡献公式）；
   - Möbius/Harsanyi dividend；
   - Shapley-from-dividends；
   - efficiency identity。
6. headline outputs 至少：
   - total GDP loss OBC；
   - total GDP loss linear；
   - ELB/OBC amplification；
   - GDP `e_u` Shapley；
   - GDP `e_i` Shapley；
   - `e_u × e_i` second-order Harsanyi；
   - total nonadditivity；
   - consumption `e_u`；
   - investment `e_i`、`e_u`；
   - risk-premium rank / rank probability（仅合格 posterior）；
   - ELB binding duration attribution to 2016Q4。
7. persistence 游戏若利用 null player reduction，必须先在当前新模型记录中证明该 player 对该 functional 为数值 null；不能因为旧 P/F 是 128 coalition 就沿用。
8. 任何失败记录不得填零。若存在 failure，先修 solver；若仍存在，报告 failure-selection audit，主结果不得静默条件于成功。

### Phase G — Learning decomposition

主三端点 trajectory：固定 `b=2019Q4`，比较 `a=2007Q4,2009Q4,2019Q4`。

两因素 decomposition 使用合法的完整危机 smoothing endpoint：

- `a0=2007Q4`, `a1=2019Q4`；
- `b0=2009Q4`, `b1=2019Q4`。

计算：

`F00 = F_{a0,b0}`  
`F10 = F_{a1,b0}`  
`F01 = F_{a0,b1}`  
`F11 = F_{a1,b1}`

对于每个 headline function 报**期望值**分解：

`D_theta = 0.5 * [(F10-F00) + (F11-F01)]`

`D_eta   = 0.5 * [(F01-F00) + (F11-F10)]`

必须验证：

`D_theta + D_eta = F11 - F00`

并单独报告：

`interaction = F11 - F10 - F01 + F00`

注意 interaction 已通过对称路径进入两项平均效应的路径依赖比较；不要把 interaction 再加一次。中位数之差不满足上述加总恒等式，不能拿中位数代替均值分解。

### Phase H — Final economic robustness

至少完成：

- main N vs higher-N sensitivity；
- main QMC scramble/randomization sensitivity；
- branch-range sensitivity；
- initial-covariance robustness（scientific F-stationary 为主，historical F.T 为明确 robustness/replication control）；
- historical pseudodensity protocol 只作 replication contrast；
- headline-function MCSE；
- draw-window / chain stability；
- failure-selection audit。

只有经济函数稳健后，target certification 才最终闭环。

### Phase I — Nonlinear simulate–estimate–recover（高水平加固）

主实证闭环后，再执行至少一个可审计的 nonlinear recovery exercise。目标不是证明所有 Bayes 计算完美，而是检查：已知 DGP 下，endpoint learning pipeline 是否能在合理 MC error 内恢复 narrative revision 的方向和 parameter-vs-smoothing decomposition。

如计算资源有限，可以先做固定真参数、多个模拟数据集、三个 endpoint 的精简版；必须明确规模和局限。不要把现有 Gaussian/QMC 单元测试改名为 simulate-estimate-recover。

## 4. 最终停止条件

只有满足 `ACCEPTANCE_GATES.json` 且以下输出全部非 null，才能在 `OUTPUT/final_status.json` 写：

```json
{"paper_ready": true}
```

必须同时有：

- 3/3 qualified nonlinear posteriors；
- target certification；
- high-fidelity correction evidence；
- smoothing/factual replay；
- new three-endpoint attribution；
- learning decomposition；
- economic robustness/MCSE；
- reproduction manifest。

如果某项科学上失败，**不要为了完成任务降低 gate**。应输出：

```json
{"paper_ready": false, "blocking_stage": "...", "blocking_evidence": "..."}
```

并尽最大努力解决阻塞项；只有无法在现有模型/数据/资源下解决时才停止。

## 5. Codex 的运行纪律

- 优先使用长时持续作业，而不是交互式 1–2 分钟调用；后台作业必须有日志、PID/job ID 和 checkpoint。
- 自动恢复只能从原子 checkpoint；失败的半步不得记入 trace。
- 每次重大阶段结束后生成机器可读 JSON + 人可读 Markdown。
- 任何新的 bug/方法改变必须版本化，旧输出隔离，不覆盖。
- 所有 random seeds、source commit、package hash、Python/package version 必须记录。
- 所有数值表必须能够由原始数组重新生成。
- 不提前查看/优化 risk-premium 的最终排名来决定 sampling policy。

## 6. 最终交付

最终至少交付：

- `OUTPUT/final_status.json`
- `OUTPUT/PREWRITING_FINAL_REPORT_CN.md`
- `OUTPUT/target_contract.json`
- `OUTPUT/target_certification.json`
- `OUTPUT/posterior_diagnostics.csv`
- `OUTPUT/endpoint_main_results.csv`
- `OUTPUT/learning_decomposition.csv`
- `OUTPUT/economic_robustness.csv`
- `OUTPUT/factual_replay_diagnostics.csv`
- 三端点 posterior draws/checkpoints；
- 三端点 smoothing records；
- coalition raw arrays and validity ledger；
- all logs；
- `OUTPUT/SHA256SUMS.json`
- 可从零重现核心表的单命令脚本。

详细字段见 `OUTPUT_REQUIREMENTS.md`。
