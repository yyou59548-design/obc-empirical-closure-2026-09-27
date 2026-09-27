# CODEX MASTER TASK — Complete the OBC endpoint-learning empirical pipeline before manuscript drafting

## Mission

Carry this project from the current audited state to a scientifically defensible **paper-ready empirical closure**. Do not draft or polish the manuscript until the empirical gates are passed. The final deliverable is a reproducible evidence package, not a persuasive narrative.

Current truth: `paper_ready=false`; qualified nonlinear endpoint posteriors = `0/3`.

Read `START_HERE_CN.md` for the full scientific design. This file is the compact execution contract.

## Authoritative inputs

1. `INPUTS/01_LATEST_PREWRITING_AUDIT.zip` — latest audit, code, linear reference runs, QMC investigations, downstream interface audit.
2. `INPUTS/02_JOINT_QMC_CHECKPOINT_AND_AUDIT.zip` — repaired joint-QMC implementation and nonlinear pilot checkpoint. Its trace is **validation only** unless independently re-qualified; use locations/proposal information only as allowed.

Historical filenames containing “final”, “production”, or “certified” are not evidence. Recompute gates.

## Scientific target candidate

Same 6420 nonlinear OBC DSGE model/prior/data/filter_R; sample starts 1964Q1. Use:

- full 8-dimensional Gaussian forecast score;
- numerically stable QR/square-root algebra;
- stationary initial covariance for the actual linear transition `F` as the scientific main specification;
- a **single joint** random/QMC design covering initial state + all structural shocks + all observation perturbations; never separately generate blocks and pair rows;
- endpoint-prefix slicing only;
- fixed equilibrium-selection semantics; expand numerical search only transparently and with sensitivity evidence;
- no wall-clock timeout-as-rejection.

The high-fidelity ensemble size N is **not yet certified**. Certify it before final production. Low-fidelity likelihoods may be used for proposal/screening only.

## Mandatory stage order

A. Target certification  
B. Qualified 2007Q4 nonlinear posterior  
C. Qualified 2009Q4 and 2019Q4 nonlinear posteriors  
D. Exact high-fidelity correction / sensitivity  
E. Theta-conditional smoothing and factual replay  
F. Exact coalition attribution  
G. Parameter-vs-smoothing learning decomposition  
H. Economic-function robustness and MCSE  
I. Nonlinear simulate-estimate-recover (top-tier robustness)

Do not jump to F/G using historical P/F shocks or parameters.

## Sampling rules

- Minimum 4 independent cold chains.
- Adaptation/tuning is warmup only and discarded.
- Production kernel frozen before draw 1.
- Replica exchange, delayed acceptance, DE, transport, surrogate screening, SMC, or other kernels are allowed if they preserve the locked exact target.
- Surrogates may screen; they may not replace the exact target.
- Every resume re-evaluates current-state exact density.
- A crashed/timed-out proposal is not automatically a rejection; rerun from last atomic checkpoint.
- If a sampler is redesigned, old production becomes pilot; do not splice incompatible kernels unless the transition/adaptation protocol makes this mathematically valid and is explicitly documented.

## Posterior inference gate

Use `ACCEPTANCE_GATES.json`. At minimum, for each endpoint:

- max rank-normalized split R-hat < 1.01;
- minimum bulk ESS >= 400;
- minimum tail ESS >= 400;
- full / last 2/3 / last 1/2 consistency;
- multimodal/round-trip evidence if tempering is used;
- exact current-state cache revalidation;
- headline economic-function MCSE acceptable once those functions exist.

Do not delete slow parameters from diagnostics.

## Empirical information design

Primary trajectory:

- parameter endpoint `a ∈ {2007Q4, 2009Q4, 2019Q4}`;
- common smoothing information endpoint `b=2019Q4`.

Crisis innovations: 2007Q4–2009Q1.  
Short outcome window: 2008Q1–2009Q2.  
Persistence outcome: through 2016Q4.

For every theta draw, re-extract states/shocks **conditional on that same theta**. Never combine a new theta with archived residuals.

## Attribution

Use actual `mod.shocks` to define K players. Run all `2^K` coalitions for short-window outputs. Use the exact same theta, inherited state, shock paths, and future innovations for OBC and linear benchmark; only the constraint solution differs.

Must output at least:

- total GDP loss OBC / linear;
- OBC amplification;
- risk-premium (`e_u`) GDP Shapley;
- investment (`e_i`) GDP Shapley;
- `e_u × e_i` Harsanyi dividend;
- total nonadditivity;
- consumption `e_u`;
- investment `e_i`, `e_u`;
- risk-premium ranking/probability;
- ELB binding-duration attribution through 2016Q4.

Cross-check Shapley using both marginal-weight and Möbius/dividend formulas. Preserve full validity ledgers.

## Learning decomposition

Main trajectory holds b=2019Q4 fixed. For two-factor decomposition use:

- a0=2007Q4, a1=2019Q4;
- b0=2009Q4, b1=2019Q4.

Compute expectation cells F00/F10/F01/F11, then:

D_theta = 0.5 * ((F10-F00) + (F11-F01))
D_eta   = 0.5 * ((F01-F00) + (F11-F10))
interaction = F11-F10-F01+F00

Verify D_theta + D_eta = F11 - F00 numerically. Use expectations for additive decomposition; report medians/quantiles separately, never force median additivity.

## Anti-hallucination / anti-shortcut rules

- Never fill a missing result with zero.
- Never rename an archival result as a new endpoint result.
- Never call a stored-state count an ESS.
- Never call good swap acceptance convergence.
- Never call linear posterior nonlinear posterior.
- Never call a fixed-point likelihood grid posterior certification.
- Never claim exact original-model likelihood when using an EnKF/QMC approximation; state the actual numerical target.
- Never optimize analysis choices based on the sign/rank of the risk-premium result.

## Finish condition

Only set `OUTPUT/final_status.json -> paper_ready=true` when every mandatory gate in `ACCEPTANCE_GATES.json` is true and all headline result cells are non-null.

If a blocker remains, keep `paper_ready=false`, identify the blocker with quantitative evidence, leave resumable checkpoints, and do not write a favorable economic conclusion.
