# V4 theta-box initial-law orientation erratum

Status: validation-only correction. The frozen V4 code and reports remain
unchanged; this note governs their interpretation.

`auxiliary_v4_theta39_source_stationary_support_arb_v4.json` proves support for
the **scientific main** initial law specified in `START_HERE_CN.md` and
`CODEX_MASTER_TASK.md`: the stationary covariance of the actual physical
transition `F=qm[1,0][:25,:25]`, with controllability columns
`[E,FE,F^2E,F^3E]`. The source row proof in that report also verifies
`qm[1,0]=Lambda_E` over the theta box. These algebraic and support results
remain valid.

The report's scope phrase `ORIGINAL_SOURCE_INITIAL_LAW`, its
`source_stationary_gaussian_support_intersects_33d_input_ball` key, and the
V4 protocol's description of an “original preprocess initial law” should
**not** be read as validating the historical executable's covariance
orientation. `locked_author_likelihood.py` explicitly computes
`P=solve_discrete_lyapunov(F.T,E diag(sigma^2) E.T)` before drawing the initial
state, whereas its transition updates the state with `F @ x`. The historical
covariance support is spanned by `[(F.T)^k E]`; V4's chosen 21 columns do not
cover the frozen input under that law (point least-squares L-infinity distance
about `3.285`).

The attached task documents designate `F`-stationary covariance as scientific
main and historical `F.T` as robustness/replication control. Therefore the V4
support certificate applies to the requested scientific main law. A separate
versioned V5 interval certificate is required before making a historical
`F.T` support claim. No endpoint or impossibility flag is changed by this
erratum.
