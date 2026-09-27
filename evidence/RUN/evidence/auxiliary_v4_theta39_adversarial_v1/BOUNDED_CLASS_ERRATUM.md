# Bounded-path scope correction to the frozen V4 adversarial report

This append-only erratum narrows one sentence in `REPORT.md`; that report and
its manifest remain unchanged.

The V4 row-sum and Fubini calculation applies to the **Green branch**
`d_t=Σ_(k≥0) B^k a λ_(t+k)`, with bounded `λ`, or to any original path for which
`B^n d_(t+n)→0` for each `t`. A bounded original `d` suffices because
`||B^64||∞<1`. Then `d` is bounded by the summable kernel and, for finite
initial `x`, the forward `F` recursion yields bounded `x`. Bounded `λ` by
itself does not force a general solution of `d_t=B d_(t+1)+a λ_t` onto this
branch: an unbounded homogeneous backward component may persist. The V4
report's phrase “also proves bounded d and x ... given finite initial x” must
be read with the Green-branch/vanishing-tail condition.

The finite-support obstruction to the locked `find_lk` eventually-slack
selector is unaffected. Its selected path imposes the all-slack terminal
branch, and the all-future dual excludes every finite-support nonnegative
policy wedge satisfying nonnegative slack on the certified positive-mass
initial-state/current-shock event. No claim about every possible unbounded
equilibrium path follows.
