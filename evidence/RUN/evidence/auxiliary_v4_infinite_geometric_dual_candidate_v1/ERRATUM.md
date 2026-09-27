# Candidate V1 output-label correction

The JSON field `asymptotic_normalized_column_coefficient_numeric` is
mislabelled. The stored value `-2.0531141792775816` is the limiting
coefficient of `c_s/y_s`, where `y_s=alpha*rho^(s-192)` and
`alpha=7.224524343770655e-05`. The limit for the JSON's stated
normalization `c_s/rho^(s-192)` is therefore approximately
`-1.483277e-4`. The sign, all stored finite-scan columns, and the
`q_y_infinite_numeric` value are unaffected. This numerical output is
still **not** an interval certificate or original-model impossibility
proof. The V1 source/manifest/candidate are preserved unchanged.
