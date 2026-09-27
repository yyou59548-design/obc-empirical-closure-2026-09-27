# V3 frozen-binary64 Green provenance erratum

Status: validation-only metadata correction; no V3 script, archive, or report
has been overwritten.

The V3 `frozen_binary64` Green report at
`auxiliary_v4_eaware_green_arb_frozen_binary64_radius1over1000000000000000000000000_v3/report.json`
records `source_sha256.source_archive` equal to the V2 **analytic** primitive
archive hash `d2313857acd06e6249b0f1fa088d81cc3decf03de7136a4496f6796fa5837ac7`.
That field is a metadata path selected at module import, not the actual
binary64 coefficient input in this branch. Line 133 onward in frozen script
`auxiliary_v4_eaware_green_arb_v3.py` instead reconstructs `Pd,Nd` with
`raw_pencil(target, system)` and encloses each float with `qfloat` (exact
dyadic). The graph `frozen_binary64` branch likewise calls `raw_pencil` and
checks its primitive rows against the captured primitive matrices. The V2
implementation-hull archive is used for the graph's checksum context; neither
archive is the Green branch's coefficient operator.

The actual binary64 reconstruction code and input hashes as of this erratum:

| Input | SHA256 |
| --- | --- |
| `auxiliary_v4_eaware_graph_arb_v3.py` | `1cb6ef3a5d7561e0244b65a851b9eb8dabbbea258db9cbb8fb62630d4829d0b0` |
| `auxiliary_v4_eaware_green_arb_v3.py` | `dc3f8158603319e223763e544ce5ea3267ce8e2ca5b7e02c877b1e4b8b7701cf` |
| `auxiliary_v4_actual_graph_pilot_v1.py` (`raw_pencil`) | `acc524776180cfe1a2385f590652a660bafde8af6940adad7d39dfce16bb9f8f` |
| `auxiliary_joint_qmc_gaussian_v4.py` | `2a61064ddb3a384ec8cb198aeeb1ef92f492665d62ad0ed5ec1e0eeaaf383470` |
| `locked_author_likelihood.py` | `fd89dfe000698fc5ce2c695ab91e5292b928e0ad4369835483f9073b76882dc5` |
| `model_engine.py` | `44e746289575b4612bfc84baf8207eeb4ba30647511b87d5927596ec8432297a` |
| frozen high-density capture | `0951396abd8cd091dd6267c836b4121ab170b13025eff0f300f927da354b5a2a` |
| source meta NPZ with YAML | `b371d3016bd46daa2de99227d80c6cef0671e93221e508ff62e2325fa3cbd3e6` |

This erratum must accompany any V3 frozen-binary64 report citation. The
independent V4 positive-theta-box analytic chain has explicit source archive
hashes and is unaffected. No endpoint status is asserted here.
