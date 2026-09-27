"""V3 directed Arb adjoint for E-aware source Green and full input-ball margin.

VALIDATION_ONLY: no endpoint, paper, or scientific-impossibility flag is changed.
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import time
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools" / "pyflint"), str(ROOT / "scripts"), str(ROOT / "methods")]

import numpy as np
from flint import arb, arb_mat, ctx, fmpq
from threadpoolctl import threadpool_limits

from auxiliary_v4_actual_arb_krawczyk_pilot_v2 import qfloat, sha, unserial
from auxiliary_v4_eaware_green_arb_v4 import atomic_json, row_sum_norm_upper

MODE = os.environ.get("OBC_EAWARE_DUAL_MODE", "analytic_theta39_box_exp120")
DEN = 10**24
GREEN = ROOT / "evidence" / f"auxiliary_v4_eaware_green_arb_{MODE}_radius1over{DEN}_v4"
CANDIDATE = ROOT / "evidence" / "auxiliary_v4_infinite_geometric_dual_candidate_v1" / "candidate.json"
OUT = ROOT / "evidence" / f"auxiliary_v4_eaware_dual_adjoint_arb_{MODE}_v4.json"
LOG = ROOT / "logs" / "auxiliary_v4_eaware_dual_adjoint_arb_v4.log"
INPUT = ROOT / "evidence" / "auxiliary_joint_qmc_2009_stress_v4" / "high_density_flag_e625_seed456789_N8192_flagged_inputs.npz"
T, S = 192, 1024


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat()} mode={MODE} {message}\n")
        f.flush()
        os.fsync(f.fileno())


def aq(value: str | Fraction) -> arb:
    f = value if isinstance(value, Fraction) else Fraction(Decimal(str(value)))
    return arb(fmpq(f.numerator, f.denominator))


def mat(row: dict) -> arb_mat:
    n, m = row["shape"]
    if len(row["entries"]) != n*m:
        raise RuntimeError("Green matrix archive length mismatch")
    result = arb_mat(n, m)
    for i in range(n):
        for j in range(m): result[i, j] = unserial(row["entries"][i*m+j])
    return result


def eye(n: int) -> arb_mat:
    result = arb_mat(n, n)
    for i in range(n): result[i, i] = 1
    return result


def row(M: arb_mat, j: int) -> arb_mat:
    out = arb_mat(1, M.ncols())
    for k in range(M.ncols()): out[0, k] = M[j, k]
    return out


def one_norm_upper(v: arb_mat) -> arb:
    total = arb(0)
    for k in range(v.ncols()): total += abs(v[0, k]).upper()
    return total.upper()


def max_abs_upper(v: arb_mat) -> arb:
    return max(abs(v[k, 0]).upper() for k in range(v.nrows()))


def main() -> None:
    threadpool_limits(1)
    ctx.prec, ctx.threads = 256, 1
    start = time.monotonic()
    log("START E-aware directed dual")
    gr = json.loads((GREEN / "report.json").read_text(encoding="utf-8"))
    archive = GREEN / "green_interval_matrices.json.gz"
    if (gr["green_interval_archive_sha256"] != sha(archive) or not gr["power_contraction_verified"]
            or not gr["E_projection_explicit"] or not gr["raw_equivalent_source_graph_residual_contains_zero"]):
        raise RuntimeError("Green certificate missing or hash mismatch")
    if gr["semantic_mode"] != ("frozen_implementation_binary64_primitive" if MODE == "frozen_binary64"
                                else "analytic_exact_decimal_YAML_theta39_box_exp120" if MODE == "analytic_theta39_box_exp120"
                                else "analytic_exact_decimal_YAML_binary64_theta"):
        raise RuntimeError("Green semantic mode mismatch")
    with gzip.open(archive, "rt", encoding="utf-8") as f:
        d = json.load(f)
    cand = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    ytxt = cand["first_192_y_binary64_repr"]
    alpha_text = cand["alpha_binary64_repr"]
    if len(ytxt) != T or any(Fraction(Decimal(str(v))) < 0 for v in ytxt) or Fraction(Decimal(alpha_text)) <= 0:
        raise RuntimeError("Frozen dual rational positivity or length fails")
    y = [aq(s) for s in ytxt]
    alpha = aq(alpha_text)
    rho = aq(Fraction(199, 200))
    lam, F, G, L, B, a = (mat(d[k]) for k in ("lambda", "F", "G", "L", "B", "a"))
    xbar = unserial(d["xbar"])
    E = eye(33)
    for i in range(25, 33): E[i, i] = 0
    rl, rg = row(lam, 4), row(G, 4)
    invF = (eye(33)-rho*F).inv()
    pcoef = rl*invF
    # Binary powers avoid the exponential interval-dependency growth of a
    # 192-step forward/backward point recurrence on a wide-norm matrix.
    Fpowers = [F**k for k in range(T+1)]
    Bpowers = [B**k for k in range(T)]
    rlf = [rl*Fk for Fk in Fpowers]
    pcf = [pcoef*Fk for Fk in Fpowers]
    p = [arb_mat(1, 33) for _ in range(T+1)]
    p[T] = alpha*pcoef
    for t in range(T-1, -1, -1):
        term = alpha*pcf[T-t]
        for k in range(t, T): term += y[k]*rlf[k-t]
        p[t] = term
    hcoef = 1+L[4, 0]+(rho*pcoef*E*L)[0, 0]
    dcoef = rg+rho*pcoef*E*G
    with np.load(INPUT, allow_pickle=False) as z:
        qin = np.asarray(z["auxiliary_inputs"][0], float)
    import hashlib
    if hashlib.sha256(np.ascontiguousarray(qin).tobytes()).hexdigest() != "418ef31d9f7bebe39f5787dc6bceabf24d7eb96583e871ccaddd7563ffa1a8d5":
        raise RuntimeError("Incoming input hash mismatch")
    x0 = arb_mat(33, 1)
    for i in range(33): x0[i, 0] = qfloat(qin[i])
    qdot = arb(0)
    qgradient = arb_mat(1, 33)
    for t in range(T):
        qdot += y[t]*((rlf[t]*x0)[0, 0]-xbar)
        qgradient += y[t]*rlf[t]
    qdot += alpha*((rl*invF*Fpowers[T]*x0)[0, 0]-xbar/(1-rho))
    qgradient += alpha*rl*invF*Fpowers[T]
    gradient_l1 = one_norm_upper(qgradient)
    ball_radius = aq(Fraction(4, 10**7))
    qdot_ball_upper = qdot.upper()+gradient_l1*ball_radius
    result = {"scope": "VALIDATION_ONLY_E_AWARE_SOURCE_TO_GEOMETRIC_DUAL_ARB_V4",
              "semantic_mode": MODE,
              "source_sha256": {"runner": sha(Path(__file__)), "green_archive": sha(archive),
                                "green_report": sha(GREEN / "report.json"),
                                "candidate": sha(CANDIDATE), "incoming_archive": sha(INPUT)},
              "q_dot_y_upper": float(qdot.upper()),
              "q_dot_y_negative": bool(qdot < 0),
              "full_33d_input_gradient_l1_upper": float(gradient_l1),
              "input_ball_linf_radius": "4e-7",
              "q_dot_y_full_input_ball_upper": float(qdot_ball_upper),
              "full_input_ball_keeps_negative": bool(qdot_ball_upper < 0),
              "prefix_columns": S,
              "prefix_all_negative": False,
              "all_later_negative": False,
              "scientific_impossibility_proven": False,
              "stage": "Q_DOT_Y_COMPLETE"}
    atomic_json(OUT, result)
    log(f"Q_DOT_Y_COMPLETE upper={float(qdot.upper())} ball_upper={float(qdot_ball_upper)}")
    worst_upper = -float("inf")
    worst_s = -1
    dpre = [y[t]*rg+p[t+1]*E*G for t in range(T)]
    hpre = [y[t]*(1+L[4, 0])+(p[t+1]*E*L)[0, 0] for t in range(T)]
    gT = arb_mat(1, 31)
    for t in range(T): gT += dpre[t]*Bpowers[T-t-1]
    for s in range(T):
        gs = arb_mat(1, 31)
        for t in range(s): gs += dpre[t]*Bpowers[s-t-1]
        upper = float((hpre[s]+(gs*a)[0, 0]).upper())
        if upper > worst_upper: worst_upper, worst_s = upper, s
    result["stage"] = "PREFIX_192_CHECKPOINT"
    result["prefix_max_upper_so_far"] = worst_upper
    result["prefix_worst_s_so_far"] = worst_s
    atomic_json(OUT, result)
    vT = gT*(1/alpha)
    vstar = dcoef*(rho*eye(31)-B).inv()
    Cstar = hcoef+(vstar*a)[0, 0]
    A = B*(1/rho)
    vdiff = vT-vstar
    tail_prefix_max = -float("inf")
    tail_prefix_worst_s = -1
    for s in range(T, S):
        normalized = Cstar+(vdiff*(A**(s-T))*a)[0, 0]
        upper = float(normalized.upper())
        if upper > tail_prefix_max: tail_prefix_max, tail_prefix_worst_s = upper, s
        if s in (511, 1023):
            result["stage"] = f"PREFIX_{s+1}_CHECKPOINT"
            result["tail_normalized_max_upper_so_far"] = tail_prefix_max
            result["tail_normalized_worst_s_so_far"] = tail_prefix_worst_s
            atomic_json(OUT, result)
    result["prefix_max_upper"] = worst_upper
    result["prefix_worst_s"] = worst_s
    result["tail_normalized_prefix_max_upper"] = tail_prefix_max
    result["tail_normalized_prefix_worst_s"] = tail_prefix_worst_s
    result["prefix_all_negative"] = worst_upper < 0 and tail_prefix_max < 0
    Ak = eye(31)
    prefix_power_max = arb(0)
    for k in range(128):
        prefix_power_max = max(prefix_power_max, row_sum_norm_upper(Ak))
        Ak = Ak*A
    contraction = row_sum_norm_upper(Ak)
    if not contraction < 1:
        raise RuntimeError("B/rho power 128 did not contract")
    blocks = (S-T)//128
    tail_error = one_norm_upper(vdiff)*prefix_power_max*(contraction**blocks)*max_abs_upper(a)
    tail_upper = Cstar.upper()+tail_error.upper()
    result.update({"asymptotic_normalized_C_upper": float(Cstar.upper()),
                   "B_over_rho_power128_norm_upper": float(contraction),
                   "B_over_rho_prefix0to127_norm_max_upper": float(prefix_power_max),
                   "tail_blocks_lower_bound": blocks,
                   "tail_normalized_error_upper": float(tail_error.upper()),
                   "all_s_ge_1024_normalized_C_upper": float(tail_upper),
                   "all_later_negative": bool(tail_upper < 0),
                   "stage": "COMPLETE_DIRECTIONAL_INTERVAL_TEST",
                   "elapsed_seconds": time.monotonic()-start})
    atomic_json(OUT, result)
    log(f"COMPLETE prefix={result['prefix_all_negative']} tail={result['all_later_negative']}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
