"""V3 interval bridge using the original source's E-projected next state."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools" / "pyflint"), str(ROOT / "scripts"), str(ROOT / "methods")]

import numpy as np
from flint import arb, arb_mat, ctx, fmpq
from threadpoolctl import threadpool_limits

from auxiliary_v4_actual_arb_krawczyk_pilot_v2 import (
    matrix_from_sparse, numeric_pencil, qfloat, select_rows,
    sha, unserial,
)
from auxiliary_joint_qmc_gaussian_v4 import AuxiliaryJointQMCGaussianV4
from auxiliary_v4_actual_graph_pilot_v1 import raw_pencil
import locked_author_likelihood as author
from auxiliary_v4_source_primitive_arb_v4 import serial
from independent_multispell_milp_v1 import load_capture
from independent_multispell_uncapped_farkas_v3 import HIGH_DENSITY

MODE = os.environ.get("OBC_EAWARE_GREEN_MODE", "analytic_exact")
RADIUS_DEN = int(os.environ.get("OBC_EAWARE_GREEN_RADIUS_DEN", "1000000000000000000000000"))
if MODE not in ("analytic_exact", "analytic_theta39_box_exp120", "frozen_binary64"):
    raise ValueError("Unsupported semantic mode")
SOURCE = ROOT / "evidence" / ("auxiliary_v4_source_primitive_arb_analytic_theta39_box_exp120_v4"
                              if MODE == "analytic_theta39_box_exp120" else
                              "auxiliary_v4_source_primitive_arb_analytic_exact_v2")
ARCHIVE = SOURCE / "primitive_interval_sparse.json.gz"
SOURCE_REPORT = SOURCE / "report.json"
ROWS = ROOT / "evidence" / "auxiliary_v4_raw_row_elimination_control_v1.json"
KRAWCZYK = ROOT / "evidence" / f"auxiliary_v4_eaware_graph_arb_{MODE}_radius1over{RADIUS_DEN}_v4.json"
INPUT = ROOT / "evidence" / "auxiliary_joint_qmc_2009_stress_v4" / "high_density_flag_e625_seed456789_N8192_flagged_inputs.npz"
OUT = ROOT / "evidence" / f"auxiliary_v4_eaware_green_arb_{MODE}_radius1over{RADIUS_DEN}_v4"
LOG = ROOT / "logs" / "auxiliary_v4_eaware_green_arb_v4.log"
H = 24
Q, P = 33, 31


def atomic_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def atomic_gzip_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6) as f:
            f.write(json.dumps(obj, separators=(",", ":"), allow_nan=False).encode())
        raw.flush()
        os.fsync(raw.fileno())
    os.replace(tmp, path)


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat()} {message}\n")
        f.flush()
        os.fsync(f.fileno())


def block(M: arb_mat, rows: range, cols: range) -> arb_mat:
    out = arb_mat(len(rows), len(cols))
    for i, ri in enumerate(rows):
        for j, cj in enumerate(cols):
            out[i, j] = M[ri, cj]
    return out


def row_sum_norm_upper(M: arb_mat) -> arb:
    best = arb(0)
    for i in range(M.nrows()):
        row = arb(0)
        for j in range(M.ncols()):
            row += abs(M[i, j]).upper()
        best = max(best, row.upper())
    return best.upper()


def row_response(M: arb_mat, j: int) -> arb_mat:
    out = arb_mat(1, M.ncols())
    for k in range(M.ncols()):
        out[0, k] = M[j, k]
    return out


def pack_mat(M: arb_mat) -> dict:
    return {"shape": [M.nrows(), M.ncols()],
            "entries": [serial(M[i, j]) for i in range(M.nrows()) for j in range(M.ncols())]}


def main() -> None:
    threadpool_limits(1)
    ctx.prec, ctx.threads = 256, 1
    start = time.monotonic()
    log(f"START E-aware source primitive to Green interval bridge mode={MODE}")
    sr = json.loads(SOURCE_REPORT.read_text(encoding="utf-8"))
    kr = json.loads(KRAWCZYK.read_text(encoding="utf-8"))
    rr = json.loads(ROWS.read_text(encoding="utf-8"))
    if sr["interval_archive_sha256"] != sha(ARCHIVE):
        raise RuntimeError("Source archive checksum mismatch")
    if not kr["krawczyk_inclusion_verified"]:
        raise RuntimeError("Matching source graph certificate missing")
    if kr["semantic_mode"] != MODE or kr["radius"] != f"1/{RADIUS_DEN}" or not kr["E_projection_explicit"]:
        raise RuntimeError("Unexpected graph certificate mode or radius")
    with gzip.open(ARCHIVE, "rt", encoding="utf-8") as f:
        data = json.load(f)
    d, _ = load_capture(HIGH_DENSITY)
    if MODE in ("analytic_exact", "analytic_theta39_box_exp120"):
        if kr["source_sha256"]["source_archive"] != sha(ARCHIVE):
            raise RuntimeError("Source certificate archive hash mismatch")
        Praw, Nraw = matrix_from_sparse(data["Praw"]), matrix_from_sparse(data["Nraw"])
        Pd, Nd = numeric_pencil(Praw, Nraw)
        xbar = unserial(data["xbar_interval"])
    else:
        author.META = Path("O:/RUN/jw_win_v1/source/OBC_CLOSURE_AUDIT_AND_REPAIR_PACKAGE/"
                           "workspace/obc_repro/OBC_frontier_run/reference/previous/raw/"
                           "rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz")
        target = AuxiliaryJointQMCGaussianV4(N=2, seed=456789, l_max=4, k_max=32)
        system = target._system(target.full_from_estimated(np.asarray(d["theta"], float)))
        Pd, Nd, *_ = raw_pencil(target, system)
        Praw, Nraw = arb_mat(64, 64), arb_mat(64, 64)
        for i in range(64):
            for j in range(64):
                Praw[i, j], Nraw[i, j] = qfloat(Pd[i, j]), qfloat(Nd[i, j])
        xbar = qfloat(float(d["x_bar"]))
    radius = arb(fmpq(1, RADIUS_DEN))
    omega = arb_mat(P, Q)
    center_archive = Path(kr["center_archive"])
    if kr["center_archive_sha256"] != sha(center_archive):
        raise RuntimeError("Certified E-aware Arb center checksum mismatch")
    with gzip.open(center_archive, "rt", encoding="utf-8") as f:
        center_data = json.load(f)
    if center_data["shape"] != [P, Q]:
        raise RuntimeError("Certified center shape mismatch")
    for i in range(P):
        for j in range(Q):
            omega[i, j] = unserial(center_data["entries"][i*Q+j]) + arb(0, radius)

    bottom = rr["bottom_row_indices"]
    top = rr["top_row_indices"]
    if len(bottom) != P or len(top) != Q or sorted(bottom+top) != list(range(64)):
        raise RuntimeError("Bad row permutation")
    if top[-1] != 63 or 63 in bottom:
        raise RuntimeError("Wedge row convention invalid")
    Pb, Nb = select_rows(Praw, bottom), select_rows(Nraw, bottom)
    Pt, Nt = select_rows(Praw, top), select_rows(Nraw, top)
    D = block(Nb, range(P), range(Q, 64))
    Dinv = D.inv()
    Pbot = Dinv*Pb
    Nbot = Dinv*Nb
    factor = block(Nt, range(Q), range(Q, 64))
    Ptop = Pt-factor*Pbot
    Ntop = Nt-factor*Nbot
    S = arb_mat(64, 64)
    T = arb_mat(64, 64)
    for i in range(64):
        for j in range(64):
            S[i, j] = Ptop[i, j] if i < Q else Pbot[i-Q, j]
            T[i, j] = Ntop[i, j] if i < Q else Nbot[i-Q, j]
    structural_identity_enclosed = (all(T[i, Q+j].contains(0) for i in range(Q) for j in range(P))
                                    and all(T[Q+i, Q+j].contains(int(i == j))
                                            for i in range(P) for j in range(P)))
    if not structural_identity_enclosed:
        raise RuntimeError("Raw row elimination structural identity excluded")
    # Forcing is -e_(raw policy row 63), and the row partition leaves this
    # row last in the top block, so no inexact row transform enters u.
    u_q = arb_mat(Q, 1)
    u_q[Q-1, 0] = -1
    u_p = arb_mat(P, 1)
    Sqq = block(S, range(Q), range(Q))
    Sqp = block(S, range(Q), range(Q, 64))
    Spq = block(S, range(Q, 64), range(Q))
    Spp = block(S, range(Q, 64), range(Q, 64))
    Tqq = block(T, range(Q), range(Q))
    Tpq = block(T, range(Q, 64), range(Q))
    E = arb_mat(Q, Q)
    for i in range(25): E[i, i] = 1
    A = Sqq+Sqp*omega*E
    Ainv = A.inv()
    lam = Ainv*Tqq
    G = -(Ainv*Sqp)
    L = -(Ainv*u_q)
    left = Spq+Spp*omega*E
    B = left*G+Spp
    a = left*L+u_p
    F = E*lam
    source_top_residual = A*lam-Tqq
    source_bottom_residual = left*lam-Tpq-omega
    source_residual_contains_zero = all(source_top_residual[i,j].contains(0)
                                        for i in range(Q) for j in range(Q)) and all(
                                        source_bottom_residual[i,j].contains(0)
                                        for i in range(P) for j in range(Q))
    if not source_residual_contains_zero:
        raise RuntimeError("E-aware raw-equivalent source graph residual excludes zero")

    report = {"scope": "VALIDATION_ONLY_E_AWARE_SOURCE_GREEN_INTERVAL_V4",
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "source_sha256": {"runner": sha(Path(__file__)), "source_archive": sha(ARCHIVE),
                                "source_report": sha(SOURCE_REPORT), "krawczyk": sha(KRAWCZYK),
                                "row_partition": sha(ROWS), "incoming": sha(INPUT)},
              "semantic_mode": ("analytic_exact_decimal_YAML_binary64_theta" if MODE == "analytic_exact"
                                else "analytic_exact_decimal_YAML_theta39_box_exp120" if MODE == "analytic_theta39_box_exp120"
                                else "frozen_implementation_binary64_primitive"),
              "krawczyk_inclusion_verified": True,
              "E_projection_explicit": True,
              "raw_row_elimination_structural_identity_enclosed": structural_identity_enclosed,
              "raw_equivalent_source_graph_residual_contains_zero": source_residual_contains_zero,
              "certified_center_archive_sha256": kr.get("center_archive_sha256"),
              "green_intervals_constructed": False,
              "power_contraction_verified": False,
              "source_binary64_semantics_included": MODE == "frozen_binary64",
              "dual_sign_certified": False,
              "scientific_impossibility_proven": False,
              "stage": "GREEN_BLOCKS_INTERVAL_COMPLETE"}
    atomic_json(OUT / "report.json", report)
    log("GREEN_BLOCKS_INTERVAL_COMPLETE")
    with np.load(INPUT, allow_pickle=False) as z:
        qin = z["auxiliary_inputs"][0].copy()
    if hashlib.sha256(np.ascontiguousarray(qin).tobytes()).hexdigest() != "418ef31d9f7bebe39f5787dc6bceabf24d7eb96583e871ccaddd7563ffa1a8d5":
        raise RuntimeError("Incoming vector checksum mismatch")
    x0 = arb_mat(Q, 1)
    for i in range(Q): x0[i, 0] = qfloat(qin[i])
    qvec = []
    x = x0
    for t in range(H):
        qvec.append((lam*x)[4, 0]-xbar)
        x = F*x
    M = arb_mat(H, H)
    for s in range(H):
        b = [arb_mat(P, 1) for _ in range(H+1)]
        for t in range(s, -1, -1):
            b[t] = B*b[t+1]+(a if t == s else arb_mat(P, 1))
        x = arb_mat(Q, 1)
        for t in range(H):
            y = lam*x+G*b[t+1]+(L if t == s else arb_mat(Q, 1))
            M[t, s] = y[4, 0]+(1 if t == s else 0)
            x = E*y
    payload = {"semantics": report["semantic_mode"], "precision_bits": 256,
               "source_sha256": report["source_sha256"],
               "omega": pack_mat(omega), "lambda": pack_mat(lam), "F": pack_mat(F),
               "G": pack_mat(G), "L": pack_mat(L), "B": pack_mat(B),
               "a": pack_mat(a), "M24": pack_mat(M),
               "q24": [serial(v) for v in qvec], "xbar": serial(xbar)}
    atomic_gzip_json(OUT / "green_interval_matrices.json.gz", payload)
    report["green_interval_archive_sha256"] = sha(OUT / "green_interval_matrices.json.gz")
    report["green_intervals_constructed"] = True
    report["stage"] = "M24_Q24_INTERVAL_COMPLETE"
    report["M24_max_radius_float"] = max(float(M[i, j].rad()) for i in range(H) for j in range(H))
    report["q24_max_radius_float"] = max(float(v.rad()) for v in qvec)
    atomic_json(OUT / "report.json", report)
    log("M24_Q24_INTERVAL_COMPLETE")
    # Matrix powers enclose all large powers via submultiplicativity.
    # This is a finite exact sequence of interval matrix multiplications.
    contraction = {}
    for name, base, powers in (("F", F, (16, 32, 64, 128, 256)),
                               ("B", B, (16, 32, 64, 128))):
        square = base
        n = 1
        norms = {}
        while n < max(powers):
            square = square*square
            n *= 2
            if n in powers:
                norms[str(n)] = float(row_sum_norm_upper(square))
                report[f"{name}_power_norm_upper"] = norms
                report["stage"] = f"{name}_POWER_{n}_CHECKPOINT"
                atomic_json(OUT / "report.json", report)
                log(f"{name}_POWER_{n} norm={norms[str(n)]}")
        contraction[name] = norms
    report["power_norm_upper"] = contraction
    report["power_contraction_verified"] = any(v < 1 for v in contraction["F"].values()) and any(v < 1 for v in contraction["B"].values())
    report["stage"] = "COMPLETE_INTERVAL_BRIDGE"
    report["elapsed_seconds"] = time.monotonic()-start
    atomic_json(OUT / "report.json", report)
    log(f"COMPLETE contraction={report['power_contraction_verified']}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
