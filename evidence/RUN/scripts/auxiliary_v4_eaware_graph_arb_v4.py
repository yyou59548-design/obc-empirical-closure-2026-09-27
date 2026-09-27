"""E-aware 1,023D Arb Krawczyk from the original source row equations.

The future state is E*qfull, where E zeroes the eight shock coordinates.
All output is validation-only and fail-closed.
"""
from __future__ import annotations

import gzip
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
from scipy.linalg import inv, solve
from threadpoolctl import threadpool_limits

from auxiliary_joint_qmc_gaussian_v4 import AuxiliaryJointQMCGaussianV4
from auxiliary_v4_actual_graph_pilot_v1 import raw_pencil
from auxiliary_v4_actual_arb_krawczyk_pilot_v2 import (
    atomic_json, atomic_gzip_json, matrix_from_sparse, numeric_pencil, qfloat,
    serial, sha, upper_abs,
)
from auxiliary_v4_source_green_arb_v2 import block
from independent_multispell_milp_v1 import load_capture
from independent_multispell_uncapped_farkas_v3 import HIGH_DENSITY
import locked_author_likelihood as author

MODE = os.environ.get("OBC_EAWARE_MODE", "analytic_theta39_box_exp120")
if MODE not in ("analytic_exact", "analytic_theta39_box_exp120", "frozen_binary64"):
    raise ValueError("Unsupported E-aware source mode")
SOURCE = ROOT / "evidence" / f"auxiliary_v4_source_primitive_arb_{MODE}_v4"
if MODE == "frozen_binary64":
    SOURCE = ROOT / "evidence" / "auxiliary_v4_source_primitive_arb_hull_implementation_v2"
ROWS = ROOT / "evidence" / "auxiliary_v4_raw_row_elimination_control_v1.json"
META = Path("O:/RUN/jw_win_v1/source/OBC_CLOSURE_AUDIT_AND_REPAIR_PACKAGE/"
            "workspace/obc_repro/OBC_frontier_run/reference/previous/raw/"
            "rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz")
DEN = int(os.environ.get("OBC_EAWARE_GRAPH_RADIUS_DEN", str(10**24)))
OUT = ROOT / "evidence" / f"auxiliary_v4_eaware_graph_arb_{MODE}_radius1over{DEN}_v4.json"
LOG = ROOT / "logs" / "auxiliary_v4_eaware_graph_arb_v4.log"
Q, P = 33, 31


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat()} mode={MODE} {message}\n")
        f.flush()
        os.fsync(f.fileno())


def select_rows(M: arb_mat, rows: list[int]) -> arb_mat:
    out = arb_mat(len(rows), M.ncols())
    for i, ri in enumerate(rows):
        for j in range(M.ncols()): out[i, j] = M[ri, j]
    return out


def interval_ST(Praw: arb_mat, Nraw: arb_mat, top: list[int], bottom: list[int]):
    Pb, Nb = select_rows(Praw, bottom), select_rows(Nraw, bottom)
    Pt, Nt = select_rows(Praw, top), select_rows(Nraw, top)
    D = block(Nb, range(P), range(Q, 64))
    Dinv = D.inv()
    Pbot, Nbot = Dinv*Pb, Dinv*Nb
    factor = block(Nt, range(Q), range(Q, 64))
    Ptop, Ntop = Pt-factor*Pbot, Nt-factor*Nbot
    S, T = arb_mat(64, 64), arb_mat(64, 64)
    for i in range(64):
        for j in range(64):
            S[i, j] = Ptop[i, j] if i < Q else Pbot[i-Q, j]
            T[i, j] = Ntop[i, j] if i < Q else Nbot[i-Q, j]
    return S, T


def np_ST(Praw: np.ndarray, Nraw: np.ndarray, top: list[int], bottom: list[int]):
    D = Nraw[np.ix_(bottom, range(Q, 64))]
    Pb, Nb = solve(D, Praw[bottom], check_finite=False), solve(D, Nraw[bottom], check_finite=False)
    factor = Nraw[np.ix_(top, range(Q, 64))]
    return np.vstack((Praw[top]-factor@Pb, Pb)), np.vstack((Nraw[top]-factor@Nb, Nb))


def np_HJ(S: np.ndarray, T: np.ndarray, omega: np.ndarray, E: np.ndarray):
    Sqq, Sqp, Spq, Spp = S[:Q, :Q], S[:Q, Q:], S[Q:, :Q], S[Q:, Q:]
    Tqq, Tpq = T[:Q, :Q], T[Q:, :Q]
    A = Sqq+Sqp@omega@E
    lam = solve(A, Tqq, check_finite=False)
    left = Spq+Spp@omega@E
    H = left@lam-Tpq-omega
    U = Spp-left@solve(A, Sqp, check_finite=False)
    J = np.kron(U, (E@lam).T)-np.eye(P*Q)
    return H, J


def arb_HJ(S: arb_mat, T: arb_mat, omega: arb_mat, E: arb_mat, need_J: bool):
    Sqq = block(S, range(Q), range(Q))
    Sqp = block(S, range(Q), range(Q, 64))
    Spq = block(S, range(Q, 64), range(Q))
    Spp = block(S, range(Q, 64), range(Q, 64))
    Tqq = block(T, range(Q), range(Q))
    Tpq = block(T, range(Q, 64), range(Q))
    A = Sqq+Sqp*omega*E
    Ainv = A.inv()
    lam = Ainv*Tqq
    left = Spq+Spp*omega*E
    H = left*lam-Tpq-omega
    if not need_J: return H, None
    U = Spp-left*Ainv*Sqp
    W = E*lam
    J = arb_mat(P*Q, P*Q)
    for i in range(P):
        for a in range(P):
            u = U[i, a]
            for k in range(Q):
                for b in range(Q):
                    J[i*Q+k, a*Q+b] = u*W[b, k]-(1 if i == a and k == b else 0)
    return H, J


def main() -> None:
    threadpool_limits(1)
    ctx.prec, ctx.threads = 256, 1
    start = time.monotonic()
    log("START")
    source_report = json.loads((SOURCE / "report.json").read_text(encoding="utf-8"))
    source_archive = SOURCE / "primitive_interval_sparse.json.gz"
    if source_report["interval_archive_sha256"] != sha(source_archive):
        raise RuntimeError("Source archive checksum mismatch")
    rr = json.loads(ROWS.read_text(encoding="utf-8"))
    top, bottom = rr["top_row_indices"], rr["bottom_row_indices"]
    if sorted(top+bottom) != list(range(64)) or top[-1] != 63 or 63 in bottom:
        raise RuntimeError("Raw row partition mismatch")
    with gzip.open(source_archive, "rt", encoding="utf-8") as f:
        source = json.load(f)
    d, _ = load_capture(HIGH_DENSITY)
    if MODE == "frozen_binary64":
        author.META = META
        target = AuxiliaryJointQMCGaussianV4(N=2, seed=456789, l_max=4, k_max=32)
        system = target._system(target.full_from_estimated(np.asarray(d["theta"], float)))
        Pd, Nd, AA, BB, CC, DD = raw_pencil(target, system)
        for name, M in (("A", AA[:55,:56]), ("B", BB[:55,:56]),
                        ("C", CC[:55,:56]), ("D", DD[:55])):
            if not np.array_equal(M, np.asarray(d["primitive_"+name], float)):
                raise RuntimeError("Frozen source primitive does not equal capture")
        Praw, Nraw = arb_mat(64, 64), arb_mat(64, 64)
        for i in range(64):
            for j in range(64):
                Praw[i, j], Nraw[i, j] = qfloat(Pd[i, j]), qfloat(Nd[i, j])
    else:
        if source_report["coefficient_mode"] != MODE:
            raise RuntimeError("Source semantics mismatch")
        Praw, Nraw = matrix_from_sparse(source["Praw"]), matrix_from_sparse(source["Nraw"])
        Pd, Nd = numeric_pencil(Praw, Nraw)
    S, T = interval_ST(Praw, Nraw, top, bottom)
    Sd, Td = np_ST(Pd, Nd, top, bottom)
    E, Ed = arb_mat(Q, Q), np.diag(np.r_[np.ones(25), np.zeros(8)])
    for i in range(25): E[i, i] = 1
    omega0 = np.asarray(d["recurrence_omega"], float)
    H0, J0 = np_HJ(Sd, Td, omega0, Ed)
    step = solve(J0, -H0.ravel(), check_finite=False)
    omega1 = omega0+step.reshape(P, Q)
    H1, J1 = np_HJ(Sd, Td, omega1, Ed)
    Cfloat = inv(J1, check_finite=False)
    Cmat = arb_mat(P*Q, P*Q)
    for i in range(P*Q):
        for j in range(P*Q): Cmat[i, j] = qfloat(Cfloat[i, j])
    report = {"scope": "VALIDATION_ONLY_E_AWARE_SOURCE_GRAPH_ARB_KRAWCZYK",
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "semantic_mode": MODE, "radius": f"1/{DEN}",
              "source_sha256": {"runner": sha(Path(__file__)), "source_archive": sha(source_archive),
                                "source_report": sha(SOURCE / "report.json"),
                                "row_partition": sha(ROWS), "capture": sha(HIGH_DENSITY)},
              "numeric_initial_residual_max_abs": float(np.max(np.abs(H0))),
              "numeric_newton_step_max_abs": float(np.max(np.abs(step))),
              "numeric_post_newton_residual_max_abs": float(np.max(np.abs(H1))),
              "E_projection_explicit": True,
              "stage": "NUMERIC_CENTER_COMPLETE", "krawczyk_inclusion_verified": False,
              "scientific_impossibility_proven": False}
    atomic_json(OUT, report)
    center = arb_mat(P, Q)
    for i in range(P):
        for j in range(Q): center[i, j] = qfloat(omega1[i, j])
    for k in range(3):
        Hc, _ = arb_HJ(S, T, center, E, False)
        hv = arb_mat(P*Q, 1)
        for i in range(P):
            for j in range(Q): hv[i*Q+j, 0] = Hc[i, j]
        corr = Cmat*hv
        for i in range(P):
            for j in range(Q): center[i, j] = (center[i, j]-corr[i*Q+j, 0]).mid()
        report[f"refinement_{k+1}_residual_upper"] = float(max(upper_abs(Hc[i,j]) for i in range(P) for j in range(Q)))
        report["stage"] = f"REFINEMENT_{k+1}_COMPLETE"
        atomic_json(OUT, report)
        log(report["stage"])
    center_path = OUT.with_suffix(".center.json.gz")
    atomic_gzip_json(center_path, {"shape": [P,Q], "entries": [serial(center[i,j]) for i in range(P) for j in range(Q)]})
    report["center_archive_sha256"] = sha(center_path)
    report["center_archive"] = str(center_path)
    radius = arb(fmpq(1, DEN))
    box = arb_mat(P, Q)
    for i in range(P):
        for j in range(Q): box[i,j] = center[i,j]+arb(0,radius)
    Hc, _ = arb_HJ(S, T, center, E, False)
    _, Jbox = arb_HJ(S, T, box, E, True)
    report["stage"] = "BOX_JACOBIAN_COMPLETE"
    atomic_json(OUT, report)
    log("BOX_JACOBIAN_COMPLETE")
    if time.monotonic()-start > 900:
        report["stage"] = "TIME_CAP_BEFORE_KRAWCZYK_PRODUCT"
        atomic_json(OUT, report)
        return
    hv = arb_mat(P*Q, 1)
    for i in range(P):
        for j in range(Q): hv[i*Q+j, 0] = Hc[i,j]
    CH, CJ = Cmat*hv, Cmat*Jbox
    b = max(upper_abs(CH[i,0]) for i in range(P*Q))
    eta = arb(0)
    for i in range(P*Q):
        rowsum = arb(0)
        for j in range(P*Q): rowsum += upper_abs((1 if i==j else 0)-CJ[i,j])
        eta = max(eta, rowsum.upper())
    condition = bool(b+eta*radius < radius)
    report.update({"preconditioned_residual_upper": float(b),
                   "preconditioned_jacobian_defect_norm_upper": float(eta),
                   "krawczyk_radius_condition": condition,
                   "krawczyk_inclusion_verified": condition,
                   "stage": "COMPLETE_INTERVAL_TEST", "elapsed_seconds": time.monotonic()-start})
    atomic_json(OUT, report)
    log(f"COMPLETE inclusion={condition}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
