"""Validation-only directed source initial-law support over a positive 39D box.

The candidate z is one fixed exact-dyadic vector obtained from the immutable
V3 analytic point archive; the V4 interval operator varies over the entire box.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools/pyflint"), str(ROOT / "scripts"), str(ROOT / "methods")]

import numpy as np
import yaml
from flint import arb, arb_mat, ctx, fmpq

from auxiliary_joint_qmc_gaussian_v4 import AuxiliaryJointQMCGaussianV4
import locked_author_likelihood as author
from auxiliary_v4_actual_arb_krawczyk_pilot_v2 import matrix_from_sparse
from auxiliary_v4_eaware_graph_arb_v4 import interval_ST
from auxiliary_v4_source_green_arb_v2 import block
from auxiliary_v4_source_stationary_support_arb_v1 import (
    COLS, atomic_json, eye, inf_norm_upper, qfloat, sha, transpose, unserial,
)
from independent_multispell_milp_v1 import load_capture
from independent_multispell_uncapped_farkas_v3 import HIGH_DENSITY

BOX = ROOT / "evidence/auxiliary_v4_theta39_box_support_v1.json"
SOURCE = ROOT / "evidence/auxiliary_v4_source_primitive_arb_analytic_theta39_box_exp120_v4"
GRAPH = ROOT / "evidence/auxiliary_v4_eaware_graph_arb_analytic_theta39_box_exp120_radius1over1000000000000000000000000_v4.json"
GREEN = ROOT / "evidence/auxiliary_v4_eaware_green_arb_analytic_theta39_box_exp120_radius1over1000000000000000000000000_v4"
POINT = ROOT / "evidence/auxiliary_v4_eaware_green_arb_analytic_exact_radius1over1000000000000000000000000_v3"
ROWS = ROOT / "evidence/auxiliary_v4_raw_row_elimination_control_v1.json"
INPUT = ROOT / "evidence/auxiliary_joint_qmc_2009_stress_v4/high_density_flag_e625_seed456789_N8192_flagged_inputs.npz"
OUT = ROOT / "evidence/auxiliary_v4_theta39_source_stationary_support_arb_v4.json"
LOG = ROOT / "logs/auxiliary_v4_theta39_source_stationary_support_arb_v4.log"


def log(s: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat()} {s}\n")
        f.flush()
        import os
        os.fsync(f.fileno())


def mat(entry: dict) -> arb_mat:
    n, m = entry["shape"]
    if len(entry["entries"]) != n*m:
        raise RuntimeError("Matrix archive length changed")
    result = arb_mat(n, m)
    for i in range(n):
        for j in range(m):
            result[i, j] = unserial(entry["entries"][i*m+j])
    return result


def control(F: arb_mat, Eshock: arb_mat) -> arb_mat:
    powers = [eye(25)]
    for _ in range(3):
        powers.append(powers[-1]*F)
    allcols = arb_mat(25, 32)
    for k in range(4):
        part = powers[k]*Eshock
        for i in range(25):
            for j in range(8):
                allcols[i, 8*k+j] = part[i, j]
    C = arb_mat(25, len(COLS))
    for i in range(25):
        for j, col in enumerate(COLS):
            C[i, j] = allcols[i, col]
    return C


def sigmas(box: dict) -> list[dict]:
    d, _ = load_capture(HIGH_DENSITY)
    author.META = Path("O:/RUN/jw_win_v1/source/OBC_CLOSURE_AUDIT_AND_REPAIR_PACKAGE/"
                       "workspace/obc_repro/OBC_frontier_run/reference/previous/raw/"
                       "rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz")
    target = AuxiliaryJointQMCGaussianV4(N=2, seed=456789, l_max=4, k_max=32)
    full = target.full_from_estimated(np.asarray(d["theta"], float))
    raw = target.raw.replace("^", "**").replace(";", "")
    raw = raw.replace("\n ~ ", "\n - ").replace("\n  ~ ", "\n  - ").replace("   ~ ", "   - ")
    source = yaml.load(raw, Loader=yaml.BaseLoader)
    result = []
    for shock in target.shocks:
        name = str(target.yy["calibration"]["covariances"][shock])
        index = target.pnames.index(name)
        center_impl = Fraction(*float(full[index]).as_integer_ratio())
        if name in target.prior_names:
            k = target.prior_names.index(name)
            check = box["coordinate_checks"][k]
            if check["prior_name"] != name or check["center_binary64_hex"] != float(full[index]).hex():
                raise RuntimeError("Sigma theta coordinate mismatch")
            lo = Fraction(*check["lower_ratio"])
            hi = Fraction(*check["upper_ratio"])
            if not (check["strict_archive_bounds"] and check["strict_natural_prior_support"]):
                raise RuntimeError("Sigma prior support not strict")
        else:
            exact = Fraction(Decimal(source["calibration"]["parameters"][name]))
            lo = hi = exact
        if lo <= 0 or center_impl <= 0:
            raise RuntimeError("Nonpositive source shock standard deviation")
        result.append({"shock": str(shock), "parameter": name,
                       "analytic_lower_ratio": [lo.numerator, lo.denominator],
                       "analytic_upper_ratio": [hi.numerator, hi.denominator],
                       "implementation_center_ratio": [center_impl.numerator, center_impl.denominator],
                       "strictly_positive_over_theta_box": True})
    if len(result) != 8:
        raise RuntimeError("Expected eight source shocks")
    return result


def main() -> None:
    ctx.prec, ctx.threads = 256, 1
    log("START")
    box = json.loads(BOX.read_text(encoding="utf-8"))
    if not (box["theta_dimension"] == 39 and box["radius_exponent"] == 120
            and box["all_strict_archive_prior_bounds"] and box["all_strict_natural_prior_support"]
            and box["same_binary64_implementation_parameter_on_box"]
            and all(c["same_binary64_for_entire_real_box"] for c in box["coordinate_checks"])):
        raise RuntimeError("Positive theta-box prior/rounding-cell certificate absent")
    sigma_rows = sigmas(box)
    log("SIGMA_SUPPORT_COMPLETE")

    sr = json.loads((SOURCE/"report.json").read_text(encoding="utf-8"))
    kr = json.loads(GRAPH.read_text(encoding="utf-8"))
    gr = json.loads((GREEN/"report.json").read_text(encoding="utf-8"))
    archive = GREEN/"green_interval_matrices.json.gz"
    point_archive = POINT/"green_interval_matrices.json.gz"
    point_report = json.loads((POINT/"report.json").read_text(encoding="utf-8"))
    source_archive = SOURCE/"primitive_interval_sparse.json.gz"
    if not (sr["interval_archive_sha256"] == sha(source_archive)
            and kr["krawczyk_inclusion_verified"] and kr["source_sha256"]["source_archive"] == sha(source_archive)
            and gr["green_interval_archive_sha256"] == sha(archive) and gr["power_contraction_verified"]
            and gr["source_sha256"]["source_archive"] == sha(source_archive)
            and point_report["green_interval_archive_sha256"] == sha(point_archive)):
        raise RuntimeError("Required V4 and fixed-point source chain hash mismatch")
    with gzip.open(source_archive, "rt", encoding="utf-8") as f:
        primitive = json.load(f)
    Praw, Nraw = matrix_from_sparse(primitive["Praw"]), matrix_from_sparse(primitive["Nraw"])
    # These eight original source rows are exact algebra, not a numerical
    # approximation: q_shock,t+1=0 and no forward-price coefficient.
    for i in range(8):
        for j in range(64):
            if Praw[55+i, j] != (-1 if j == 25+i else 0) or Nraw[55+i, j] != 0:
                raise RuntimeError("Original shock reset row is not exact")
    with gzip.open(archive, "rt", encoding="utf-8") as f:
        gd = json.load(f)
    with gzip.open(point_archive, "rt", encoding="utf-8") as f:
        point = json.load(f)
    omega, lam = mat(gd["omega"]), mat(gd["lambda"])
    point_lam = mat(point["lambda"])
    rr = json.loads(ROWS.read_text(encoding="utf-8"))
    top, bottom = rr["top_row_indices"], rr["bottom_row_indices"]
    if top[24:32] != list(range(55, 63)):
        raise RuntimeError("Shock rows changed in original row partition")
    S, T = interval_ST(Praw, Nraw, top, bottom)
    Eproj = arb_mat(33, 33)
    for i in range(25):
        Eproj[i, i] = 1
    Sqq, Sqp = block(S, range(33), range(33)), block(S, range(33), range(33, 64))
    Tqq = block(T, range(33), range(33))
    Aplain = Sqq+Sqp*omega
    Aprojected = Sqq+Sqp*omega*Eproj
    for i in range(8):
        for j in range(33):
            expected = -1 if j == 25+i else 0
            if Aplain[24+i, j] != expected or Aprojected[24+i, j] != expected or Tqq[24+i, j] != 0:
                raise RuntimeError("Unprojected/projected shock row bridge fails")
    for i in range(33):
        for j in range(25):
            if not (Aplain[i, j]-Aprojected[i, j]).contains(0):
                raise RuntimeError("Projected and unprojected physical block differ")
    det_plain = Aplain.det()
    det_proj = Aprojected.det()
    if det_plain.contains(0) or det_proj.contains(0):
        raise RuntimeError("Original preprocess qmat inverse not interval-certified")
    # Exact shock rows imply (I-E)*Lambda_E=0 for the true operator. Its
    # upper physical block is shared by Aplain and Aprojected, so both qmat
    # and Lambda_E are the unique solution of the same source row equations.
    bridge_residual = Aplain*lam-Tqq
    if not all(bridge_residual[i, j].contains(0) for i in range(33) for j in range(33)):
        raise RuntimeError("Interval qmat=Lambda source residual excluded zero")
    log("ORIGINAL_QMAT_BRIDGE_COMPLETE")

    with np.load(INPUT, allow_pickle=False) as z:
        qin = np.asarray(z["auxiliary_inputs"][0], float)
    if hashlib.sha256(np.ascontiguousarray(qin).tobytes()).hexdigest() != "418ef31d9f7bebe39f5787dc6bceabf24d7eb96583e871ccaddd7563ffa1a8d5":
        raise RuntimeError("Incoming input hash changed")
    def blocks(L: arb_mat):
        F, Sh = arb_mat(25, 25), arb_mat(25, 8)
        for i in range(25):
            for j in range(25):
                F[i, j] = L[i, j]
            for j in range(8):
                Sh[i, j] = L[i, 25+j]
        return F, Sh
    F, Sh = blocks(lam)
    point_F, point_Sh = blocks(point_lam)
    C = control(F, Sh)
    point_C = control(point_F, point_Sh)
    gram = transpose(C)*C
    det_gram = gram.det()
    if det_gram.contains(0) or not det_gram > 0:
        raise RuntimeError("Theta-box source controllability Gram not positive")
    point_C_float = np.asarray([[float(point_C[i, j].mid()) for j in range(len(COLS))] for i in range(25)], float)
    z_float, *_ = np.linalg.lstsq(point_C_float, qin[:25], rcond=None)
    Z = arb_mat(len(COLS), 1)
    for i, value in enumerate(z_float):
        Z[i, 0] = qfloat(float(value))
    X = C*Z
    radius = max(abs(X[i, 0]-qfloat(qin[i])).upper() for i in range(25))
    if not radius < arb(fmpq(4, 10**7)):
        raise RuntimeError("Fixed source support point outside input ball")
    contraction = inf_norm_upper(F**256)
    if not contraction < 1:
        raise RuntimeError("Theta-box physical F contraction fails")
    result = {
        "scope": "VALIDATION_ONLY_THETA39_E_AWARE_ORIGINAL_SOURCE_INITIAL_LAW_SUPPORT_V4",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": {"runner": sha(Path(__file__)), "box": sha(BOX),
                          "source_report": sha(SOURCE/"report.json"), "source_archive": sha(source_archive),
                          "graph_report": sha(GRAPH), "green_report": sha(GREEN/"report.json"),
                          "green_archive": sha(archive), "point_green_report": sha(POINT/"report.json"),
                          "point_green_archive": sha(point_archive), "row_partition": sha(ROWS),
                          "input": sha(INPUT), "support_helper": sha(ROOT/"scripts/auxiliary_v4_source_stationary_support_arb_v1.py")},
        "theta_dimension": 39, "theta_radius_exponent": 120,
        "prior_open_box_positive_measure": True,
        "binary64_implementation_constant_on_box": True,
        "source_shock_standard_deviation_intervals": sigma_rows,
        "all_eight_source_sigmas_strictly_positive": True,
        "raw_shock_reset_rows_exact": True,
        "unprojected_A_determinant_encloses_zero": bool(det_plain.contains(0)),
        "projected_A_determinant_encloses_zero": bool(det_proj.contains(0)),
        "unprojected_A_determinant_interval": str(det_plain),
        "projected_A_determinant_interval": str(det_proj),
        "original_preprocess_qmat_equals_E_aware_Lambda_by_source_algebra": True,
        "interval_bridge_residual_contains_zero": True,
        "selected_controllability_column_indices": COLS,
        "fixed_z_binary64_hex": [float(v).hex() for v in z_float],
        "selected_21d_gram_determinant_lower": str(det_gram.lower()),
        "selected_21d_gram_positive": True,
        "fixed_source_support_point_to_qin25_linf_upper": str(radius),
        "source_Fphysical_power256_inf_norm_upper": str(contraction),
        "source_stationary_gaussian_support_intersects_33d_input_ball": True,
        "complete_valid_structural_history_proven": False,
        "scientific_impossibility_proven": False,
        "endpoint_result": False,
    }
    atomic_json(OUT, result)
    log("COMPLETE")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
