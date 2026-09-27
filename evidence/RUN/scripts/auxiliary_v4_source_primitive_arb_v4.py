"""Arb interval pilot from locked source expressions to 55 primitive rows.

This produces source-analytic interval coefficients and hulls each one with
the frozen binary64 implementation. It does not verify the 64D Klein graph.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "tools" / "pyflint"), str(ROOT / "scripts"),
                 str(ROOT / "methods")]

import numpy as np
import sympy as sp
import yaml
from flint import arb, ctx, fmpq
from threadpoolctl import threadpool_limits

from auxiliary_joint_qmc_gaussian_v4 import AuxiliaryJointQMCGaussianV4
from auxiliary_v4_actual_graph_pilot_v1 import raw_pencil
import locked_author_likelihood as author
import model_engine as engine
from independent_multispell_milp_v1 import load_capture
from independent_multispell_uncapped_farkas_v3 import HIGH_DENSITY

LOG = ROOT / "logs" / "auxiliary_v4_source_primitive_arb_v4.log"
META = Path("O:/RUN/jw_win_v1/source/OBC_CLOSURE_AUDIT_AND_REPAIR_PACKAGE/"
            "workspace/obc_repro/OBC_frontier_run/reference/previous/raw/"
            "rank_spreads_exo_ztrend_SW_BAA_6420_ninit_0_meta.npz")
PREC = 256
MODE = os.environ.get("OBC_SOURCE_INTERVAL_MODE", "hull_implementation")
if MODE not in ("hull_implementation", "analytic_exact", "analytic_theta39_box_exp120"):
    raise ValueError(f"Unknown source interval mode {MODE}")
OUT = ROOT / "evidence" / f"auxiliary_v4_source_primitive_arb_{MODE}_v4"
THETA_BOX = ROOT / "evidence" / "auxiliary_v4_theta39_box_support_v1.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def log(line: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now(timezone.utc).isoformat()} {line}\n")
        f.flush()
        os.fsync(f.fileno())


def atomic_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, sort_keys=True, indent=2, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def atomic_gzip_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with gzip.open(temp, "wt", encoding="utf-8") as f:
        json.dump(obj, f, sort_keys=True)
    os.replace(temp, path)


def aq(frac: Fraction) -> arb:
    return arb(fmpq(frac.numerator, frac.denominator))


def qfloat(value: float) -> arb:
    n, d = float(value).as_integer_ratio()
    return arb(fmpq(n, d))


def eval_ast(text: str, values: dict[str, arb]) -> arb:
    tree = ast.parse(text.replace("^", "**"), mode="eval")

    def walk(node):
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant):
            token = ast.get_source_segment(text.replace("^", "**"), node)
            return aq(Fraction(Decimal(token)))
        if isinstance(node, ast.Name):
            return values[node.id]
        if isinstance(node, ast.UnaryOp):
            val = walk(node.operand)
            if isinstance(node.op, ast.USub):
                return -val
            if isinstance(node.op, ast.UAdd):
                return val
        if isinstance(node, ast.BinOp):
            a, b = walk(node.left), walk(node.right)
            if isinstance(node.op, ast.Add): return a+b
            if isinstance(node.op, ast.Sub): return a-b
            if isinstance(node.op, ast.Mult): return a*b
            if isinstance(node.op, ast.Div): return a/b
            if isinstance(node.op, ast.Pow): return a**b
        if isinstance(node, ast.Call) and len(node.args) == 1 and isinstance(node.func, ast.Name):
            v = walk(node.args[0])
            if node.func.id == "sqrt": return v.sqrt()
            if node.func.id == "exp": return v.exp()
            if node.func.id == "log": return v.log()
            if node.func.id == "abs": return abs(v)
        raise ValueError(f"Unsupported source expression {ast.dump(node)}")

    return walk(tree)


def eval_sympy(expr, values: dict[str, arb]) -> arb:
    if expr == 0: return arb(0)
    if expr.is_Integer: return arb(int(expr))
    if expr.is_Rational: return arb(fmpq(int(expr.p), int(expr.q)))
    if expr.is_Float: return aq(Fraction(Decimal(str(expr))))
    if expr.is_Symbol: return values[str(expr.name)]
    if expr.is_Add:
        out = arb(0)
        for arg in expr.args: out += eval_sympy(arg, values)
        return out
    if expr.is_Mul:
        out = arb(1)
        for arg in expr.args: out *= eval_sympy(arg, values)
        return out
    if expr.is_Pow:
        return eval_sympy(expr.base, values) ** eval_sympy(expr.exp, values)
    if expr.func == sp.Abs: return abs(eval_sympy(expr.args[0], values))
    raise ValueError(f"Unsupported symbolic node {expr}")


def symbolic_matrix(mat, values):
    return [[eval_sympy(mat[i, j], values) for j in range(mat.cols)]
            for i in range(mat.rows)]


def serial(x: arb) -> list[int]:
    m, e = x.mid().man_exp()
    r, re = x.rad().man_exp()
    return [int(m), int(e), int(r), int(re)]


def main() -> None:
    threadpool_limits(1)
    ctx.prec = PREC
    ctx.threads = 1
    log("START source primitive Arb 256-bit interval pilot")
    d, _ = load_capture(HIGH_DENSITY)
    author.META = META
    target = AuxiliaryJointQMCGaussianV4(N=2, seed=456789, l_max=4, k_max=32)
    theta = np.asarray(d["theta"], float)
    box_report = None
    if MODE == "analytic_theta39_box_exp120":
        box_report = json.loads(THETA_BOX.read_text(encoding="utf-8"))
        if (box_report["theta_dimension"] != 39 or box_report["radius_exponent"] != 120 or
                box_report["source_sha256"]["capture"] != sha(HIGH_DENSITY)):
            raise RuntimeError("Frozen theta box report invalid")
    full = target.full_from_estimated(theta)
    system = target._system(full)
    raw = target.raw.replace("^", "**").replace(";", "")
    raw = raw.replace("\n ~ ", "\n - ").replace("\n  ~ ", "\n  - ").replace("   ~ ", "   - ")
    source = yaml.load(raw, Loader=yaml.BaseLoader)
    params = source["calibration"]["parameters"]
    values = {}
    for i, name in enumerate(target.pnames):
        analytic = aq(Fraction(Decimal(params[name])))
        implemented = qfloat(full[i])
        if MODE == "analytic_theta39_box_exp120" and name in target.prior_names:
            ix = target.prior_names.index(name)
            check = box_report["coordinate_checks"][ix]
            if check["prior_name"] != name or check["center_binary64_hex"] != float(full[i]).hex():
                raise RuntimeError("Theta box prior order or value changed")
            rn, rd = check["radius_ratio"]
            values[name] = implemented + arb(0, aq(Fraction(rn, rd)))
        else:
            values[name] = (analytic.union(implemented) if MODE == "hull_implementation"
                            and name not in target.prior_names else
                            implemented if name in target.prior_names else analytic)
    for name, expression in source["calibration"]["parafunc"].items():
        values[name] = eval_ast(expression, values)
    log("PARAM_INTERVALS_COMPLETE")
    _, mats = engine.build_lambdas(target.m, target.yy)
    names = ("A", "B", "C", "PSI")
    symbolic = {}
    for name, mat in zip(names, mats[:4]):
        symbolic[name] = symbolic_matrix(mat, values)
        log(f"SYMBOLIC_{name}_INTERVALS_COMPLETE shape={mat.shape}")
    # The source adds lag identities. Reproduce that algebra on interval rows.
    nbase = len(target.m["var_ordering"])
    af = np.asarray(target.funs["AA"](system.ppar), float)
    cf = np.asarray(target.funs["CC"](system.ppar), float)
    inall = (np.abs(af) >= 1e-8).any(axis=0) & (np.abs(cf) >= 1e-8).any(axis=0)
    lag_n = int(inall.sum())
    if lag_n != 6 or nbase + lag_n != 56:
        raise RuntimeError(f"Unexpected source lag expansion base={nbase}, n={lag_n}")
    cols = np.flatnonzero(inall)
    A0 = [[arb(0) for _ in range(nbase+lag_n)] for _ in range(len(symbolic["A"])+lag_n)]
    B0 = [[arb(0) for _ in range(nbase+lag_n)] for _ in range(len(symbolic["B"])+lag_n)]
    C0 = [[arb(0) for _ in range(nbase+lag_n)] for _ in range(len(symbolic["C"])+lag_n)]
    D0 = [[arb(0) for _ in range(8)] for _ in range(len(symbolic["PSI"])+lag_n)]
    for i in range(len(symbolic["A"])):
        for j in range(nbase):
            A0[i][j] = symbolic["A"][i][j]
            B0[i][j] = symbolic["B"][i][j]
            C0[i][j] = symbolic["C"][i][j]
        for j in range(8):
            D0[i][j] = -symbolic["PSI"][i][j]
    for k, j in enumerate(cols):
        B0[-lag_n+k][-lag_n+k] = arb(1)
        B0[-lag_n+k][int(j)] = arb(-1)
        for i in range(len(symbolic["C"])):
            C0[i][-lag_n+k] = C0[i][int(j)]
            C0[i][int(j)] = arb(0)
    arrays = {"A": A0, "B": B0, "C": C0, "D": D0}
    before_union_misses = 0
    max_abs_capture_gap = 0.0
    payload = {}
    for name, rows in arrays.items():
        captured = np.asarray(d["primitive_"+name], float)
        if captured.shape != (len(rows), len(rows[0])):
            raise RuntimeError(f"{name} primitive interval shape mismatch")
        entries = []
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                implementation = qfloat(captured[i, j])
                if not value.contains(implementation):
                    before_union_misses += 1
                    max_abs_capture_gap = max(max_abs_capture_gap,
                                              abs(float(value.mid())-float(captured[i, j])))
                joined = (value.union(implementation) if MODE == "hull_implementation"
                          else value)
                row[j] = joined
                if joined != 0:
                    entries.append([i, j, *serial(joined)])
        payload[name] = {"shape": list(captured.shape), "nonzero_interval_entries": entries}
        log(f"PRIMITIVE_{name}_HULL_COMPLETE nonzero={len(entries)}")

    # Source constraint coefficients and the raw pencil; all source intervals
    # are finally hulled with the exact-dyadic float implementation.
    bb = symbolic_matrix(mats[4], values)[0]
    bbpsi = symbolic_matrix(mats[5], values)[0]
    fb = [-v for v in bb[:nbase]]
    fc = [-v for v in bb[nbase:]]
    c_arg = [v.name for v in target.m["var_ordering"]].index(str(system.const_var))
    den = fb[c_arg]
    if den.contains(0):
        raise RuntimeError("Source policy normalization denominator contains zero")
    fc = [-v/den for v in fc]
    fb = [-v/den for v in fb]
    fb += [arb(0)] * (lag_n + 8)
    fc += [arb(0)] * lag_n
    fd = [-v for v in bbpsi]
    fc = [-v for v in fc + fd]
    AAe = [r + [arb(0)]*8 for r in A0] + [[arb(0)]*64 for _ in range(8)]
    BBe = [r + [arb(0)]*8 for r in B0]
    BBe += [[arb(0)]*56 + [arb(int(i == j)) for j in range(8)] for i in range(8)]
    CCe = [C0[i] + D0[i] for i in range(55)] + [[arb(0)]*64 for _ in range(8)]
    AAe += [[arb(0)]*64]
    BBe += [fb]
    CCe += [fc]
    inq = np.flatnonzero(system.inq)
    inp = np.flatnonzero(system.inp)
    if len(inq) != 33 or len(inp) != 31:
        raise RuntimeError("Frozen graph variable masks changed")
    raw_analytic = {
        "Praw": [[-BBe[i][int(j)] for j in inq] +
                 [-AAe[i][int(j)] for j in inp] for i in range(64)],
        "Nraw": [[CCe[i][int(j)] for j in inq] +
                 [BBe[i][int(j)] for j in inp] for i in range(64)],
    }
    Pfloat, Nfloat, *_ = raw_pencil(target, system)
    raw_float = {"Praw": Pfloat, "Nraw": Nfloat}
    raw_misses = 0
    raw_max_gap = 0.0
    for name, rows in raw_analytic.items():
        captured = raw_float[name]
        entries = []
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                implementation = qfloat(captured[i, j])
                if not value.contains(implementation):
                    raw_misses += 1
                    raw_max_gap = max(raw_max_gap,
                                      abs(float(value.mid())-float(captured[i, j])))
                joined = (value.union(implementation) if MODE == "hull_implementation"
                          else value)
                row[j] = joined
                if joined != 0:
                    entries.append([i, j, *serial(joined)])
        payload[name] = {"shape": [64, 64], "nonzero_interval_entries": entries}
        log(f"RAW_{name}_HULL_COMPLETE nonzero={len(entries)}")
    # Arb comparisons are definite for the frozen 1e-8 source classification.
    tau = arb("1e-8")
    mask_ambiguous = []
    mask_mismatches = []
    for col in range(64):
        def classify(items):
            definitely_zero = all(abs(v) < tau for v in items)
            definitely_nonzero = any(abs(v) > tau for v in items)
            return (False if definitely_zero else True if definitely_nonzero else None)
        qflag = classify([CCe[i][col] for i in range(63)] + [fc[col]])
        pflag = classify([AAe[i][col] for i in range(63)] +
                         [BBe[i][col] for i in range(63)])
        if qflag is None or pflag is None:
            mask_ambiguous.append(col)
        elif qflag != bool(system.inq[col]) or (pflag and not qflag) != bool(system.inp[col]):
            mask_mismatches.append(col)
    payload["source_sha256"] = {"meta": sha(META), "capture": sha(HIGH_DENSITY),
                                "runner": sha(Path(__file__))}
    if box_report is not None:
        payload["source_sha256"]["theta_box_report"] = sha(THETA_BOX)
    payload["precision_bits"] = PREC
    xbar_analytic = values["x_bar"]
    xbar_implemented = qfloat(float(d["x_bar"]))
    payload["xbar_interval"] = serial(
        xbar_analytic.union(xbar_implemented)
        if MODE == "hull_implementation" else xbar_analytic)
    payload["semantics"] = (
        "Arb hull of exact-decimal source-analytic evaluation and exact-dyadic frozen binary64 implementation"
        if MODE == "hull_implementation" else
        "Exact-decimal YAML baseline and 39 exact rational theta boxes of radius 2^-120 max(1,abs(theta0)), Arb-evaluated original source"
        if MODE == "analytic_theta39_box_exp120" else
        "Exact-decimal YAML baseline, exact-dyadic frozen estimated theta, Arb-evaluated analytic parafunc and symbolic source equations")
    atomic_gzip_json(OUT / "primitive_interval_sparse.json.gz", payload)
    report = {
        "scope": "VALIDATION_ONLY_SOURCE_TO_PRIMITIVE_ARB_INTERVAL_PILOT",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": payload["source_sha256"],
        "python_flint_version": "0.9.0",
        "arb_precision_bits": PREC,
        "coefficient_mode": MODE,
        "primitive_shapes": {k: v["shape"] for k, v in payload.items() if k in arrays},
        "source_analytic_intervals_missing_capture_before_hull": before_union_misses,
        "max_abs_analytic_mid_to_capture_on_misses": max_abs_capture_gap,
        "analytic_xbar_contains_implementation": bool(xbar_analytic.contains(xbar_implemented)),
        "raw_pencil_analytic_intervals_missing_implementation_before_hull": raw_misses,
        "raw_pencil_max_abs_analytic_mid_to_implementation_on_misses": raw_max_gap,
        "frozen_fast0_mask_ambiguous_columns": mask_ambiguous,
        "frozen_fast0_mask_mismatch_columns": mask_mismatches,
        "interval_archive": str(OUT / "primitive_interval_sparse.json.gz"),
        "interval_archive_sha256": sha(OUT / "primitive_interval_sparse.json.gz"),
        "exact_source_literal_semantics_proven": False,
        "constraint_and_raw_pencil_interval_enclosed_given_source_semantics":
            not mask_ambiguous and not mask_mismatches,
        "krawczyk_inclusion_verified": False,
        "scientific_impossibility_proven": False,
    }
    atomic_json(OUT / "report.json", report)
    log(f"COMPLETE misses={before_union_misses} gap={max_abs_capture_gap:.3e}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
