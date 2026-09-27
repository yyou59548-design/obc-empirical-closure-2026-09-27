"""Read-only independent hashes for the V4 positive-volume chain."""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("O:/RUN")
E = ROOT / "evidence"
S = ROOT / "scripts"
HERE = Path(__file__).resolve().parent
M = E / "auxiliary_v4_theta39_eaware_chain_manifest_v4.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda: f.read(1 << 20), b""):
            h.update(buf)
    return h.hexdigest()


def same(path: Path, digest: str) -> bool:
    got = sha(path)
    if got != digest:
        raise RuntimeError(f"SHA drift: {path}: {got} != {digest}")
    return True


def main() -> None:
    d = json.loads(M.read_text(encoding="utf-8"))
    prefix = "auxiliary_v4_eaware_"
    b = "analytic_theta39_box_exp120"
    g = E / f"{prefix}graph_arb_{b}_radius1over{10**24}_v4.json"
    green = E / f"{prefix}green_arb_{b}_radius1over{10**24}_v4"
    source = E / f"auxiliary_v4_source_primitive_arb_{b}_v4"
    paths = {
        "primitive_v4": S / "auxiliary_v4_source_primitive_arb_v4.py",
        "graph_v4": S / "auxiliary_v4_eaware_graph_arb_v4.py",
        "green_v4": S / "auxiliary_v4_eaware_green_arb_v4.py",
        "dual_v4": S / "auxiliary_v4_eaware_dual_adjoint_arb_v4.py",
        "support_v4": S / "auxiliary_v4_theta39_source_stationary_support_arb_v4.py",
    }
    script_hashes = {}
    for key, path in paths.items():
        same(path, d["source_sha256"][key])
        script_hashes[key] = sha(path)
    paths2 = {
        "source_archive_sha256": source / "primitive_interval_sparse.json.gz",
        "source_report_sha256": source / "report.json",
        "graph_report_sha256": g,
        "graph_center_sha256": g.with_suffix(".center.json.gz"),
        "green_report_sha256": green / "report.json",
        "green_archive_sha256": green / "green_interval_matrices.json.gz",
        "dual_report_sha256": E / f"{prefix}dual_adjoint_arb_{b}_v4.json",
        "support_report_sha256": E / "auxiliary_v4_theta39_source_stationary_support_arb_v4.json",
        "theta_box_sha256": E / "auxiliary_v4_theta39_box_support_v1.json",
    }
    archive_hashes = {}
    for key, path in paths2.items():
        same(path, d["analytic_theta_box"][key])
        archive_hashes[key] = sha(path)
    # The V3 frozen branch is supported separately; verify the erratum's
    # existence and the exact referenced dual/support artifacts.
    erratum = E / "auxiliary_v4_eaware_v3_binary_provenance_erratum_v1.md"
    if not erratum.is_file():
        raise RuntimeError("Frozen-binary64 provenance erratum missing")
    same(E / "auxiliary_v4_eaware_dual_adjoint_arb_frozen_binary64_v3.json",
         d["implementation_binary64"]["v3_dual_report_sha256"])
    same(E / "auxiliary_v4_source_stationary_support_arb_v1/summary.json",
         d["implementation_binary64"]["v3_source_support_sha256"])
    out = {"scope": "READ_ONLY_INDEPENDENT_V4_HASH_AUDIT",
           "generated_utc": datetime.now(timezone.utc).isoformat(),
           "v4_manifest_sha256": sha(M),
           "all_direct_hash_checks_pass": True,
           "v4_script_sha256": script_hashes,
           "v4_archive_sha256": archive_hashes,
           "v3_binary_provenance_erratum_sha256": sha(erratum),
           "scientific_impossibility_proven": False}
    dst = HERE / "manifest_independent_summary.json"
    tmp = dst.with_name(dst.name + f".tmp.{os.getpid()}")
    tmp.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    tmp.replace(dst)
    print(json.dumps(out, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
