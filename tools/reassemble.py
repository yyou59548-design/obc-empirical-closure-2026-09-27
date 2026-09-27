#!/usr/bin/env python3
"""Verify GitHub Release parts and reconstruct the OBC full ZIP atomically."""

import argparse
import hashlib
import json
import os
from pathlib import Path


def digest_file(path: Path) -> tuple[int, str]:
    h = hashlib.sha256()
    count = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
            count += len(block)
    return count, h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parts_dir", type=Path, help="Directory containing every downloaded Release part listed in the manifest")
    parser.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parents[1] / "release" / "OBC_SPLIT_MANIFEST.json")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    output = args.output or args.parts_dir / manifest["archive_name"]
    output = output.resolve()
    temporary = output.with_name(output.name + ".assembling")
    if output.exists() or temporary.exists():
        parser.error(f"Output or temporary file already exists: {output} / {temporary}")

    parts = manifest["parts"]
    if [part["name"] for part in parts] != manifest["restore_order"]:
        parser.error("Manifest part order mismatch")
    for part in parts:
        path = args.parts_dir / part["name"]
        actual_bytes, actual_sha = digest_file(path)
        if actual_bytes != part["bytes"] or actual_sha.lower() != part["sha256"].lower():
            parser.error(f"Part failed verification: {path}")
        print(f"Verified {path.name}: {actual_bytes} bytes, SHA-256 {actual_sha}", flush=True)

    full_hash = hashlib.sha256()
    total = 0
    with temporary.open("xb") as target:
        for part in parts:
            with (args.parts_dir / part["name"]).open("rb") as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    target.write(block)
                    full_hash.update(block)
                    total += len(block)
        target.flush()
        os.fsync(target.fileno())
    if total != manifest["archive_bytes"] or full_hash.hexdigest().lower() != manifest["archive_sha256"].lower():
        parser.error(f"Reconstructed archive failed verification; inspect {temporary}")
    os.replace(temporary, output)
    print(f"Restored {output}: {total} bytes, SHA-256 {full_hash.hexdigest()}")


if __name__ == "__main__":
    main()
