#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from ttf.phylogatr_compact_authorization import augment_compact_phase3_authorization


BASE_AUTHORIZER = Path("scripts/authorize_phylogatr_phase3_gate_d.py")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Run the frozen fresh phylogatR Gate-D authorizer and then bind the exact "
            "compact execution code hashes without changing scientific fields."
        )
    )
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase1-manifest", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    base_authorizer = repo_root / BASE_AUTHORIZER
    if not base_authorizer.is_file():
        raise FileNotFoundError(f"base Gate-D authorizer missing: {BASE_AUTHORIZER}")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(
            "compact Phase-3 authorization output already exists; use a fresh path"
        )

    with tempfile.TemporaryDirectory(prefix="ttf_phylogatr_compact_auth_") as temporary:
        base_output = Path(temporary) / "base_authorization.json"
        command = [
            sys.executable,
            str(base_authorizer),
            "--geometry",
            str(args.geometry.resolve()),
            "--phase1-manifest",
            str(args.phase1_manifest.resolve()),
            "--phase2-manifest",
            str(args.phase2_manifest.resolve()),
            "--phase3-rule",
            str(args.phase3_rule.resolve()),
            "--repo-root",
            str(repo_root),
            "--output",
            str(base_output),
        ]
        completed = subprocess.run(
            command,
            cwd=repo_root,
            check=False,
            text=True,
            capture_output=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "base fresh Phase-3 Gate-D authorization failed before compact provenance binding.\n"
                + completed.stderr
            )
        base = json.loads(base_output.read_text())

    augmented = augment_compact_phase3_authorization(base, repo_root=repo_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(augmented, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": augmented["status"],
                "geometry_fingerprint_sha256": augmented["geometry_fingerprint_sha256"],
                "master_seed": augmented["master_seed"],
                "compact_execution_authorized": True,
                "scientific_result_seen_before_amendment": False,
                "fresh_nucleotide_identity_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
