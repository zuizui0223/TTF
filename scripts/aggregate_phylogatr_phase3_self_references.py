#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.self_rule.read_text())
    auth = json.loads(args.phase3_authorization.read_text())
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_self_detectability_rule_v0.1":
        raise RuntimeError("fresh self rule schema drift")
    if auth.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1":
        raise RuntimeError("fresh Phase-3 authorization schema drift")
    expected = int(rule["synthetic_worlds"]["independent_null_reference_worlds"])
    rows: dict[int, float] = {}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_genetic_phylogatr_phase3_self_reference_shard_v0.1":
            continue
        if payload.get("geometry_fingerprint_sha256") != auth["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fresh self-reference geometry drift in {path}")
        if payload.get("confirmatory_sequence_identity_opened") is not False:
            raise RuntimeError(f"fresh identity firewall drift in {path}")
        for row in payload["rows"]:
            replicate = int(row["replicate"])
            if replicate in rows:
                raise RuntimeError(f"duplicate fresh self-reference replicate {replicate}")
            rows[replicate] = float(row["statistic"])
    if sorted(rows) != list(range(expected)):
        missing = sorted(set(range(expected)) - set(rows))[:10]
        raise RuntimeError(f"incomplete fresh self references; first missing={missing}")

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_self_references_v0.1",
        "status": "complete_independent_null_reference",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "n_worlds": expected,
        "statistics": [rows[index] for index in range(expected)],
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "statistics"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
