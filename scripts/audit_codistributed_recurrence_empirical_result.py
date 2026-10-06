#!/usr/bin/env python3
"""Audit the frozen one-shot codistributed-recurrence empirical summary.

This script is intentionally interpretation-only. It never reads sequence
identity, edge genetic distances, or dyad-level T_st.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from ttf.codistributed_post_result import audit_empirical_result


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--result",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--interpretation",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    result=json.loads(args.result.read_text())
    interpretation=json.loads(args.interpretation.read_text())
    if interpretation.get("schema")!="ttf_genetic_codistributed_recurrence_interpretation_v0.1":
        raise RuntimeError("unexpected interpretation contract")
    if interpretation.get("status")!="FROZEN_BEFORE_DEVELOPMENT_SCREEN_RESULT":
        raise RuntimeError("interpretation contract was not prospectively frozen")

    audit=audit_empirical_result(
        result,
        expected_authorization_sha256=sha256_path(args.authorization),
    )
    branch=(
        interpretation["result_branches"]["confirmatory_positive"]
        if audit.positive
        else interpretation["result_branches"]["confirmatory_qualified_null"]
    )
    payload={
        "schema":"ttf_genetic_codistributed_recurrence_post_result_audit_v0.1",
        "status":"PASS_FROZEN_ONE_SHOT_RESULT_AUDIT",
        "empirical_result_sha256":sha256_path(args.result),
        "authorization_sha256":sha256_path(args.authorization),
        "interpretation_contract_sha256":sha256_path(args.interpretation),
        "audit":asdict(audit),
        "frozen_headline":branch["headline"],
        "frozen_biological_reading":branch["biological_reading"],
        "allowed_claims":branch["allowed_claims"],
        "prohibited_claims":branch["prohibited_claims"],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "decision":audit.decision,
        "beta_G":audit.beta_G,
        "primary_p":audit.primary_p_value,
        "headline":payload["frozen_headline"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
