#!/usr/bin/env python3
"""Freeze the pre-result one-shot empirical execution receipt.

This command is response-blind. It refuses to create a receipt unless the
prospectively frozen authoritative survivor qualification passed and the
hash-bound identity-opening authorization is still in a fully closed state.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.codistributed_execution import build_execution_receipt, sha256_path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--qualification",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--empirical-result",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    qualification=json.loads(args.qualification.read_text())
    authorization=json.loads(args.authorization.read_text())
    payload=build_execution_receipt(
        contract,
        qualification,
        authorization,
        contract_sha256=sha256_path(args.contract),
        qualification_sha256=sha256_path(args.qualification),
        authorization_sha256=sha256_path(args.authorization),
        empirical_result_exists=(
            args.empirical_result is not None and args.empirical_result.exists()
        ),
    )
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "survivor_run_id":payload["survivor_qualification_workflow_run_id"],
        "qualification_sha256":payload["survivor_qualification_sha256"],
        "authorization_sha256":payload["identity_opening_authorization_sha256"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
