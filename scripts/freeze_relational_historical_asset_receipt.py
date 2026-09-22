#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--asset", type=Path, action="append", required=True)
    ap.add_argument("--resolved-url", action="append", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_historical_climate_exposure_rule_v0.1":
        raise RuntimeError("unexpected Study-C historical rule")
    required = list(rule["historical_climate"]["asset_contract"]["required_basenames"])
    if len(args.asset) != len(required) or len(args.resolved_url) != len(required):
        raise RuntimeError("exactly eight Study-C historical assets and URLs are required")

    by_lower = {Path(path).name.casefold(): path for path in args.asset}
    assets = []
    for basename in required:
        key = basename.casefold()
        if key not in by_lower:
            raise RuntimeError(f"missing frozen logical historical asset: {basename}")
        path = by_lower[key]
        idx = next(i for i, p in enumerate(args.asset) if p == path)
        url = str(args.resolved_url[idx])
        if Path(url.split("?", 1)[0]).name.casefold() != key:
            raise RuntimeError(f"resolved URL basename does not match frozen logical asset: {basename}")
        assets.append({
            "logical_basename": basename,
            "resolved_basename": Path(path).name,
            "resolved_url": url,
            "sha256": sha256_path(path),
            "size_bytes": path.stat().st_size,
        })

    payload = {
        "schema": "ttf_relational_historical_climate_asset_receipt_v0.1",
        "status": "FROZEN_RESPONSE_BLIND_HISTORICAL_ASSETS",
        "dataset": rule["historical_climate"]["dataset"],
        "version": rule["historical_climate"]["version"],
        "assets": assets,
        "substitution_used": False,
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": payload["status"], "assets": len(assets)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
