#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_optional(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _precision_pass(payload: dict, *, lane: str) -> bool:
    precision = payload.get("precision_qualification")
    if not isinstance(precision, dict) or "passed" not in precision:
        raise RuntimeError(f"{lane} result lacks precision qualification")
    return bool(precision["passed"])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Jointly reduce TTF Gate-I-A and Gate-I-B high-precision results."
    )
    parser.add_argument("--gate-i-a", type=Path, required=True)
    parser.add_argument("--gate-i-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    a = _load_optional(args.gate_i_a)
    b = _load_optional(args.gate_i_b)
    available = {
        "I-A": a is not None,
        "I-B": b is not None,
    }

    if a is not None:
        if a.get("schema") != "ttf_gate_i_geometry_precision_v0.1":
            raise RuntimeError("unexpected Gate-I-A schema")
        if a.get("gate") != "I-A external non-flower geometry stress":
            raise RuntimeError("Gate-I-A identity drift")
        if a.get("gate_i_b_required") is not True:
            raise RuntimeError("Gate-I-A result no longer requires Gate-I-B")
    if b is not None:
        if b.get("schema") != "ttf_gate_i_geometry_batch_precision_v0.1":
            raise RuntimeError("unexpected Gate-I-B schema")
        if b.get("gate") != "I-B intended empirical TTF sampling geometry":
            raise RuntimeError("Gate-I-B identity drift")
        source = b.get("source") or {}
        if source.get("source_commit") != "e6ad7aa76261906b6f757ffc324ea6cd062a0917":
            raise RuntimeError("Gate-I-B did not use the frozen pre-inference RGFCA source")
        if source.get("source_g1_g3_inference_run") is not False:
            raise RuntimeError("Gate-I-B source is not pre-inference")
        if source.get("external_overlay_opened") is True or source.get("source_external_overlay_opened") is True:
            raise RuntimeError("Gate-I-B source had an opened external overlay")
        if source.get("empirical_colour_columns_read") is not False:
            raise RuntimeError("Gate-I-B source read empirical colour columns")
        design = b.get("design") or {}
        if design.get("execution") != "exact_dense_projection_batched_worlds":
            raise RuntimeError("Gate-I-B execution changed the exact dense projection")
        if design.get("batching_changes_estimand") is not False:
            raise RuntimeError("Gate-I-B does not certify batching as estimand-preserving")

    complete = bool(a is not None and b is not None)
    a_pass = None if a is None else _precision_pass(a, lane="I-A")
    b_pass = None if b is None else _precision_pass(b, lane="I-B")
    overall = None if not complete else bool(a_pass and b_pass)
    status = (
        "pending"
        if not complete
        else ("passed" if overall else "failed")
    )

    payload = {
        "schema": "ttf_gate_i_joint_qualification_v0.1",
        "gate": "I actual empirical sampling geometry qualification",
        "status": status,
        "available_lanes": available,
        "gate_i_a_pass": a_pass,
        "gate_i_b_pass": b_pass,
        "overall_gate_i_pass": overall,
        "claim_ready_for_gate_i_geometry_layer": bool(overall) if complete else False,
        "downstream_gates_still_required": [
            "Gate J predictor attribution and IBD residualization where relevant",
            "Gate K observation-bias stress tests",
            "cross-domain empirical detection floor",
        ],
        "claim_boundary": (
            "A joint Gate-I pass qualifies the fixed-geometry sharedness inference layer "
            "against both an independent non-flower stress geometry and the intended "
            "pre-inference RGFCA sampling geometry. It does not establish an empirical "
            "flower-colour shared-boundary result, causal barrier mechanism, predictor "
            "attribution, or observation-bias robustness."
        ),
        "lanes": {
            "I-A": None if a is None else {
                "path": str(args.gate_i_a),
                "sha256": sha256_file(args.gate_i_a),
                "geometry_fingerprint": a["geometry"]["fingerprint_sha256"],
                "precision_qualification": a["precision_qualification"],
                "lineage": a.get("lineage"),
            },
            "I-B": None if b is None else {
                "path": str(args.gate_i_b),
                "sha256": sha256_file(args.gate_i_b),
                "geometry_fingerprint": b["geometry"]["fingerprint_sha256"],
                "precision_qualification": b["precision_qualification"],
                "source": b.get("source"),
                "lineage": b.get("lineage"),
            },
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": status,
        "gate_i_a_pass": a_pass,
        "gate_i_b_pass": b_pass,
        "overall_gate_i_pass": overall,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
