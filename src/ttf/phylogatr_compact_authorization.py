from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Mapping

from .genetic_geometry_io import sha256_path


COMPACT_EXECUTION_RULE = Path(
    "docs/supporting/genetic_phylogatr_phase3_compact_execution_v0.1.json"
)
COMPACT_EXECUTION_CODE_PATHS = (
    "src/ttf/phylogatr_compact_ibd.py",
    "src/ttf/phylogatr_compact_execution.py",
    "src/ttf/phylogatr_compact_authorization.py",
    "docs/supporting/genetic_phylogatr_phase3_compact_execution_v0.1.json",
    "scripts/authorize_phylogatr_phase3_compact_gate_d.py",
    "pyproject.toml",
)


def _load_compact_rule(repo_root: Path) -> dict:
    path = repo_root / COMPACT_EXECUTION_RULE
    if not path.is_file():
        raise FileNotFoundError(f"compact execution rule missing: {COMPACT_EXECUTION_RULE}")
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_genetic_phylogatr_phase3_compact_execution_v0.1":
        raise RuntimeError("unexpected compact Phase-3 execution rule schema")
    if payload.get("status") != (
        "FROZEN_BEFORE_ANY_SUCCESSFUL_FRESH_PHASE3_SYNTHETIC_REFERENCE_RESULT"
    ):
        raise RuntimeError("compact execution rule is not frozen at the pre-result state")
    firewall = payload.get("firewall")
    if not isinstance(firewall, dict) or not firewall:
        raise RuntimeError("compact execution rule firewall missing")
    if any(value is not False for value in firewall.values()):
        raise RuntimeError("compact execution rule firewall is open")
    return payload


def augment_compact_phase3_authorization(
    authorization: Mapping[str, object],
    *,
    repo_root: Path = Path("."),
) -> dict:
    """Add execution-only compact-code provenance to a valid Gate-D authorization.

    Scientific fields are copied unchanged.  The only additions are exact hashes
    for the compact execution surface and a receipt describing the pre-result
    amendment that justified using it.
    """
    if authorization.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1":
        raise RuntimeError("compact augmentation requires a Phase-3 Gate-D authorization")
    if authorization.get("status") != "authorize_frozen_fresh_phylogatr_phase3_gate_d":
        raise RuntimeError("compact augmentation requires an authorized Phase-3 Gate-D payload")
    if "execution_amendment" in authorization:
        raise RuntimeError("Phase-3 authorization is already execution-amended")

    root = Path(repo_root).resolve()
    rule = _load_compact_rule(root)
    result = deepcopy(dict(authorization))
    frozen = result.get("frozen_code_sha256")
    if not isinstance(frozen, dict) or not frozen:
        raise RuntimeError("base Phase-3 authorization lacks frozen code provenance")
    frozen = dict(frozen)
    for relative in COMPACT_EXECUTION_CODE_PATHS:
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"compact execution code file missing: {relative}")
        frozen[relative] = sha256_path(path)
    result["frozen_code_sha256"] = frozen
    result["execution_amendment"] = {
        "schema": str(rule["schema"]),
        "rule_sha256": sha256_path(root / COMPACT_EXECUTION_RULE),
        "scientific_result_seen_before_amendment": bool(
            rule["trigger"]["scientific_result_seen_before_amendment"]
        ),
        "execution_change_only": True,
        "legacy_gate_d_authorization_preserved": True,
        "compact_code_paths": list(COMPACT_EXECUTION_CODE_PATHS),
        "scope": dict(rule["scope"]),
    }
    return result


__all__ = [
    "COMPACT_EXECUTION_CODE_PATHS",
    "COMPACT_EXECUTION_RULE",
    "augment_compact_phase3_authorization",
]
