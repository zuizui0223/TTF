from __future__ import annotations

import hashlib
import numpy as np

HISTORICAL_VARIABLES = ("bio01", "bio07", "bio12", "bio15")
LGM_TIME_INDEX = -190
PRESENT_TIME_INDEX = 20


def logical_asset_basename(variable: str, time_index: int) -> str:
    variable = str(variable)
    if variable not in HISTORICAL_VARIABLES:
        raise ValueError(f"unsupported historical variable: {variable}")
    if int(time_index) not in {LGM_TIME_INDEX, PRESENT_TIME_INDEX}:
        raise ValueError(f"unsupported historical time index: {time_index}")
    return f"CHELSA_TraCE21K_{variable}_{int(time_index)}_V1.0.tif"


def historical_delta(past: np.ndarray, present: np.ndarray) -> np.ndarray:
    a = np.asarray(past, dtype=float)
    b = np.asarray(present, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 4:
        raise ValueError("past and present historical climate must be equal n x 4 matrices")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("historical climate matrices must be finite")
    return a - b


def hash_key(tag: str, source_sha256: str, species: str) -> tuple[str, str]:
    return (
        hashlib.sha256(f"{tag}|{source_sha256}|{species}".encode("utf-8")).hexdigest(),
        str(species),
    )


def panel_roles(names: list[str], source_sha256: str) -> tuple[list[str], list[str]]:
    ordered = sorted(
        map(str, names),
        key=lambda name: hash_key("relational-history-role-v0.1", source_sha256, name),
    )
    n_source = len(ordered) // 2
    return sorted(ordered[:n_source]), sorted(ordered[n_source:])


def panel_order(names: list[str], source_sha256: str) -> list[str]:
    return sorted(
        map(str, names),
        key=lambda name: hash_key("relational-history-panel-v0.1", source_sha256, name),
    )
