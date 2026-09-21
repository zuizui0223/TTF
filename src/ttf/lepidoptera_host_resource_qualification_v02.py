from __future__ import annotations

import hashlib
import numpy as np

from .lepidoptera_host_resource_qualification import nuisance_envelope_pvalues


REFERENCE_TAG="lepidoptera-host-resource-v02-reference"
EVALUATION_TAG="lepidoptera-host-resource-v02-evaluation"
ALLOWED_CELLS=("private","geometry_confounded_trap","host_breadth_confounded_trap","host_resource_gradient_positive")


def frozen_v02_seed(master_seed:int,tag:str,cell:str,replicate:int)->int:
    if tag not in {REFERENCE_TAG,EVALUATION_TAG}:
        raise ValueError("unauthorized v0.2 host-resource seed namespace")
    if cell not in ALLOWED_CELLS:
        raise ValueError("unknown v0.2 host-resource cell")
    if replicate<0:
        raise ValueError("replicate must be non-negative")
    payload=f"{int(master_seed)}|{tag}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")


__all__=[
    "ALLOWED_CELLS",
    "EVALUATION_TAG",
    "REFERENCE_TAG",
    "frozen_v02_seed",
    "nuisance_envelope_pvalues",
]
