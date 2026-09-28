from __future__ import annotations

import numpy as np


def climate_validity_mask(
    columns: np.ndarray,
    *,
    nodata_values: list[float | None],
) -> np.ndarray:
    """Return rows usable across every supplied climate raster.

    This is intentionally response blind. It mirrors the finite/nodata checks
    in the frozen relation builders, but applies them before historical
    displacement arithmetic so invalid raster cells never reach subtraction.
    """
    x = np.asarray(columns, dtype=float)
    if x.ndim != 2 or x.shape[1] != len(nodata_values):
        raise ValueError("columns and nodata_values must align")
    valid = np.isfinite(x).all(axis=1)
    for j, nodata in enumerate(nodata_values):
        if nodata is not None:
            valid &= ~np.isclose(
                x[:, j],
                float(nodata),
                rtol=0.0,
                atol=0.0,
            )
    return valid
