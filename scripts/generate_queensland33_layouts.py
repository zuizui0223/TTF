#!/usr/bin/env python3
from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


def placement_options(counts: list[int], n_columns: int) -> list[tuple[int, ...]]:
    return list(combinations(range(n_columns), len(counts)))


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate multiple basin layouts consistent with the published Queensland Table 1 margins.")
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--n-layouts", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.support.read_text())
    if payload.get("schema") != "ttf_queensland33_table1_support_v0.1":
        raise RuntimeError("unexpected support schema")
    basins = list(payload["basin_column_order"])
    rows = list(payload["rows"])
    basin_totals = np.asarray([payload["basin_totals"][b] for b in basins], dtype=float)
    if len(rows) != 33 or len(basins) != 8:
        raise RuntimeError("Queensland Table 1 dimension drift")

    variable_meta: list[tuple[int, tuple[int, ...]]] = []
    row_var_ids: list[list[int]] = []
    contributions: list[np.ndarray] = []
    for row_id, row in enumerate(rows):
        counts = [int(v) for v in row["ordered_nonzero_basin_counts"]]
        ids = []
        for placement in placement_options(counts, len(basins)):
            vector = np.zeros(len(basins), dtype=float)
            for value, col in zip(counts, placement):
                vector[col] = value
            ids.append(len(variable_meta))
            variable_meta.append((row_id, placement))
            contributions.append(vector)
        row_var_ids.append(ids)

    n_vars = len(variable_meta)
    base_A = np.zeros((len(rows) + len(basins), n_vars), dtype=float)
    base_b = np.concatenate([np.ones(len(rows), dtype=float), basin_totals])
    for row_id, ids in enumerate(row_var_ids):
        base_A[row_id, ids] = 1.0
    for var_id, vector in enumerate(contributions):
        base_A[len(rows):, var_id] = vector

    rng = np.random.default_rng(int(args.seed))
    selected_solutions: list[np.ndarray] = []
    layouts: list[dict] = []
    attempts = 0
    max_attempts = max(100, int(args.n_layouts) * 30)
    while len(layouts) < int(args.n_layouts) and attempts < max_attempts:
        attempts += 1
        A_parts = [base_A]
        lo_parts = [base_b]
        hi_parts = [base_b]
        # No-good cuts guarantee every retained layout differs in at least one row placement.
        for previous in selected_solutions:
            row = np.zeros(n_vars, dtype=float)
            row[previous] = 1.0
            A_parts.append(row[None, :])
            lo_parts.append(np.asarray([-np.inf]))
            hi_parts.append(np.asarray([len(rows) - 1.0]))
        A = np.vstack(A_parts)
        lower = np.concatenate(lo_parts)
        upper = np.concatenate(hi_parts)
        objective = rng.normal(size=n_vars) + 1e-9 * np.arange(n_vars, dtype=float)
        result = milp(
            c=objective,
            integrality=np.ones(n_vars, dtype=int),
            bounds=Bounds(np.zeros(n_vars), np.ones(n_vars)),
            constraints=LinearConstraint(A, lower, upper),
            options={"presolve": True, "time_limit": 120.0},
        )
        if not result.success or result.x is None:
            break
        chosen = np.flatnonzero(result.x > 0.5)
        if len(chosen) != len(rows):
            raise RuntimeError("MILP solution did not select exactly one placement per species")
        selected_solutions.append(chosen.copy())

        reconstructed = np.zeros((len(rows), len(basins)), dtype=int)
        for var_id in chosen:
            row_id, placement = variable_meta[int(var_id)]
            counts = [int(v) for v in rows[row_id]["ordered_nonzero_basin_counts"]]
            for value, col in zip(counts, placement):
                reconstructed[row_id, col] = value
        if not np.array_equal(reconstructed.sum(axis=0), basin_totals.astype(int)):
            raise RuntimeError("column totals drifted after MILP reconstruction")
        if not np.array_equal(reconstructed.sum(axis=1), np.asarray([row["N"] for row in rows], dtype=int)):
            raise RuntimeError("row totals drifted after MILP reconstruction")

        layout_rows = []
        for row_id, row in enumerate(rows):
            layout_rows.append({
                "code": row["code"],
                "species": row["species"],
                "eligible": bool(row["eligible"]),
                "basin_counts": {basin: int(reconstructed[row_id, j]) for j, basin in enumerate(basins)},
                "supported_basins_ge3": [basin for j, basin in enumerate(basins) if reconstructed[row_id, j] >= 3],
            })
        layouts.append({"layout_index": len(layouts), "rows": layout_rows})

    if len(layouts) != int(args.n_layouts):
        raise RuntimeError(f"requested {args.n_layouts} unique layouts but found {len(layouts)}")

    output = {
        "schema": "ttf_queensland33_feasible_layouts_v0.1",
        "source_schema": payload["schema"],
        "basin_column_order": basins,
        "basin_totals": payload["basin_totals"],
        "n_layouts": len(layouts),
        "solver": {
            "engine": "scipy.optimize.milp",
            "selection": "random linear objective plus no-good cuts",
            "master_seed": int(args.seed),
            "attempts": attempts,
            "uses_genetic_outcomes": False,
            "uses_named_boundary": False,
        },
        "eligibility": {
            "eligible_codes": payload["eligible_codes"],
            "n_species": payload["n_species_eligible"],
            "rule": payload["support_rule"]["eligibility"],
        },
        "layouts": layouts,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"n_layouts": len(layouts), "n_variables": n_vars, "attempts": attempts, "eligible_species": payload["n_species_eligible"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
