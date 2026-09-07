#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ortools.sat.python import cp_model


def build_model(rows: list[dict], basin_totals: list[int], *, preserve_order: bool):
    n_columns = len(basin_totals)
    model = cp_model.CpModel()
    x: dict[tuple[int, int, int], cp_model.IntVar] = {}
    positions: dict[tuple[int, int], cp_model.IntVar] = {}

    for r, row in enumerate(rows):
        counts = [int(v) for v in row["ordered_nonzero_basin_counts"]]
        for j, _value in enumerate(counts):
            vars_j = []
            for c in range(n_columns):
                var = model.new_bool_var(f"x_r{r}_j{j}_c{c}")
                x[(r, j, c)] = var
                vars_j.append(var)
            model.add(sum(vars_j) == 1)
            pos = model.new_int_var(0, n_columns - 1, f"pos_r{r}_j{j}")
            positions[(r, j)] = pos
            model.add(pos == sum(c * x[(r, j, c)] for c in range(n_columns)))

        # A basin may contain at most one of this row's transcribed non-zero counts.
        for c in range(n_columns):
            model.add(sum(x[(r, j, c)] for j in range(len(counts))) <= 1)

        if preserve_order:
            for j in range(len(counts) - 1):
                model.add(positions[(r, j)] < positions[(r, j + 1)])
        else:
            # Equal values are exchangeable; ordering identical entries removes
            # pure label symmetry without adding biological information.
            by_value: dict[int, list[int]] = {}
            for j, value in enumerate(counts):
                by_value.setdefault(value, []).append(j)
            for ids in by_value.values():
                for a, b in zip(ids, ids[1:]):
                    model.add(positions[(r, a)] < positions[(r, b)])

    for c, total in enumerate(basin_totals):
        terms = []
        for r, row in enumerate(rows):
            counts = [int(v) for v in row["ordered_nonzero_basin_counts"]]
            for j, value in enumerate(counts):
                terms.append(value * x[(r, j, c)])
        model.add(sum(terms) == int(total))

    return model, x


def solve_layouts(
    rows: list[dict],
    basin_totals: list[int],
    *,
    preserve_order: bool,
    n_layouts: int,
    seed: int,
    per_solve_seconds: float,
):
    model, x = build_model(rows, basin_totals, preserve_order=preserve_order)
    n_columns = len(basin_totals)
    layouts: list[list[list[int]]] = []
    statuses: list[str] = []

    for solution_index in range(int(n_layouts)):
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = float(per_solve_seconds)
        solver.parameters.num_search_workers = 8
        solver.parameters.random_seed = int(seed + 104729 * solution_index)
        solver.parameters.cp_model_presolve = True
        solver.parameters.symmetry_level = 2
        status = solver.solve(model)
        status_name = solver.status_name(status)
        statuses.append(status_name)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            break

        reconstructed: list[list[int]] = []
        selected_vars = []
        for r, row in enumerate(rows):
            counts = [int(v) for v in row["ordered_nonzero_basin_counts"]]
            vector = [0] * n_columns
            for j, value in enumerate(counts):
                chosen_cols = [c for c in range(n_columns) if solver.value(x[(r, j, c)]) == 1]
                if len(chosen_cols) != 1:
                    raise RuntimeError("CP-SAT item assignment drift")
                c = chosen_cols[0]
                vector[c] = value
                selected_vars.append(x[(r, j, c)])
            if sum(vector) != int(row["N"]):
                raise RuntimeError("row total drift after CP-SAT reconstruction")
            reconstructed.append(vector)

        column_sums = [sum(row[c] for row in reconstructed) for c in range(n_columns)]
        if column_sums != [int(v) for v in basin_totals]:
            raise RuntimeError("column total drift after CP-SAT reconstruction")
        if reconstructed in layouts:
            raise RuntimeError("CP-SAT returned a duplicate reconstructed layout")
        layouts.append(reconstructed)

        # Exclude this exact assignment before solving again.
        model.add(sum(selected_vars) <= len(selected_vars) - 1)

    return layouts, statuses


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate basin layouts consistent with published Queensland Table 1 sampling margins.")
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--n-layouts", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--per-solve-seconds", type=float, default=30.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.support.read_text())
    if payload.get("schema") != "ttf_queensland33_table1_support_v0.1":
        raise RuntimeError("unexpected support schema")
    basins = list(payload["basin_column_order"])
    rows = list(payload["rows"])
    basin_totals = [int(payload["basin_totals"][b]) for b in basins]
    if len(rows) != 33 or len(basins) != 8:
        raise RuntimeError("Queensland Table 1 dimension drift")

    ordered_layouts, ordered_statuses = solve_layouts(
        rows,
        basin_totals,
        preserve_order=True,
        n_layouts=int(args.n_layouts),
        seed=int(args.seed),
        per_solve_seconds=float(args.per_solve_seconds),
    )

    if len(ordered_layouts) >= int(args.n_layouts):
        chosen_layouts = ordered_layouts[: int(args.n_layouts)]
        layout_mode = "ordered_nonzero_counts"
        fallback_reason = None
        fallback_statuses: list[str] = []
    else:
        unordered_layouts, fallback_statuses = solve_layouts(
            rows,
            basin_totals,
            preserve_order=False,
            n_layouts=int(args.n_layouts),
            seed=int(args.seed) + 1000003,
            per_solve_seconds=float(args.per_solve_seconds),
        )
        if len(unordered_layouts) < int(args.n_layouts):
            raise RuntimeError(
                "could not generate requested layouts under either ordered or conservative multiset constraints; "
                f"ordered_statuses={ordered_statuses}, multiset_statuses={fallback_statuses}, "
                f"ordered_n={len(ordered_layouts)}, multiset_n={len(unordered_layouts)}"
            )
        chosen_layouts = unordered_layouts[: int(args.n_layouts)]
        layout_mode = "unordered_nonzero_count_multiset"
        fallback_reason = (
            "ordered model did not yield the requested layout envelope within the frozen solve budget; "
            "the fallback preserves every row count multiset, every row total, and all eight published basin totals, "
            "but makes no claim that plain-text extraction preserved within-row column order"
        )

    layouts = []
    for layout_index, reconstructed in enumerate(chosen_layouts):
        layout_rows = []
        for row_id, row in enumerate(rows):
            vector = reconstructed[row_id]
            layout_rows.append({
                "code": row["code"],
                "species": row["species"],
                "eligible": bool(row["eligible"]),
                "basin_counts": {basin: int(vector[j]) for j, basin in enumerate(basins)},
                "supported_basins_ge3": [basin for j, basin in enumerate(basins) if int(vector[j]) >= 3],
            })
        layouts.append({"layout_index": layout_index, "rows": layout_rows})

    output = {
        "schema": "ttf_queensland33_feasible_layouts_v0.2",
        "source_schema": payload["schema"],
        "basin_column_order": basins,
        "basin_totals": payload["basin_totals"],
        "n_layouts": len(layouts),
        "layout_mode": layout_mode,
        "solver": {
            "engine": "ortools.cp_sat",
            "master_seed": int(args.seed),
            "per_solve_seconds": float(args.per_solve_seconds),
            "ordered_statuses": ordered_statuses,
            "fallback_statuses": fallback_statuses,
            "fallback_reason": fallback_reason,
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
    print(json.dumps({
        "n_layouts": len(layouts),
        "layout_mode": layout_mode,
        "eligible_species": payload["n_species_eligible"],
        "ordered_statuses": ordered_statuses,
        "fallback_statuses": fallback_statuses,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
