#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

AMPLITUDES = (2.0, 2.5, 3.0, 4.0)
NEW_AMPLITUDES = (2.5, 3.0, 4.0)
N_LAYOUTS = 12


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate the frozen-layout Queensland detection-floor scan.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text())
    if baseline.get("schema") != "ttf_queensland33_topology_qualification_frozen_v0.2":
        raise RuntimeError("unexpected baseline schema")

    cells: dict[int, dict[tuple[float, float], dict]] = {i: {} for i in range(N_LAYOUTS)}
    designs: dict[int, dict] = {}
    for path in sorted(args.input_dir.rglob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_queensland33_detection_floor_cell_v0.1":
            continue
        layout = int(payload["layout_index"])
        cell = payload["cell"]
        key = (float(cell["shared_fraction"]), float(cell["amplitude"]))
        if key in cells[layout]:
            raise RuntimeError(f"duplicate scan cell: layout={layout}, key={key}")
        cells[layout][key] = cell
        designs.setdefault(layout, payload["design"])

    expected = {(0.0, 4.0), (1.0, 2.5), (1.0, 3.0), (1.0, 4.0)}
    for layout in range(N_LAYOUTS):
        if set(cells[layout]) != expected:
            raise RuntimeError(f"layout {layout} scan grid drift: {sorted(cells[layout])}")

    per_layout = {}
    worst_power = {amp: 1.0 for amp in AMPLITUDES}
    worst_type1_A4 = 0.0
    for layout in range(N_LAYOUTS):
        baseline_layout = baseline["per_layout"][str(layout)]
        power = {2.0: float(baseline_layout["full_shared_A2_power"])}
        for amp in NEW_AMPLITUDES:
            power[amp] = float(cells[layout][(1.0, amp)]["rejection_rate"])
        private_A4 = float(cells[layout][(0.0, 4.0)]["rejection_rate"])
        worst_type1_A4 = max(worst_type1_A4, private_A4)
        for amp in AMPLITUDES:
            worst_power[amp] = min(worst_power[amp], power[amp])
        per_layout[str(layout)] = {
            "private_A4_rejection": private_A4,
            "full_shared_power": {str(amp): power[amp] for amp in AMPLITUDES},
        }

    qualified = [amp for amp in AMPLITUDES if worst_power[amp] >= 0.80]
    detection_floor = min(qualified) if qualified else None
    payload = {
        "schema": "ttf_queensland33_detection_floor_v0.1",
        "baseline": {
            "path": str(args.baseline),
            "baseline_schema": baseline["schema"],
            "baseline_A2_worst_power": float(baseline["qualification"]["worst_layout_full_shared_A2_power"]),
            "baseline_worst_type1_up_to_A3": float(baseline["qualification"]["worst_layout_max_zero_shared_rejection"]),
        },
        "design": {
            "exact_frozen_layout_envelope": True,
            "n_layouts": N_LAYOUTS,
            "amplitudes": list(AMPLITUDES),
            "new_worlds_per_cell_per_layout": 50,
            "bootstrap_resamples": 999,
            "type1_control_amplitude": 4.0,
            "power_floor": 0.80,
            "type1_ceiling": 0.10,
            "estimator_changed_after_baseline_failure": False,
            "split_changed_after_baseline_failure": False,
            "genetic_outcomes_used": False,
            "named_Mary_Brisbane_boundary_used": False,
        },
        "worst_layout_curve": {str(amp): worst_power[amp] for amp in AMPLITUDES},
        "worst_layout_private_A4_rejection": worst_type1_A4,
        "detection_floor_amplitude": detection_floor,
        "qualification": {
            "type1_A4_pass": bool(worst_type1_A4 <= 0.10),
            "has_detectable_amplitude_in_scan": detection_floor is not None,
            "passed_for_empirical_opening": False,
            "reason": "This scan estimates an applicability/detection floor only. Empirical DNA remains sealed until the chosen floor is validated with an independent world/bootstrap seed and higher replicate count."
        },
        "per_layout": per_layout,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "worst_layout_curve": payload["worst_layout_curve"],
        "private_A4": worst_type1_A4,
        "detection_floor_amplitude": detection_floor,
        "type1_A4_pass": payload["qualification"]["type1_A4_pass"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
