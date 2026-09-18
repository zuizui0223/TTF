#!/usr/bin/env python3
"""Render frozen descriptive genetic Phase-4 figures without new inference."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


SCORES = Path("manuscript/generated/genetic_ttf_phase4_v0.1/species_scores.csv")
HANDOFF = Path("benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json")
OUTPUT = Path("manuscript/figures/genetic_phase4_v0.1")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["svg.hashsalt"] = "ttf-genetic-phase4-v0.1"

    with SCORES.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    handoff = json.loads(HANDOFF.read_text())
    empirical = handoff["empirical_result"]

    primary = sorted(float(row["post_ibd_transfer"]) for row in rows)
    mean_primary = sum(primary) / len(primary)
    expected_primary = float(empirical["primary_place_beyond_ibd"]["statistic"])
    if abs(mean_primary - expected_primary) > 1e-12:
        raise RuntimeError("species table does not reproduce frozen primary statistic")

    OUTPUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.scatter(range(1, len(primary) + 1), primary, s=18)
    ax.axhline(0, linewidth=1)
    ax.axhline(mean_primary, linestyle="--", linewidth=1.2)
    ax.set_xlabel("Evaluation species, ordered by post-IBD transfer coefficient")
    ax.set_ylabel("Held-out post-IBD transfer coefficient")
    ax.set_title(f"Held-out post-IBD transfer across {len(primary)} evaluation species")
    ax.text(
        0.99,
        0.04,
        f"species-equal mean T = {mean_primary:.4f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )
    fig.tight_layout()
    species_figure = OUTPUT / "figure2_post_ibd_species_scores.svg"
    fig.savefig(species_figure, format="svg", metadata={"Date": None})
    plt.close(fig)

    cross = handoff["formal_phase3_gate_d"]
    self_gate = handoff["self_detectability"]
    labels = [
        "Cross-species Type-I upper bound",
        "Within-species self Type-I upper bound",
        "Cross-species shared-A2 power lower bound",
        "Within-species self private-A2 power lower bound",
    ]
    observed = [
        float(cross["max_private_null_wilson95_upper"]),
        float(self_gate["null_wilson95_upper"]),
        float(cross["shared_A2_wilson95_lower"]),
        float(self_gate["private_A2_wilson95_lower"]),
    ]
    thresholds = [
        float(cross["type1_ceiling"]),
        float(cross["type1_ceiling"]),
        float(cross["power_floor"]),
        float(cross["power_floor"]),
    ]

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    y = list(range(len(labels)))
    ax.scatter(observed, y, s=55, label="Observed Wilson bound")
    ax.scatter(thresholds, y, marker="x", s=65, label="Frozen threshold")
    for index, (value, threshold) in enumerate(zip(observed, thresholds)):
        relation = "≤ ceiling" if index < 2 else "≥ floor"
        ax.text(
            max(value, threshold) + 0.025,
            index,
            f"{value:.3f}  ({relation})",
            va="center",
            fontsize=9,
        )
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Probability bound")
    ax.set_title("Exact-geometry qualification margins before empirical opening")
    ax.legend(loc="lower center")
    ax.invert_yaxis()
    fig.tight_layout()
    qualification_figure = OUTPUT / "figure1_qualification_margins.svg"
    fig.savefig(qualification_figure, format="svg", metadata={"Date": None})
    plt.close(fig)

    manifest = {
        "schema": "ttf_genetic_phase4_figure_manifest_v0.1",
        "status": "DESCRIPTIVE_VISUALIZATION_ONLY",
        "evaluation_species_count": len(rows),
        "primary_species_equal_mean": mean_primary,
        "terminal_empirical_result_sha256": empirical["result_json_sha256"],
        "qualification": {
            "cross_species_type1_wilson_upper": observed[0],
            "self_type1_wilson_upper": observed[1],
            "cross_species_power_wilson_lower": observed[2],
            "self_power_wilson_lower": observed[3],
            "type1_ceiling": thresholds[0],
            "power_floor": thresholds[2],
        },
        "new_inference_performed": False,
        "post_result_retuning_performed": False,
        "input_sha256": {
            "species_scores.csv": _sha256(SCORES),
            "phase4_empirical_handoff.json": _sha256(HANDOFF),
        },
        "output_sha256": {
            qualification_figure.name: _sha256(qualification_figure),
            species_figure.name: _sha256(species_figure),
        },
    }
    (OUTPUT / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({"status": manifest["status"], "output": str(OUTPUT)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
