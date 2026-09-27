#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STRATA = ["1_family", "2_families", "3_to_5_families", "6plus_families"]
STRATA_LABELS = {
    "1_family": "1 family",
    "2_families": "2 families",
    "3_to_5_families": "3–5 families",
    "6plus_families": "6+ families",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def deterministic_jitter(label: str, width: float = 0.12) -> float:
    import hashlib
    value = int(hashlib.sha256(label.encode("utf-8")).hexdigest()[:8], 16)
    return ((value / 0xFFFFFFFF) - 0.5) * 2 * width


def spearman(a, b) -> float:
    def ranks(x):
        x = np.asarray(x, float)
        order = np.argsort(x, kind="mergesort")
        out = np.empty(len(x), float)
        start = 0
        while start < len(x):
            stop = start + 1
            while stop < len(x) and x[order[stop]] == x[order[start]]:
                stop += 1
            out[order[start:stop]] = 0.5 * ((start + 1) + stop)
            start = stop
        return out

    x, y = ranks(a), ranks(b)
    return float(np.corrcoef(x, y)[0, 1])


def figure1(descriptors: list[dict[str, str]], outdir: Path) -> None:
    rows = [
        row
        for row in descriptors
        if float(row["host_family_count"]) > 0
        and int(row["host_wgsrpd3_unit_count"]) > 0
    ]
    adequate = [
        row
        for row in rows
        if int(row["resolved_host_species"]) >= float(row["host_family_count"])
    ]
    x = np.asarray([float(r["host_family_count"]) for r in adequate])
    y = np.asarray([int(r["host_wgsrpd3_unit_count"]) for r in adequate])
    rho = spearman(x, y)

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for row in adequate:
        xv = float(row["host_family_count"])
        yv = int(row["host_wgsrpd3_unit_count"])
        ax.scatter(
            xv + deterministic_jitter(row["species"]),
            yv,
            s=20,
            alpha=0.55,
        )

    q75 = float(np.quantile(y, 0.75))
    highlights = sorted(
        [
            row
            for row in adequate
            if float(row["host_family_count"]) == 1
            and int(row["host_wgsrpd3_unit_count"]) >= q75
        ],
        key=lambda r: -int(r["host_wgsrpd3_unit_count"]),
    )[:6]
    for row in highlights:
        xv = 1 + deterministic_jitter(row["species"])
        yv = int(row["host_wgsrpd3_unit_count"])
        ax.scatter(xv, yv, s=45, marker="o")
        ax.annotate(
            row["species"],
            (xv, yv),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_yscale("log")
    ax.set_xlabel("Larval host-plant families")
    ax.set_ylabel("Native host-resource breadth (WGSRPD3 units)")
    ax.set_title(f"Taxonomic and geographic specialization are only partly coupled\nSpearman ρ = {rho:.3f}, n = {len(adequate)}")
    save(fig, outdir / "Figure1_taxonomic_vs_geographic_specialization")


def boxplot_by_stratum(ax, rows, field, ylabel, zero_line=False):
    values = [
        [float(r[field]) for r in rows if r["host_breadth_stratum"] == stratum]
        for stratum in STRATA
    ]
    ax.boxplot(values, tick_labels=[STRATA_LABELS[s] for s in STRATA], showfliers=False)
    for i, stratum in enumerate(STRATA, start=1):
        subset = [r for r in rows if r["host_breadth_stratum"] == stratum]
        for row in subset:
            ax.scatter(
                i + deterministic_jitter(row["species"], 0.18),
                float(row[field]),
                s=14,
                alpha=0.45,
            )
    if zero_line:
        ax.axhline(0, linewidth=1, linestyle="--")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=20)


def figure2(anth: list[dict[str, str]], outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))

    native = np.asarray([int(r["native_resource_units"]) for r in anth], float)
    contemporary = np.asarray([int(r["contemporary_resource_units"]) for r in anth], float)
    axes[0].scatter(native, contemporary, s=18, alpha=0.5)
    lo = max(1, min(native.min(), contemporary.min()))
    hi = max(native.max(), contemporary.max())
    axes[0].plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Native resource breadth")
    axes[0].set_ylabel("Contemporary resource breadth")
    expanded = sum(int(r["introduced_added_units"]) > 0 for r in anth)
    axes[0].set_title(f"Introduced host ranges expand opportunity\n{expanded}/{len(anth)} species expanded")

    boxplot_by_stratum(
        axes[1],
        anth,
        "log_resource_expansion",
        "Log proportional resource expansion",
        zero_line=True,
    )
    rho = spearman(
        [float(r["host_family_count"]) for r in anth],
        [float(r["log_resource_expansion"]) for r in anth],
    )
    axes[1].set_title(f"No proportional generalist advantage\nSpearman ρ = {rho:.3f}")

    fig.tight_layout()
    save(fig, outdir / "Figure2_anthropogenic_resource_expansion")


def figure3(mech: list[dict[str, str]], outdir: Path) -> None:
    rows = [
        r
        for r in mech
        if r["host_taxonomy_lower_bound_adequate"] == "True"
        and int(r["introduced_added_units"]) > 0
        and r["effective_contributor_number"]
        and r["maximum_single_host_fractional_share"]
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    boxplot_by_stratum(
        axes[0],
        rows,
        "effective_contributor_number",
        "Effective contributing host species",
    )
    axes[0].set_title("Generalists distribute expansion across more hosts")

    boxplot_by_stratum(
        axes[1],
        rows,
        "maximum_single_host_fractional_share",
        "Maximum single-host share",
    )
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_title("Specialists are more often dominated by one host")

    fig.tight_layout()
    save(fig, outdir / "Figure3_host_contribution_architecture")


def figure4(mech: list[dict[str, str]], outdir: Path) -> None:
    rows = [
        r
        for r in mech
        if r["host_taxonomy_lower_bound_adequate"] == "True"
        and int(r["introduced_added_units"]) > 0
        and float(r["host_family_count"]) == 1
        and r["effective_contributor_number"]
        and r["maximum_single_host_fractional_share"]
    ]
    host_species = np.asarray([int(r["resolved_host_species"]) for r in rows], float)
    effective = np.asarray([float(r["effective_contributor_number"]) for r in rows])
    dominance = np.asarray([float(r["maximum_single_host_fractional_share"]) for r in rows])

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    axes[0].scatter(host_species, effective, s=22, alpha=0.6)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Resolved host species within the single host family")
    axes[0].set_ylabel("Effective contributing host species")
    axes[0].set_title(f"Portfolio richness within one-family specialists\nSpearman ρ = {spearman(host_species, effective):.3f}")

    axes[1].scatter(host_species, dominance, s=22, alpha=0.6)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Resolved host species within the single host family")
    axes[1].set_ylabel("Maximum single-host share")
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_title(f"Single-host dominance weakens with portfolio richness\nSpearman ρ = {spearman(host_species, dominance):.3f}")

    fig.tight_layout()
    save(fig, outdir / "Figure4_within_family_specialization_hierarchy")


def figure5(primary: dict, outdir: Path) -> None:
    rows = primary["species"]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))

    score_values = []
    for i, stratum in enumerate(STRATA, start=1):
        subset = [r for r in rows if r["host_breadth_stratum"] == stratum]
        vals = [float(r["climate_filtering_score"]) for r in subset]
        score_values.append(vals)
        for row in subset:
            axes[0].scatter(
                i + deterministic_jitter(row["species"], 0.16),
                float(row["climate_filtering_score"]),
                s=28,
                alpha=0.7,
            )
        if vals:
            axes[0].plot([i - 0.22, i + 0.22], [np.median(vals)] * 2, linewidth=2)
    axes[0].axhline(0.5, linestyle="--", linewidth=1)
    axes[0].set_xticks(range(1, 5), [STRATA_LABELS[s] for s in STRATA], rotation=20)
    axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("Climate-filtering score")
    axes[0].set_title("Climate filtering is strong across host-breadth classes")

    x = np.asarray([float(r["host_family_count"]) for r in rows])
    y = np.asarray([float(r["climate_filtering_score"]) for r in rows])
    sizes = np.asarray([float(r["contemporary_host_resource_units"]) for r in rows])
    sizes = 20 + 80 * (np.log1p(sizes) - np.log1p(sizes).min()) / max(
        1e-9,
        np.log1p(sizes).max() - np.log1p(sizes).min(),
    )
    axes[1].scatter(
        [v + deterministic_jitter(r["species"], 0.12) for v, r in zip(x, rows)],
        y,
        s=sizes,
        alpha=0.65,
    )
    pt = primary["primary_test"]
    axes[1].set_xlabel("Larval host-plant families")
    axes[1].set_ylabel("Climate-filtering score")
    axes[1].set_ylim(0, 1.02)
    axes[1].set_title(
        "Broader diets do not weaken climate filtering\n"
        f"partial ρ = {float(pt['observed_partial_spearman']):.3f}, "
        f"one-sided p = {float(pt['one_sided_p_value']):.4f}"
    )

    fig.tight_layout()
    save(fig, outdir / "Figure5_independent_climate_filtering")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--descriptors-csv", type=Path, required=True)
    ap.add_argument("--anthropogenic-csv", type=Path, required=True)
    ap.add_argument("--mechanism-csv", type=Path, required=True)
    ap.add_argument("--climate-primary-json", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    descriptors = read_csv(args.descriptors_csv)
    anth = read_csv(args.anthropogenic_csv)
    mech = read_csv(args.mechanism_csv)
    primary = json.loads(args.climate_primary_json.read_text(encoding="utf-8"))
    if primary.get("status") != "INDEPENDENT_TEST_COMPLETE":
        raise RuntimeError("climate primary result is not complete")

    figure1(descriptors, args.output_dir)
    figure2(anth, args.output_dir)
    figure3(mech, args.output_dir)
    figure4(mech, args.output_dir)
    figure5(primary, args.output_dir)

    print(json.dumps({
        "figures": 5,
        "output_dir": str(args.output_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
