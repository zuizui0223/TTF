from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


FIGURES = Path("manuscript/figures/genetic_phase4_v0.1")
SCORES = Path("manuscript/generated/genetic_ttf_phase4_v0.1/species_scores.csv")
HANDOFF = Path("benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json")
MANUSCRIPT = Path("manuscript/genetic_ttf_flagship_v0.2.md")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_phase4_figures_are_receipt_backed_descriptions() -> None:
    manifest = json.loads((FIGURES / "figure_manifest.json").read_text())
    handoff = json.loads(HANDOFF.read_text())
    rows = list(csv.DictReader(SCORES.open(newline="")))

    assert manifest["schema"] == "ttf_genetic_phase4_figure_manifest_v0.1"
    assert manifest["status"] == "DESCRIPTIVE_VISUALIZATION_ONLY"
    assert manifest["new_inference_performed"] is False
    assert manifest["post_result_retuning_performed"] is False
    assert manifest["evaluation_species_count"] == len(rows) == 108

    assert manifest["input_sha256"]["species_scores.csv"] == _sha256(SCORES)
    assert manifest["input_sha256"]["phase4_empirical_handoff.json"] == _sha256(HANDOFF)
    assert (
        manifest["terminal_empirical_result_sha256"]
        == handoff["empirical_result"]["result_json_sha256"]
    )

    primary = [float(row["post_ibd_transfer"]) for row in rows]
    expected_primary = handoff["empirical_result"]["primary_place_beyond_ibd"]["statistic"]
    assert math.isclose(
        math.fsum(primary) / len(primary),
        manifest["primary_species_equal_mean"],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert math.isclose(
        manifest["primary_species_equal_mean"],
        expected_primary,
        rel_tol=1e-12,
        abs_tol=1e-12,
    )

    cross = handoff["formal_phase3_gate_d"]
    self_gate = handoff["self_detectability"]
    qualification = manifest["qualification"]
    assert qualification["cross_species_type1_wilson_upper"] == cross["max_private_null_wilson95_upper"]
    assert qualification["self_type1_wilson_upper"] == self_gate["null_wilson95_upper"]
    assert qualification["cross_species_power_wilson_lower"] == cross["shared_A2_wilson95_lower"]
    assert qualification["self_power_wilson_lower"] == self_gate["private_A2_wilson95_lower"]
    assert qualification["type1_ceiling"] == cross["type1_ceiling"] == 0.10
    assert qualification["power_floor"] == cross["power_floor"] == 0.80

    for name, expected in manifest["output_sha256"].items():
        path = FIGURES / name
        assert path.is_file()
        assert _sha256(path) == expected
        text = path.read_text()
        assert "<svg" in text
        assert "matplotlib" in text.lower()


def test_empirical_manuscript_uses_figures_without_new_inference() -> None:
    text = MANUSCRIPT.read_text()
    assert "figures/genetic_phase4_v0.1/figure1_qualification_margins.svg" in text
    assert "figures/genetic_phase4_v0.1/figure2_post_ibd_species_scores.svg" in text
    assert "no species-level significance tests" in text
    assert "display ordering" in text
    assert "new inferential analysis" not in text
