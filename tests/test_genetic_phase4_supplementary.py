from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


EXPORT_DIR = Path("manuscript/generated/genetic_ttf_phase4_v0.1")
HANDOFF = Path("benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_phase4_species_table_matches_terminal_handoff() -> None:
    manifest = json.loads((EXPORT_DIR / "export_manifest.json").read_text())
    handoff = json.loads(HANDOFF.read_text())
    rows = list(csv.DictReader((EXPORT_DIR / "species_scores.csv").open(newline="")))

    assert manifest["schema"] == "ttf_genetic_manuscript_export_v0.1"
    assert manifest["mode"] == "empirical_summary"
    assert manifest["decision"] == handoff["empirical_result"]["decision"]
    assert manifest["evaluation_species_count"] == len(rows) == 108
    assert manifest["training_species_count"] == 103
    assert manifest["geometry_fingerprint_sha256"] == handoff["phase4_authorization"]["geometry_fingerprint_sha256"]
    assert manifest["receipt_consistency_checked"] is True
    assert manifest["statistical_analysis_rerun"] is False
    assert manifest["sequence_identity_read"] is False
    assert manifest["input_interpretation_text_used"] is False
    assert manifest["source_sha256"]["result"] == handoff["empirical_result"]["result_json_sha256"]
    assert manifest["source_sha256"]["authorization"] == handoff["phase4_authorization"]["authorization_json_sha256"]

    for name, expected in manifest["output_sha256"].items():
        assert _sha256(EXPORT_DIR / name) == expected

    species = [row["species"] for row in rows]
    assert len(species) == len(set(species)) == 108
    assert species == sorted(species)
    assert all(row["total_available"] == "True" for row in rows)

    primary = [float(row["post_ibd_transfer"]) for row in rows]
    self_score = [float(row["within_species_self"]) for row in rows]
    total = [float(row["total_transfer_descriptive"]) for row in rows]
    assert all(math.isfinite(value) for value in primary + self_score + total)

    empirical = handoff["empirical_result"]
    assert math.isclose(
        math.fsum(primary) / len(primary),
        empirical["primary_place_beyond_ibd"]["statistic"],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert math.isclose(
        math.fsum(self_score) / len(self_score),
        empirical["within_species_self_diagnostic"]["statistic"],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert math.isclose(
        math.fsum(total) / len(total),
        empirical["secondary_total_genetic_transfer"]["statistic"],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )


def test_frozen_phase4_generated_results_keep_interpretation_boundary() -> None:
    text = (EXPORT_DIR / "results.md").read_text()
    assert "profiled-private p = 0.739261" in text
    assert "upper-tail p = 0.00899101" in text
    assert "raw numerical sign is not interpreted against zero" in text
    assert "It has no significance test and cannot rescue the primary result." in text
    assert "zero transfer" in text
