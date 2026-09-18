from __future__ import annotations

import json
from pathlib import Path


LEDGER = Path("manuscript/genetic_ttf_submission_readiness_v0.1.json")
MANUSCRIPT = Path("manuscript/genetic_ttf_flagship_v0.2.md")


def test_genetic_submission_readiness_separates_science_from_human_metadata() -> None:
    ledger = json.loads(LEDGER.read_text())
    manuscript = MANUSCRIPT.read_text()

    assert ledger["schema"] == "ttf_genetic_submission_readiness_v0.1"
    assert ledger["status"] == "SCIENCE_COMPLETE_HUMAN_METADATA_AND_JOURNAL_FORMATTING_PENDING"

    science = ledger["frozen_science"]
    assert science["decision"] == "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN"
    assert science["survivor_species"] == 211
    assert science["training_species"] == 103
    assert science["evaluation_species"] == 108
    assert science["primary_place_beyond_ibd"]["statistic"] == 0.03256548471362711
    assert science["primary_place_beyond_ibd"]["profiled_private_p_value"] == 0.7392607392607392
    assert science["primary_place_beyond_ibd"]["positive"] is False
    assert science["within_species_self"]["upper_tail_p_value"] == 0.008991008991008992
    assert science["within_species_self"]["raw_sign_interpreted_against_zero"] is False
    assert science["secondary_total_genetic_transfer"]["inferentially_qualified"] is False
    assert science["one_shot_completed"] is True
    assert science["empirical_rerun_allowed"] is False
    assert science["post_result_retuning_allowed"] is False

    package = ledger["frozen_artifacts"]
    assert package["figure_count"] == 2
    assert package["supplementary_evaluation_species"] == 108
    for path in package["figures"]:
        assert Path(path).is_file()
    for key in ("supplementary_results", "supplementary_species_table", "supplementary_manifest", "terminal_handoff"):
        assert Path(package[key]).is_file()

    manuscript_state = ledger["manuscript"]
    assert manuscript_state["author_list_finalized"] is False
    assert manuscript_state["target_journal_selected"] is False
    assert manuscript_state["abstract_word_count"] == 260
    assert "Author list: TBD" in manuscript
    assert "Target journal: deferred" in manuscript
    assert "EMPIRICAL RESULT SLOT" not in manuscript
    assert "FINAL CONCLUSION SLOT" not in manuscript

    assert "proof of zero transfer" in ledger["claim_boundary"]
