import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def git_blob_sha(path: Path) -> str:
    raw = path.read_bytes()
    payload = f"blob {len(raw)}\0".encode() + raw
    return hashlib.sha1(payload).hexdigest()


def test_ttf_q_submission_scientific_freeze_binds_current_science():
    freeze = json.loads(
        (ROOT / "benchmarks/frozen/ttf_q_submission_scientific_freeze_v1.json").read_text()
    )
    assert freeze["status"] == "SCIENTIFIC_CONTENT_FROZEN_AWAITING_AUTHOR_INPUT_ONLY"
    assert freeze["source_main_sha"] == "944d6a06311360107f85f4c62c39ab5fd549ff40"

    bindings = {}
    bindings.update(freeze["manuscript_bindings"])
    bindings.update(freeze["method_and_case_bindings"])
    for rel, meta in bindings.items():
        path = ROOT / rel
        assert path.is_file(), rel
        assert git_blob_sha(path) == meta["git_blob_sha"], rel

    science = freeze["scientific_state"]
    assert science["prospective_palearctic_route"] == "CLOSED_NOT_EVALUABLE_WITHOUT_GENETIC_RESPONSE"
    assert science["subpanel_nucleotide_identity_opened"] is False
    assert science["pairwise_T_st_computed"] is False
    assert science["beta_LGM_computed"] is False
    assert science["replacement_subgroup_or_predictor_authorized"] is False

    validated = freeze["validated_execution"]
    assert validated["tests"]["pytest"] == "485 passed"
    assert validated["figures"]["figure_count"] == 5
    assert validated["submission_readiness"]["scientific_package_ready"] is True
    assert validated["submission_readiness"]["structural_failures"] == []
    assert validated["current_manuscript_renderer_smoke"]["conclusion"] == "success"

    remaining = freeze["remaining_author_input"]
    assert remaining["hard_blockers"] == [
        "choose open-source license",
        "complete title-page author/affiliation/correspondence/acknowledgement/contribution/conflict fields",
    ]

    policy = freeze["freeze_policy"]
    assert "new analysis" in policy["requires_explicit_scientific_unfreeze"]
    assert "new subgroup" in policy["requires_explicit_scientific_unfreeze"]
    assert "new predictor" in policy["requires_explicit_scientific_unfreeze"]
    assert "opening a previously closed outcome" in policy["requires_explicit_scientific_unfreeze"]
