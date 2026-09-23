import json
import subprocess
import sys
from pathlib import Path


def test_relational_v03_closure_verifier_passes(tmp_path: Path):
    output = tmp_path / "closure.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/verify_relational_program_v03_closure.py",
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        "closure verifier failed\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    payload = json.loads(output.read_text())
    assert payload["schema"] == "ttf_relational_program_v03_closure_audit_v0.1"
    assert payload["status"] == "PASS_FROZEN_V03_CONTRACT_CLOSURE"
    assert payload["prospective_slots"] == ["B", "C"]
    assert payload["slot_alpha"] == {"B": 0.025, "C": 0.025}
    assert payload["current_S3"]["p_one_sided"] == 0.6376882940222616
    assert payload["current_S3"]["species_union"] == 399
    assert payload["current_S3"]["dyads"] == 11696
    assert payload["study_B_full_downstream_contract_frozen"] is True
    assert payload["study_B_final_transport_retry_round"] == 4
    assert payload["study_B_no_fifth_transport_retry"] is True
    assert payload["study_B_final_relation_producer_prefrozen"] is True
    assert payload["study_B_final_transport_audit_only_prefrozen"] is True
    assert payload["study_C_full_downstream_contract_frozen"] is True
    assert payload["all_future_genetic_response_firewalls_closed"] is True
    assert "PASS_FROZEN_V03_CONTRACT_CLOSURE" in completed.stdout
