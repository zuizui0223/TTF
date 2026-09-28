import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_empirical_bc_detectability_is_calibration_qualified():
    result=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_bc_calibrated_detectability_result_v0.2.json").read_text()
    )
    final=json.loads(
        (ROOT/"benchmarks/frozen/relational_program_final_state_v0.3.json").read_text()
    )

    assert result["status"]=="FROZEN_RESPONSE_BLIND_CALIBRATED_BC_DETECTABILITY"
    assert final["final_state"]=="CLOSED_NO_EVALUABLE_TEST"
    assert result["boundary"]["relational_v03_final_state_unchanged"]==final["final_state"]
    assert result["boundary"]["relational_v03_reopened"] is False
    assert result["boundary"]["genetic_response_used"] is False

    current=result["current_climate"]["calibrated_evaluable_envelope"]
    historical=result["historical_given_current"]["calibrated_evaluable_envelope"]
    assert all(cell["calibration_pass"] for cell in current.values())
    assert all(cell["calibration_pass"] for cell in historical.values())
    assert all(cell["null_wilson95_upper"] <= 0.05 for cell in current.values())
    assert all(cell["null_wilson95_upper"] <= 0.05 for cell in historical.values())

    assert [current[k]["evaluable_grid_mde"] for k in ("A0","A1","A2","A3")] == [
        0.02,0.02,0.05,0.05
    ]
    assert [historical[k]["evaluable_grid_mde"] for k in ("A0","A1","A2","A3")] == [
        0.02,0.03,0.05,0.08
    ]
    assert result["interpretation"]["legacy_beta_0p03_A3_below_evaluable_mde_for_both"] is True
