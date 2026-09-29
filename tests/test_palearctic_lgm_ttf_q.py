import runpy
from dataclasses import make_dataclass
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
NS=runpy.run_path(str(ROOT/"scripts/run_palearctic_lgm_ttf_q.py"))


Info=make_dataclass("Info",[("total_unique_variance_fraction",float)])
Support=make_dataclass(
    "Support",
    [
        ("signal_effective_sources",float),
        ("signal_effective_targets",float),
        ("max_source_signal_share",float),
        ("max_target_signal_share",float),
    ],
)


RULE={
    "minimum_total_unique_variance_fraction":0.15,
    "minimum_signal_effective_sources":10,
    "minimum_signal_effective_targets":10,
    "maximum_single_source_signal_share":0.15,
    "maximum_single_target_signal_share":0.15,
    "required_null_qualified_amplitudes":[0,1,2],
    "required_A2_evaluable_mde_max":0.10,
}


def env(mde2=0.08,pass2=True):
    return {
        "0":{"calibration_pass":True,"evaluable_grid_mde":0.05},
        "1":{"calibration_pass":True,"evaluable_grid_mde":0.08},
        "2":{"calibration_pass":pass2,"evaluable_grid_mde":mde2 if pass2 else None},
        "3":{"calibration_pass":False,"evaluable_grid_mde":None},
    }


def test_opening_rule_passes_without_requiring_extreme_A3():
    gates=NS["evaluate_opening"](
        Info(0.20),
        Support(14,13,0.11,0.12),
        env(),
        RULE,
    )
    assert gates["overall_pass"] is True


def test_opening_rule_fails_for_redundant_or_unqualified_relation():
    assert not NS["evaluate_opening"](
        Info(0.10),Support(14,13,0.11,0.12),env(),RULE
    )["overall_pass"]
    assert not NS["evaluate_opening"](
        Info(0.20),Support(14,13,0.11,0.12),env(pass2=False),RULE
    )["overall_pass"]
    assert not NS["evaluate_opening"](
        Info(0.20),Support(14,13,0.11,0.12),env(mde2=0.15),RULE
    )["overall_pass"]
