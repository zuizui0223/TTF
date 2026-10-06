from pathlib import Path
import subprocess
import sys

import numpy as np

from ttf.codistributed_survivor import evaluate_survivor_information_gate


ROOT=Path(__file__).resolve().parents[1]


def _balanced_panel(seed=7):
    rng=np.random.default_rng(seed)
    ns=20
    nt=20
    s=np.repeat(np.arange(ns),nt)
    t=np.tile(np.arange(nt),ns)
    g=rng.normal(size=len(s))
    controls=np.column_stack([
        rng.normal(size=len(s)),
        rng.normal(size=len(s)),
        rng.normal(size=len(s)),
        rng.integers(0,2,size=len(s)).astype(float),
    ])
    return s,t,g,controls


def test_survivor_information_gate_passes_broad_pair_specific_signal():
    s,t,g,c=_balanced_panel()
    result=evaluate_survivor_information_gate(s,t,g,c)
    assert result.total_unique_variance_fraction >= 0.15
    assert result.max_source_signal_share <= 0.15
    assert result.max_target_signal_share <= 0.15
    assert result.overall_pass is True


def test_survivor_information_gate_rejects_endpoint_concentrated_signal():
    s,t,g,c=_balanced_panel()
    g=np.zeros_like(g)
    rng=np.random.default_rng(9)
    g[s==0]=rng.normal(size=np.count_nonzero(s==0))
    g[s!=0]=1e-6*rng.normal(size=np.count_nonzero(s!=0))
    result=evaluate_survivor_information_gate(s,t,g,c)
    assert result.max_source_signal_share > 0.15
    assert result.source_concentration_pass is False
    assert result.overall_pass is False


def test_survivor_information_runner_cli_imports_cleanly():
    runner=ROOT/"scripts"/"run_codistributed_recurrence_survivor_information.py"
    completed=subprocess.run(
        [sys.executable,str(runner),"--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode==0,completed.stderr
    assert "--mask-result" in completed.stdout
    assert "--survivors" in completed.stdout
