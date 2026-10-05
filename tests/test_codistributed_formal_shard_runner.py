from pathlib import Path
import subprocess
import sys


ROOT=Path(__file__).resolve().parents[1]


def test_formal_shard_runner_cli_imports_cleanly():
    runner=ROOT/"scripts"/"run_codistributed_recurrence_fixed_tst_formal_shard.py"
    completed=subprocess.run(
        [sys.executable,str(runner),"--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode==0, completed.stderr
    assert "--component" in completed.stdout
    assert "--sharding" in completed.stdout
