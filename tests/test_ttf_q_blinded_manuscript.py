import runpy
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/build_ttf_q_blinded_manuscript.py"
blind=runpy.run_path(str(SCRIPT))["blind"]


def test_ttf_q_blinded_manuscript_is_deterministic_and_has_no_public_repo_identity():
    source=(ROOT/"manuscript/ttf_q_methods_v0.1.md").read_text()
    committed=(ROOT/"manuscript/ttf_q_methods_blinded_v0.1.md").read_text()
    assert committed == blind(source)
    lower=committed.lower()
    assert "github.com/zuizui0223" not in lower
    assert "zuizui0223/ttf" not in lower
    assert "omitted for double-anonymous review" in lower
    assert "an anonymized peer-review archive" in lower
