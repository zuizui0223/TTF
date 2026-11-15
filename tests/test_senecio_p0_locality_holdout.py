"""Offline synthetic-only safety checks for exploratory Senecio P0 scorer."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "research" / "prospective" / "senecio_p0_pair_holdout_20261115.py"
sp = importlib.util.spec_from_file_location("senecio_p0", PATH)
m = importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def fake_pair_means(missing=False, duplicate=False, nonfinite=False) -> bytes:
    fields = ["Group.1", *m.TRAITS, "Ecotype", "Pair"]
    lines = ["\t".join(fields)]
    for i in range(8):
        pair = f"D{i:02d}-H{i:02d}"
        for eco in ("Dune", "Headland"):
            vals = ["0" if eco == "Dune" else str(1.0 + i / 100) for _ in m.TRAITS]
            if nonfinite and i == 0 and eco == "Dune":
                vals[0] = "nan"
            if missing and i == 7 and eco == "Headland":
                continue
            lines.append("\t".join([str(i), *vals, eco, pair]))
            if duplicate and i == 0 and eco == "Headland":
                lines.append(lines[-1])
    return ("\n".join(lines) + "\n").encode()


def run_with_mock_sha(raw: bytes):
    old = m.SHA
    m.SHA = m.git_blob_sha(raw)
    try:
        return m.score(raw)
    finally:
        m.SHA = old


def test_full_eight_pair_source_grain_and_per_trait_scoring():
    r = run_with_mock_sha(fake_pair_means())
    assert r["n_independent_evaluation_units"] == 8
    assert r["n_phenotype_traits"] == 9
    assert len(r["trait_results"]) == 9
    assert r["summary_descriptive"]["traits_with_positive_heldout_skill"] == 9


def test_unpinned_source_rejected():
    try:
        m.score(fake_pair_means())
    except ValueError as e:
        assert "SHA mismatch" in str(e)
    else:
        raise AssertionError("wrong source SHA accepted")


def test_missing_ecotype_and_duplicate_pairs_rejected():
    for raw in (fake_pair_means(missing=True), fake_pair_means(duplicate=True)):
        try:
            run_with_mock_sha(raw)
        except ValueError:
            pass
        else:
            raise AssertionError("bad pair/ecotype keys accepted")


def test_nonfinite_values_rejected():
    try:
        run_with_mock_sha(fake_pair_means(nonfinite=True))
    except ValueError as e:
        assert "non-finite" in str(e)
    else:
        raise AssertionError("non-finite outcome accepted")
