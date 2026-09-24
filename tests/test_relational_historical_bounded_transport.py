from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FALSE_FIREWALL = {
    "Study_C_sequence_identity_opened": False,
    "Study_C_pairwise_genetic_distances_opened": False,
    "Study_C_T_st_computed": False,
    "Study_C_beta_hist_computed": False,
}
FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def load_script(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def candidate_names() -> list[str]:
    return [f"Species {i:04d}" for i in range(1000)]


def write_candidates(path: Path, names: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["species"])
        writer.writeheader()
        for name in names:
            writer.writerow({"species": name})


def write_initial(input_dir: Path, names: list[str], statuses: dict[str, str]) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    for batch in range(250):
        selected = [name for i, name in enumerate(names) if i % 250 == batch]
        ledger = {
            "schema": "ttf_relational_historical_occurrence_bounded_batch_v0.2",
            "batch_index": batch,
            "batches": 250,
            "candidate_species_total": 1000,
            "species_requested": selected,
            "species": [
                {"species": name, "status": statuses.get(name, "FAIL_OCCURRENCE_GEOMETRY")}
                for name in selected
            ],
            "scientific_query_change": False,
            "response_firewall": dict(FALSE_FIREWALL),
        }
        (input_dir / f"ledger-bounded-{batch}.json").write_text(
            json.dumps(ledger), encoding="utf-8"
        )
        with (input_dir / f"occurrences-bounded-{batch}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            for ordinal, name in enumerate(selected):
                if statuses.get(name) == "PASS_OCCURRENCE_GEOMETRY":
                    writer.writerow(
                        {
                            "species": name,
                            "source_key": batch * 10 + ordinal,
                            "latitude": 0.0,
                            "longitude": 0.0,
                            "priority_rank": 0,
                            "priority_sha256": "0" * 64,
                        }
                    )


def write_plan(path: Path, request_errors: list[str]) -> None:
    payload = {
        "schema": "ttf_relational_historical_occurrence_retry_plan_v0.2",
        "status": "RETRY_EXACT_INITIAL_REQUEST_ERRORS" if request_errors else "NO_RETRY_NEEDED",
        "candidate_species": 1000,
        "request_error_count": len(request_errors),
        "request_error_species": request_errors,
        "retry_batch_size": 4,
        "retry_batch_count": (len(request_errors) + 3) // 4,
        "maximum_retry_rounds": 1,
        "scientific_query_change": False,
        "response_firewall": dict(FALSE_FIREWALL),
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_retry(input_dir: Path, species_rows: list[tuple[str, str]]) -> None:
    input_dir.mkdir(parents=True, exist_ok=True)
    ledger = {
        "schema": "ttf_relational_historical_occurrence_bounded_retry_batch_v0.2",
        "retry_round": 1,
        "batch_index": 0,
        "batch_size": 4,
        "species_requested": [name for name, _ in species_rows],
        "species": [{"species": name, "status": status} for name, status in species_rows],
        "retained_occurrence_rows": 0,
        "scientific_query_change": False,
        "response_firewall": dict(FALSE_FIREWALL),
    }
    (input_dir / "ledger-retry-0.json").write_text(json.dumps(ledger), encoding="utf-8")
    with (input_dir / "occurrences-retry-0.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        csv.DictWriter(handle, fieldnames=FIELDS).writeheader()


def run_finalize(monkeypatch, base: Path, statuses: dict[str, str], request_errors: list[str], retry_rows: list[tuple[str, str]]):
    module = load_script("hist_finalize_bounded", "scripts/finalize_relational_historical_occurrence_bounded.py")
    names = candidate_names()
    candidates = base / "candidates.csv"
    initial = base / "initial"
    retry = base / "retry"
    plan = base / "plan.json"
    output_csv = base / "out.csv"
    output_ledger = base / "out.json"
    base.mkdir(parents=True, exist_ok=True)
    write_candidates(candidates, names)
    write_initial(initial, names, statuses)
    write_plan(plan, request_errors)
    if retry_rows:
        write_retry(retry, retry_rows)
    else:
        retry.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "finalize_relational_historical_occurrence_bounded.py",
            "--initial-dir",
            str(initial),
            "--retry-dir",
            str(retry),
            "--retry-plan",
            str(plan),
            "--candidates",
            str(candidates),
            "--output-csv",
            str(output_csv),
            "--output-ledger",
            str(output_ledger),
        ],
    )
    module.main()
    return json.loads(output_ledger.read_text())


def test_bounded_partition_matches_frozen_modulo_250_contract():
    module = load_script("hist_bounded", "scripts/acquire_relational_historical_occurrences_bounded.py")
    names = candidate_names()
    seen: list[str] = []
    for batch in range(250):
        selected = list(module.shard_items(names, shard_index=batch, shards=250))
        assert selected == [names[batch], names[batch + 250], names[batch + 500], names[batch + 750]]
        seen.extend(selected)
    assert len(seen) == 1000
    assert set(seen) == set(names)


def test_retry_planner_selects_only_request_errors_in_frozen_candidate_order(tmp_path, monkeypatch):
    module = load_script("hist_retry_plan", "scripts/plan_relational_historical_occurrence_retry.py")
    names = candidate_names()
    candidates = tmp_path / "candidates.csv"
    initial = tmp_path / "initial"
    output = tmp_path / "plan.json"
    write_candidates(candidates, names)
    statuses = {name: "PASS_OCCURRENCE_GEOMETRY" for name in names}
    statuses[names[501]] = "REQUEST_ERROR"
    statuses[names[999]] = "REQUEST_ERROR"
    write_initial(initial, names, statuses)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plan_relational_historical_occurrence_retry.py",
            "--input-dir",
            str(initial),
            "--candidates",
            str(candidates),
            "--output",
            str(output),
        ],
    )
    module.main()
    plan = json.loads(output.read_text())
    assert plan["request_error_species"] == [names[501], names[999]]
    assert plan["request_error_count"] == 2
    assert plan["maximum_retry_rounds"] == 1


def test_unresolved_request_error_is_terminal_technical_not_a_biological_null(tmp_path, monkeypatch):
    names = candidate_names()
    statuses = {name: "PASS_OCCURRENCE_GEOMETRY" for name in names}
    statuses[names[1]] = "REQUEST_ERROR"
    result = run_finalize(
        monkeypatch,
        tmp_path / "technical",
        statuses,
        [names[1]],
        [(names[1], "REQUEST_ERROR")],
    )
    assert result["species_passing_ge_30"] == 999
    assert result["final_request_error_count"] == 1
    assert result["status"] == "NOT_EVALUABLE_HISTORICAL_TECHNICAL_TRANSPORT"
    assert result["additional_retry_authorized"] is False


@pytest.mark.parametrize(
    ("n_pass", "expected"),
    [
        (499, "NOT_EVALUABLE_HISTORICAL_OCCURRENCE_GEOMETRY"),
        (500, "PASS_TO_HISTORICAL_ASSET_EXTRACTION"),
    ],
)
def test_resolved_transport_uses_exact_500_species_gate(tmp_path, monkeypatch, n_pass, expected):
    names = candidate_names()
    statuses = {
        name: ("PASS_OCCURRENCE_GEOMETRY" if i < n_pass else "FAIL_OCCURRENCE_GEOMETRY")
        for i, name in enumerate(names)
    }
    result = run_finalize(
        monkeypatch,
        tmp_path / f"gate-{n_pass}",
        statuses,
        [],
        [],
    )
    assert result["final_request_error_count"] == 0
    assert result["species_passing_ge_30"] == n_pass
    assert result["status"] == expected


def test_finalizer_rejects_retrying_species_outside_initial_request_error_set(tmp_path, monkeypatch):
    names = candidate_names()
    statuses = {name: "PASS_OCCURRENCE_GEOMETRY" for name in names}
    statuses[names[1]] = "REQUEST_ERROR"
    with pytest.raises(RuntimeError, match="retry species set differs"):
        run_finalize(
            monkeypatch,
            tmp_path / "wrong-retry",
            statuses,
            [names[1]],
            [(names[2], "PASS_OCCURRENCE_GEOMETRY")],
        )
