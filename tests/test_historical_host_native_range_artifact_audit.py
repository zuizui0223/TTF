"""Response-blind exact native-range identity binding (does not open genetics)."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FREEZE=ROOT/"benchmarks/frozen"


def test_exact_135_species_222_host_id_binding():
    panel=FREEZE/"genetic_historical_host_connectivity_panel_v0.2.csv"
    receipt=json.loads((FREEZE/"genetic_historical_host_native_range_artifact_audit_v0.1.json").read_text())
    with panel.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    ids={host for row in rows for host in row["accepted_host_ids"].split(";") if host}

    assert len(rows)==135
    assert len(ids)==222
    assert hashlib.sha256(panel.read_bytes()).hexdigest()==receipt["panel"]["sha256"]
    digest=hashlib.sha256(("\n".join(sorted(ids))+"\n").encode()).hexdigest()
    assert digest==receipt["audit"]["panel_host_id_list_sha256"]
    assert digest==receipt["audit"]["host_id_list_sha256"]
    assert receipt["audit"]["distinct_geojson_accepted_host_ids"]==222
    assert receipt["audit"]["empty_geometry_count"]==0
    assert receipt["frozen_transport"]["artifact_id"]==11519874583
    assert receipt["obsolete_artifact_warning"]["artifact_id"]!=receipt["frozen_transport"]["artifact_id"]
    assert receipt["can_authorize_nucleotide_opening"] is False
    assert all(x is False for x in receipt["response_firewall"].values())
