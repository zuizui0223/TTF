import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from aggregate_historical_host_occurrences import sha256_path
from acquire_historical_host_occurrences import panel_hosts


def test_panel_hosts_deduplicates_hosts(tmp_path):
    p=tmp_path/"panel.csv"
    p.write_text(
        "species,n_hosts,accepted_host_ids,accepted_host_names\n"
        "A a,2,1;2,Plant one;Plant two\n"
        "B b,1,1,Plant one\n"
    )
    rows,hosts=panel_hosts(p)
    assert len(rows)==2
    assert hosts=={"1":"Plant one","2":"Plant two"}
