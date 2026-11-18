"""Synthetic-only NCBI original source submitter alias-voucher matching tests."""
import importlib.util
from pathlib import Path
import json
import pytest

PATH=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_ncbi_sourcesample_alias_pilot_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_5_voucher_alias",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def source():
    samn=[f"SAMN6099293{i}" for i in range(1,6)]
    aliases=["cardBL_1_S443","CAW123_1_S454","XX0","",None]
    uids=list(map(str,range(1,6)))
    items={"uids":uids}
    for j,(k,s) in enumerate(zip(uids,samn)):
        items[k]={"accession":s,"organism":"Penstemon barbatus",
                  "sourcesample":aliases[j],"identifiers":{"other":"data"},
                  "sampledata":"SENTINEL_UNREQUESTED_BIOLOGICAL_VALUES"}
    return samn,{"result":items}


def test_only_short_submitter_original_key_parsed_and_exact_tip_matched():
    ids,raw=source()
    r=m.inspect_alias_fields(json.dumps(raw).encode(),ids,
        {"cardBL_1_S443":"cardinalis","CAW123_1_S454":"eatonii"})
    assert r["n_fixed_source_ids"]==5
    assert r["n_source_submitter_keys_available"]==3
    assert r["n_exact_old_2023_mapped_tips"]==2
    assert "SENTINEL_UNREQUESTED" not in repr(r)
    assert r["not_full_2026_to_2023_genotype_identity"] is True


def test_foreign_biosample_uid_join_fails_closed():
    ids,raw=source()
    raw["result"]["5"]["accession"]="SAMN00000"
    with pytest.raises(ValueError,match="outside fixed five source IDs"):
        m.inspect_alias_fields(json.dumps(raw).encode(),ids,{})


def test_older_groupmap_contract_and_fixed_five_samples():
    assert m.SOURCE_BLOB=="360d904346db4ebc17b13d559aed5cd13dde2867"
    assert m.N==5
