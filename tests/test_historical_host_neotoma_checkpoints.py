"""Response-blind fossil-query checkpoint contract; never inspect genetics."""
import sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from census_historical_host_neotoma import (
    checkpoint_name,read_completed_checkpoint,write_completed_checkpoint
)


def test_checkpoint_replays_exact_genus_dataset_age(tmp_path):
    result={
        "genus":"Quercus","dataset_type":"pollen",
        "records":2,"pages":1,"pairs":[["Quercus",10],["Quercus robur",20]],
    }
    written=write_completed_checkpoint(tmp_path,"Quercus","pollen",19000,26500,result)
    assert written.name==checkpoint_name("Quercus","pollen")
    restored=read_completed_checkpoint(tmp_path,"Quercus","pollen",19000,26500)
    assert restored==result
    assert read_completed_checkpoint(tmp_path,"Pinus","pollen",19000,26500) is None
    with pytest.raises(RuntimeError,match="source query drift"):
        read_completed_checkpoint(tmp_path,"Quercus","pollen",18000,26500)


def test_checkpoint_refuses_wrong_result_identity(tmp_path):
    result={"genus":"Betula","dataset_type":"pollen","pairs":[]}
    with pytest.raises(RuntimeError,match="does not match"):
        write_completed_checkpoint(tmp_path,"Quercus","pollen",19000,26500,result)


def test_corrupted_checkpoint_is_error_not_zero_records(tmp_path):
    result={"genus":"Quercus","dataset_type":"pollen","pairs":[]}
    written=write_completed_checkpoint(tmp_path,"Quercus","pollen",19000,26500,result)
    written.write_text('{"schema":"incorrect"}')
    with pytest.raises(RuntimeError,match="source query drift"):
        read_completed_checkpoint(tmp_path,"Quercus","pollen",19000,26500)
