"""Mask-only inherited phylogatR parser tests; no empirical response is used."""
from pathlib import Path

import pytest

from scripts.freeze_historical_host_memory_character_masks import (
    read_frozen_geo, safe_raw_path, score_species
)


HEADER="phylogatr_id\taccession\tsource_id\tlatitude\tlongitude\tbasis_of_record\tcoordinate_uncertainty_in_meters\tissue\tflag\n"


def fixture(tmp_path: Path, mutate=None):
    d=tmp_path/"Insecta"/"Taxon"
    d.mkdir(parents=True)
    seq=["ACGTACGT","ACGTACGT","ACGTACGT","ACGTACGT"]
    if mutate is not None:
        seq[2]=mutate
    fasta=d/"Taxon-COI.afa"
    fasta.write_text("".join(f">h{i}\n{x}\n" for i,x in enumerate(seq)))
    (d/"occurrences.txt").write_text(HEADER+"".join(
        f"h{i}\ta{i}\tsrc\t{i}.0\t{i+2}.0\tOBS\t0\t\t\n"
        for i in range(4)
    ))
    row={"species":"Taxon","raw_dir":"Insecta/Taxon","raw_gene":"Taxon-COI",
         "n_localities":"4","n_headers":"4","edges":"3"}
    ll={i:(float(i),float(i+2)) for i in range(4)}
    edges={0:(0,1),1:(1,2),2:(2,3)}
    return row,ll,edges


def test_frozen_mask_parser_retains_only_boolean_validity(tmp_path):
    row,ll,edges=fixture(tmp_path)
    got=score_species(row,tmp_path,ll,edges)
    assert got["status"]=="PASS_MASK"
    assert got["valid_edges"]==3
    assert got["invalid_edges"]==0
    assert got["minimum_comparable_columns"]==4
    assert got["nucleotide_identity_persisted"] is False
    assert "ACGT" not in repr(got)


def test_failed_edge_drops_species_not_edges(tmp_path):
    row,ll,edges=fixture(tmp_path,mutate="NNNNNNNN")
    got=score_species(row,tmp_path,ll,edges)
    assert got["status"]=="NOT_EVALUABLE_CHARACTER_SUPPORT"
    assert got["frozen_edges"]==3
    assert got["invalid_edges"]>=1
    assert got["survives"] is False


def test_missing_file_and_path_traversal_are_fatal(tmp_path):
    with pytest.raises(RuntimeError,match="path traversal"):
        safe_raw_path(tmp_path,"../outside","x.afa")
    with pytest.raises(FileNotFoundError,match="exact raw source absent"):
        score_species({"species":"X y","raw_dir":"none","raw_gene":"fake",
                       "n_localities":"2","edges":"1","n_headers":"2"},
                      tmp_path,{0:(0.,0.),1:(1.,1.)},{0:(0,1)})


def test_duplicate_header_fails_species_only(tmp_path):
    row,ll,edges=fixture(tmp_path)
    path=tmp_path/"Insecta"/"Taxon"/"Taxon-COI.afa"
    path.write_text(">h0\nACGT\n>h0\nACGT\n")
    got=score_species(row,tmp_path,ll,edges)
    assert got["status"]=="NOT_EVALUABLE_CHARACTER_SUPPORT"
    assert got["survives"] is False


def test_unfrozen_geo_structure_is_rejected(tmp_path):
    row,ll,edges=fixture(tmp_path)
    with pytest.raises(RuntimeError,match="noncontiguous"):
        score_species(row,tmp_path,{2:(0.,1.),3:(2.,3.)},edges)
    with pytest.raises(RuntimeError,match="graph geometry changed"):
        score_species(row,tmp_path,ll,{0:(0,1)})
