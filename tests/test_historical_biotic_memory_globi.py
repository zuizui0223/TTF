from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts"/"summarize_historical_biotic_memory_globi.py"
spec=importlib.util.spec_from_file_location("hibm_globi",SCRIPT)
m=importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(m)


RULE={
    "accepted_orientations":{
        "focal_as_source":{"interaction_types":["hasHost","eats"]},
        "focal_as_target":{"interaction_types":["hostOf","eatenBy"]},
    }
}


def row(**kw):
    out={
        "interaction_type":"",
        "source_taxon_name":"",
        "target_taxon_name":"",
        "source_taxon_path":"",
        "target_taxon_path":"",
        "study_title":"study",
    }
    out.update(kw)
    return out


def test_species_level_name_boundary():
    assert m.species_level_name("Quercus robur")
    assert m.species_level_name("Capparis spinosa subsp. cordifolia")
    assert not m.species_level_name("Quercus")
    assert not m.species_level_name("Quercus sp.")
    assert not m.species_level_name("Fabaceae host")


def test_forward_and_inverse_host_relations_collapse_to_same_plant():
    forward=row(
        interaction_type="eats",
        source_taxon_name="Butterfly alpha",
        target_taxon_name="Quercus robur",
        target_taxon_path="Eukaryota | Plantae | Fagaceae | Quercus | Quercus robur",
    )
    inverse=row(
        interaction_type="eatenBy",
        source_taxon_name="Quercus robur",
        target_taxon_name="Butterfly alpha",
        source_taxon_path="Eukaryota | Plantae | Fagaceae | Quercus | Quercus robur",
    )
    assert m.classify_row(forward,focal="Butterfly alpha",family="Fagaceae",rule=RULE)==(
        "accepted","Quercus robur"
    )
    assert m.classify_row(inverse,focal="Butterfly alpha",family="Fagaceae",rule=RULE)==(
        "accepted","Quercus robur"
    )


def test_host_family_mismatch_and_visit_are_rejected():
    mismatch=row(
        interaction_type="hasHost",
        source_taxon_name="Butterfly alpha",
        target_taxon_name="Salix alba",
        target_taxon_path="Eukaryota | Plantae | Salicaceae | Salix | Salix alba",
    )
    visit=row(
        interaction_type="visitsFlowersOf",
        source_taxon_name="Butterfly alpha",
        target_taxon_name="Quercus robur",
        target_taxon_path="Eukaryota | Plantae | Fagaceae | Quercus | Quercus robur",
    )
    assert m.classify_row(mismatch,focal="Butterfly alpha",family="Fagaceae",rule=RULE)==(
        "host_family_mismatch",None
    )
    assert m.classify_row(visit,focal="Butterfly alpha",family="Fagaceae",rule=RULE)==(
        "orientation_or_verb",None
    )


def test_genus_only_host_is_audit_only():
    genus=row(
        interaction_type="hasHost",
        source_taxon_name="Butterfly alpha",
        target_taxon_name="Quercus",
        target_taxon_path="Eukaryota | Plantae | Fagaceae | Quercus",
    )
    assert m.classify_row(genus,focal="Butterfly alpha",family="Fagaceae",rule=RULE)==(
        "not_species_level_plant",None
    )
