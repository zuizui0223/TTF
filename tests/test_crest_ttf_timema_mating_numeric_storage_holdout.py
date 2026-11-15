"""Synthetic-only original-data storage semantics; no Timema outcome opened."""
from pathlib import Path
import importlib.util

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_timema_mating_numeric_storage_holdout_20261115.py"
spec=importlib.util.spec_from_file_location("ttf_timema_storage",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fake_rows():
    return [
      {"pair_id":"1","species_1":"A","species_2":"A","population_pair":"a1-a2",
       "sexual_cell_storage_class":"numeric_literal_storage_not_authenticated_trials"},
      {"pair_id":"2","species_1":"A","species_2":"B","population_pair":"a1-b1",
       "sexual_cell_storage_class":"shared_string_storage_not_decoded"},
      {"pair_id":"3","species_1":"B","species_2":"B","population_pair":"b1-b2",
       "sexual_cell_storage_class":"numeric_literal_storage_not_authenticated_trials"},
      {"pair_id":"4","species_1":"B","species_2":"C","population_pair":"b1-c1",
       "sexual_cell_storage_class":"numeric_literal_storage_not_authenticated_trials"},
      {"pair_id":"5","species_1":"C","species_2":"C","population_pair":"c1-c2",
       "sexual_cell_storage_class":"text_storage_not_decoded"},
      {"pair_id":"6","species_1":"C","species_2":"D","population_pair":"c1-d1",
       "sexual_cell_storage_class":"numeric_literal_storage_not_authenticated_trials"},
      {"pair_id":"7","species_1":"D","species_2":"D","population_pair":"d1-d2",
       "sexual_cell_storage_class":"numeric_literal_storage_not_authenticated_trials"},
    ]


def test_source_numeric_storage_is_not_all_serialized_cells():
    arm=m.count_by_arm(fake_rows(),["A","B"],["C","D"])
    assert arm["train"]["source_dyads"]==3
    assert arm["train"]["numeric_literal_dyads"]==2
    assert arm["train"]["numeric_literal_between"]==0
    assert arm["evaluation"]["numeric_literal_dyads"]==2
    assert arm["evaluation"]["numeric_literal_between"]==1
    assert arm["cross_boundary_dyads"]==1


def test_restricted_support_not_conflated_with_species_sample_size():
    result=m.evaluate_support(fake_rows())
    assert result["literal_numeric_cell_gate_partitions"]==0
    assert result["no_response_values_decoded"] is True
    assert result["numeric_storage_does_not_validate_an_original_mating_experiment"] is True
