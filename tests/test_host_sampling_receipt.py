"""Public summary guards; no original genetic response or per-species coefficient."""
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'benchmarks/exploratory/host_memory_sampling_overlap_v0.1.json'


def test_public_641_external_predictor_geography_domain_and_consumed_original():
    x=json.loads(REPORT.read_text())
    assert x['status']=='POST_OUTCOME_RESPONSE_BLIND_PREDICTOR_GEOMETRY_AUDIT'
    assert x['source_bound_species']==641
    assert x['response_blind_external_edge_rows']==71520
    assert x['panels']['development']['nearest_neighbor_partial_r2']['species']==320
    assert x['panels']['confirmatory']['nearest_neighbor_partial_r2']['species']==321
    assert sum(v['species'] for v in x['original_graph_edge_count_strata'])==641
    interp=x['interpretation']
    assert interp['original_259_empirical_status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert interp['coefficient_or_p_value_tested'] is False
    assert interp['original_259_genetic_response_not_reopened'] is True
    assert interp['independently_confirmatory_genetic_panel_used'] is False
    assert 'species_name' not in str(x)


def test_r2_and_adjustment_are_distinguished_not_ecological_effects():
    x=json.loads(REPORT.read_text())
    a=x['panels']['all']
    assert a['nearest_neighbor_partial_r2']['median']==pytest.approx(0.03006705485148804,abs=2e-11)
    assert a['nearest_neighbor_partial_r2']['species_above_0_10']==147
    assert a['nearest_neighbor_partial_r2']['species_above_0_25']==49
    assert a['nearest_neighbor_df_adjusted_partial_r2']['median']==pytest.approx(0.007089410985149591,abs=2e-11)
    assert a['joint_sampling_geometry_partial_r2']['median']==pytest.approx(0.10846522741595432,abs=2e-11)
    assert a['joint_sampling_geometry_df_adjusted_partial_r2']['median']==pytest.approx(0.06171870664382151,abs=2e-11)
    b=x['original_graph_edge_count_strata']
    assert [s['species'] for s in b]==[185,300,156]
    assert b[0]['median_nearest_neighbor_df_adjusted_partial_r2'] < 0
    assert b[-1]['median_nearest_neighbor_df_adjusted_partial_r2'] > 0
