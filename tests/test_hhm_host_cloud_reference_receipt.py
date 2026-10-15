"""Verify published host-cloud rarefaction aggregates without genetic responses."""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'benchmarks/exploratory/historical_host_lcp_reference_size_26of665_result_v0.1.json'
CONTRACT=ROOT/'docs/supporting/historical_host_lcp_reference_size_rarefaction_contract_v0.1.json'
ORIGINAL=ROOT/'benchmarks/frozen/historical_host_memory_259_one_shot_empirical_result_v0.1.json'


def load():
    return json.loads(RECEIPT.read_text())


def test_original_64_edge_replay_and_sha_bound_exploratory_scope():
    o=load()
    q=json.loads(CONTRACT.read_text())
    assert o['contract_sha256']==hashlib.sha256(CONTRACT.read_bytes()).hexdigest()
    assert q['status']=='POST_LCP_OUTCOME_EXPLORATORY_FIXED_BEFORE_FIRST_HOST_RAREFACTION_RESULT'
    assert o['status']=='EXPLORATORY_HOST_REFERENCE_SIZE_SENSITIVITY_ONLY'
    assert o['species']=='Bungalotis erythus' and o['sites']==26 and o['edges']==64
    assert o['original_reference_cloud_sizes']=={'host':665,'self':26,'global_climate_standardization':6343}
    assert o['same_pixel_edges_intact']==1
    assert o['maximum_original_LCP_replay_abs_error']<=1e-9
    assert o['source_seed_namespace']=='ttf-hhm-host-reference-26-of-665-50draws-v0.1'
    assert o['original_pilot_sha256']['edges']=='7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b'


def test_fifty_draw_means_quantiles_and_climate_difference_recalculate():
    o=load();x=np.asarray(o['rarefaction']['host_delta_mean_by_draw'],float)
    assert len(x)==50
    assert o['rarefaction']['draw_count']==50 and o['rarefaction']['sample_size']==26
    assert x.size==50 and np.isfinite(x).all()
    assert np.mean(x)==o['rarefaction']['host_delta_mean_across_draws']
    assert np.median(x)==o['rarefaction']['host_delta_median_across_draws']
    assert np.min(x)==o['rarefaction']['host_delta_min_across_draws']
    assert np.max(x)==o['rarefaction']['host_delta_max_across_draws']
    assert np.quantile(x,0.1)==o['rarefaction']['host_delta_quantile_10']
    assert np.quantile(x,0.9)==o['rarefaction']['host_delta_quantile_90']
    self_d=o['original_self_delta_cost_km_mean']
    host_d=o['original_full_host_delta_cost_km_mean']
    assert abs(host_d-1.665783407966145)<1e-12
    assert abs(self_d-17.794186622317866)<1e-12
    assert abs((host_d-self_d)-o['original_full_host_minus_self_delta_cost_km_mean'])<1e-12
    dif=x-self_d
    for key,value in [('mean_across_draws',np.mean(dif)),('median_across_draws',np.median(dif)),('min_across_draws',np.min(dif)),('max_across_draws',np.max(dif)),('fraction_positive_draws_descriptive_only',np.mean(dif>0))]:
        assert abs(o['rarefied_host_minus_fixed_self'][key]-value)<1e-12
    assert int(np.sum(x>0))==48
    assert int(np.sum(dif>0))==0


def test_old_terminal_coi_result_unmodified_and_no_replacement_test():
    o=load();terminal=json.loads(ORIGINAL.read_text())
    assert terminal['status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert terminal['monte_carlo_p'] is None
    assert o['scientific_boundary']['original_259_result']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert o['scientific_boundary']['no_genetic_response_or_coi_distance_reopened'] is True
    assert o['scientific_boundary']['rarefaction_quantiles_are_not_CI'] is True
    assert o['rarefaction']['new_empirical_genetic_test'] is False
    assert o['rarefaction']['new_p_values']==0
