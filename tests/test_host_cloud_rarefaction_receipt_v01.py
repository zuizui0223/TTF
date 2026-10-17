"""Public, no-genetics regression audit for the exact CHELSA host cloud sensitivity."""
from pathlib import Path
import hashlib
import importlib.util
import json

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'benchmarks/exploratory/historical_host_memory_bungalotis_26point_host_rarefaction_v0.1.json'
ORIGINAL_REPORT_SHA='c12d2dbcf4601c42ce5cfac5c887c80d004417d0705fbd65bb5a9c19f3621a6b'


def report():
    return json.loads(REPORT.read_text())


def test_exact_public_receipt_bytes_match_original_successful_actions_artifact():
    assert hashlib.sha256(REPORT.read_bytes()).hexdigest()==ORIGINAL_REPORT_SHA
    j=report()
    assert j['status']=='EXPLORATORY_CLIMATE_ONLY_REFERENCE_SIZE_SENSITIVITY_NOT_GENETIC_TEST'
    assert j['contract_sha256']=='8f9ce1fa83a5a0ead8ab68a7278bd6b6e4d6695241fd515bdfd31bc062da3c72'
    assert j['frozen_replay_max_absolute_lcp_error']==0.0
    assert j['study']['frozen_edges']==64
    assert j['study']['localities']==26
    assert j['study']['draws']==64
    assert j['study']['original_host_reference_n']==665
    assert j['study']['fixed_insect_self_reference_n']==26
    assert j['input_sha256']['original_pilot_edges']=='7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b'


def test_all_fixed_draws_and_quantiles_recompute_independently():
    j=report()['rarefaction']
    host=np.asarray(j['rarefied_host_delta_cost_km']['means_per_draw'])
    diff=np.asarray(j['rarefied_host_minus_self_delta_cost_km']['means_per_draw'])
    selfmean=j['fixed_self_mean_delta_cost_km']
    assert len(host)==len(diff)==64
    np.testing.assert_allclose(host-selfmean,diff,rtol=0,atol=1e-13)
    for name,values in [('rarefied_host_delta_cost_km',host),('rarefied_host_minus_self_delta_cost_km',diff)]:
        q=j[name]
        for key,val in [('mean',np.mean(values)),('median',np.median(values)),('min',np.min(values)),('max',np.max(values)),
                        ('q05',np.quantile(values,.05)),('q25',np.quantile(values,.25)),
                        ('q75',np.quantile(values,.75)),('q95',np.quantile(values,.95))]:
            assert abs(q[key]-val)<1e-12
    assert np.count_nonzero(host>0)==63
    assert np.count_nonzero(host<0)==1
    assert np.count_nonzero(diff>0)==0
    assert j['rarefied_host_minus_self_delta_cost_km']['positive_draws']==0
    assert np.mean(host)>1.665783407966145
    assert np.max(host)<selfmean


def test_seed_reflects_64_predeclared_host_cloud_subsamples_not_response():
    path=ROOT/'scripts/pilot_host_cloud_size_rarefaction.py'
    spec=importlib.util.spec_from_file_location('rarefaction',path)
    impl=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(impl)
    ix=impl.fixed_subsample_indices()
    expected=hashlib.sha256(ix.astype('<i8').tobytes()).hexdigest()
    assert report()['source_row_order_bounded_sample_indices_sha256']==expected


def test_result_is_explicitly_descriptive_without_p_value_or_land_reconstruction():
    j=report()
    boundary=j['boundary']
    assert boundary['source_panel_already_genetically_opened'] is True
    assert boundary['previous_confirmatory_genetic_status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert boundary['new_genetic_sequences_opened'] is False
    assert boundary['gene_distance_or_p_value_computed'] is False
    assert boundary['rarefaction_quantiles_not_confidence_intervals'] is True
    assert boundary['same_number_of_reference_points_not_same_biological_niche_domain'] is True
    assert boundary['no_verified_21ka_land_or_ice_mask'] is True
    assert boundary['not_a_new_confirmatory_hypothesis'] is True
