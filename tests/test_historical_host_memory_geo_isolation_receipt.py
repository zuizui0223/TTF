"""Public aggregate fixture; no private per-species or per-locality response is loaded."""
from pathlib import Path
import json
import math

PUBLIC=Path(__file__).resolve().parents[1] / 'benchmarks/exploratory/historical_host_memory_259_geo_isolation_by_signal_v0.1.json'


def test_original_population_and_frozen_response_firewall():
    r=json.loads(PUBLIC.read_text())
    assert r['schema']=='ttf_hhm259_fragility_by_semipartial_and_sampling_isolation_v0.1'
    assert r['original_confirmatory_status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert (r['n_frozen_survivors'],r['n_nonzero_defined'],r['n_invariant'])==(259,255,4)
    assert (r['frozen_original_edges'],r['attempted_site_deletions'],r['defined_site_deletions'],r['undefined_deletion_variants'],r['flip_event_total'])==(20362,6374,6373,1,420)
    assert r['new_p_values']==0 and r['new_hypothesis_tests']==0
    assert r['original_panel_changed'] is False and r['genetic_identity_reopened'] is False
    assert 'Abagrotis' not in PUBLIC.read_text()


def test_quartiles_reconcile_exact_prepublished_fragility_counts():
    q=json.loads(PUBLIC.read_text())['semipartial_quartile_geography']
    assert [x['species'] for x in q]==[64,64,63,64]
    assert [x['flip_events'] for x in q]==[303,79,30,8]
    assert [x['species_with_any_reversal'] for x in q]==[56,44,21,8]
    assert [x['defined_deletions'] for x in q]==[1843,1577,1586,1302]
    assert sum(x['flip_events'] for x in q)==420
    assert sum(x['nonflip_defined_events'] for x in q)==6308-420
    assert sum(x['undefined_deletions'] for x in q)==1
    assert all(x['species_with_both_types']==x['species_with_any_reversal'] for x in q)


def test_overall_nearest_rank_reproduces_independent_geo_audit():
    r=json.loads(PUBLIC.read_text())
    assert math.isclose(r['overall_flip_site_nearest_neighbor_rank'],0.5722075791028662,abs_tol=1e-11)
    assert math.isclose(r['overall_other_site_nearest_neighbor_rank'],0.5164893823783056,abs_tol=1e-11)
    assert [round(q['mean_flip_nearest_rank'],3) for q in r['semipartial_quartile_geography']]==[.545,.612,.695,.744]
    assert [round(q['mean_other_nearest_rank'],3) for q in r['semipartial_quartile_geography']]==[.512,.515,.516,.523]
    assert [q['species_with_both_types'] for q in r['semipartial_quartile_geography']]==[56,44,21,8]
