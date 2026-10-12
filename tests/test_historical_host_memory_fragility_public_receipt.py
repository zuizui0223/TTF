"""Public aggregate snapshot QA; no private species names, raw COI or edge distances."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
RESULT=HERE / 'benchmarks' / 'exploratory' / 'historical_host_memory_259_fragility_drivers_v0.1.json'


def load():
    return json.loads(RESULT.read_text())


def test_original_stopped_one_shot_stays_the_only_confirmatory_result():
    a=load()
    assert a['original_confirmatory_statistic_status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert a['status']=='EXPLORATORY_AFTER_EMPIRICAL_OUTCOME_EXPOSURE'
    assert a['new_p_values']==0
    assert a['limitations']['original_primary_estimate_still_undefined'] is True
    assert a['limitations']['no_new_confirmatory_panel'] is True
    assert a['limitations']['no_raw_nucleotide_sequences_or_edge_p_distances_read_by_this_script'] is True


def test_source_provenance_and_deletion_accounting():
    a=load()
    assert all(len(v)==64 for v in a['source_sha256'].values())
    p=a['panel']
    assert (p['n_all'],p['n_nonzero_defined'],p['n_baseline_invariant']) == (259,255,4)
    assert (p['n_deletions_attempted'],p['n_deletions_defined'],p['n_deletions_undefined']) == (6374,6373,1)
    assert p['n_defined_deletions_among_nonzero_species']==6308
    assert p['n_with_any_sign_reversal']==129
    assert p['n_reversal_events']==420
    assert abs(p['fraction_nonzero_species_deletions_reversing']-420/6308)<1e-15


def test_all_quartiles_and_host_counts_are_public_aggregates():
    a=load()
    b=a['absolute_standardized_beta_quartiles']
    s=a['absolute_semipartial_quartiles']
    assert [q['n'] for q in b] == [64,64,63,64]
    assert [q['n_with_any_reversal'] for q in b] == [57,40,19,13]
    assert [q['n_with_any_reversal'] for q in s] == [56,44,21,8]
    assert sum(q['reversal_events'] for q in s)==420
    assert sum(q['defined_locality_deletions'] for q in s)==6308
    assert [(q['number_of_frozen_recorded_host_species'],q['total_species']) for q in a['recorded_host_metadata_by_count']] == [(1,143),(2,70),(3,46)]
    assert a['identity_verification']['r_semi_equals_beta_times_sqrt_unique_fraction_max_absolute_error']<1e-12
    content=RESULT.read_text()
    assert 'Abagrotis nefascia' not in content
    assert 'Holomelina laeta' not in content
    assert 'locality_index' not in content


def test_family_held_out_prediction_is_exploratory_and_not_mechanism_proof():
    a=load()
    m=a['family_held_out_internal_prediction']
    assert m['cv']=='five_fold_family_held_out_GroupKFold'
    scores=m['models']
    beta=scores['absolute_standardized_beta_and_locality_count']['brier_score']
    semi=scores['absolute_semipartial_and_locality_count']['brier_score']
    host=scores['semipartial_plus_recorded_host_count']['brier_score']
    assert 0.139<beta<0.142
    assert 0.120<semi<0.123
    assert semi<beta
    assert host>=semi
    assert m['output_not_a_genetic_effect_p_value'] is True
    assert a['limitations']['CV_score_is_not_out_of_sample_validation_of_an_ecological_host_memory_mechanism'] is True


def test_original_259_confirmatory_terminal_receipt_is_byte_pinned_and_undefined():
    import hashlib
    report=load()
    original=HERE/'benchmarks'/'frozen'/'historical_host_memory_259_one_shot_empirical_result_v0.1.json'
    digest=hashlib.sha256(original.read_bytes()).hexdigest()
    assert digest==report['source_sha256']['primary_result']
    original_data=json.loads(original.read_text())
    assert original_data['status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert original_data['n_species_frozen']==259
    assert original_data['undefined_species']==4
    assert original_data['primary_statistic'] is None
    assert original_data['monte_carlo_p'] is None
