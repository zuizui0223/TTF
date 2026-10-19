"""Regression tests for paired species / genus audit logic; not genetic validation."""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_lgm_host_identity_v05 import classify_transition, audit_pair_records, genus_cluster_interval

def test_transition_is_defined_on_same_species_only():
    prior={'A alpha':0.02,'B beta':-0.01,'C gamma':0.03}
    current={'A alpha':-0.004,'B beta':0.002,'C gamma':0.001}
    trans=classify_transition(prior,current)
    assert trans=={'positive_to_positive':1,'positive_to_negative':1,
                   'negative_to_positive':1,'negative_to_negative':0}
    with pytest.raises(ValueError):
        classify_transition(prior,{'A alpha':0.1})

def test_real_alt_pair_counts_not_comparison_of_global_medians():
    rows=[{'species':'A alpha','fraction_diff_of_medians':-.1,
           'raw_residual_paired_diff_median':-.2,'actual_fraction_median':.1,
           'alternate_fraction_median':.2,'current_mismatch_gap':-.01},
          {'species':'B beta','fraction_diff_of_medians':.01,
           'raw_residual_paired_diff_median':.1,'actual_fraction_median':.01,
           'alternate_fraction_median':0,'current_mismatch_gap':.00},
          {'species':'C gamma','fraction_diff_of_medians':.02,
           'raw_residual_paired_diff_median':.2,'actual_fraction_median':.01,
           'alternate_fraction_median':0,'current_mismatch_gap':.03}]
    q=audit_pair_records(rows)
    assert q['positive_species_fraction_count']==2
    assert q['positive_species_raw_variance_count']==2
    assert q['species_median_fraction_difference']==pytest.approx(.01)
    assert q['median_actual_fraction']==pytest.approx(.01)
    assert q['median_alt_fraction']==pytest.approx(0.)

def test_genus_cluster_interval_is_seed_deterministic():
    rows=[{'species':'A a','fraction_diff_of_medians':.1},
          {'species':'A b','fraction_diff_of_medians':.2},
          {'species':'B c','fraction_diff_of_medians':-.05},
          {'species':'C d','fraction_diff_of_medians':-.1}]
    x=genus_cluster_interval(rows,draws=40,seed=123)
    y=genus_cluster_interval(rows,draws=40,seed=123)
    assert x==y and x['genus_clusters']==3
    assert -1<=x['median_diff_limits'][0]<=x['median_diff_limits'][1]<=1

def test_v05_source_has_exactly_333_scored_species_and_eight_groups():
    from audit_lgm_host_identity_v05 import audit_sources
    root=Path(__file__).resolve().parents[1]/'benchmarks'/'exploratory'
    data=audit_sources(root/'lgm_same_niche_split_null_v04',
                       root/'lgm_modern_matched_alternative_plant_v05')
    assert data['species_compared']==333
    assert data['group_rows']==2664
    assert data['v04_actual_positive_on_shared_species']==269
    assert data['v05_actual_positive_on_shared_species']==147
    assert data['v05_raw_positive_on_shared_species']==142
    assert data['transition_counts']['positive_to_negative']==132
    assert data['transition_counts']['negative_to_positive']==10
