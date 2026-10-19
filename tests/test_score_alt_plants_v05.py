import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from alternative_plant_match_v05 import (equal_reference_indices,observation_signature,genus_of)
from same_niche_null_v04 import score_real_and_placebo

def test_identical_modern_and_historical_alt_plant_climate_gives_no_unique_effect():
    rng=np.random.default_rng(9)
    current=rng.normal(size=(75,5,4))
    past=current+rng.normal(size=(75,5,4),scale=.27)
    plant=rng.normal(size=(60,4));insect=rng.normal(size=(60,4))
    x=score_real_and_placebo(current,past,plant,insect,plant,np.log1p(np.arange(75)))
    assert x[0]['unique_fraction']==pytest.approx(x[1]['unique_fraction'],abs=1e-12)
    assert x[0]['residual_variance']==pytest.approx(x[1]['residual_variance'],abs=1e-12)
    assert x[0]['host_current_mean']==pytest.approx(x[1]['host_current_mean'],abs=1e-12)

def test_equal_reference_draws_do_not_introduce_any_count_advantage():
    for ns in [60,61,75,150,200]:
        si,real,alt=equal_reference_indices(200,165,ns,'Butterfly a',0,0)
        target=min(100,200,165,ns//2)
        assert len(si)==len(real)==len(alt)==target

def test_same_name_and_draw_produce_same_samples():
    a=equal_reference_indices(150,150,150,'Moth a',1,7)
    b=equal_reference_indices(150,150,150,'Moth a',1,7)
    assert all(np.array_equal(x,y) for x,y in zip(a,b))

def test_repeated_donor_ids_cannot_pass_observation_hashing():
    with pytest.raises(ValueError):
        observation_signature([{'gbif_id':123},{'gbif_id':123}])

def test_genus_hides_plant_taxonomic_binomial():
    assert genus_of('Ulmus minor subsp. minor')=='ulmus'
    assert genus_of('Ulmus americana')=='ulmus'
    assert genus_of('Picea abies')=='picea'

def test_full_scoring_runner_import_compiles_before_opening_rasters():
    # Regression: pure helper tests passed while production scoring script had
    # SyntaxError in output-row construction. Import must check its syntax.
    import importlib
    runner=importlib.import_module('score_matched_alternative_plants_v05')
    assert callable(runner.main)
