import numpy as np

from ttf.ecological_information_ladder import ecological_information_ladder


def make_panel(ns=12,nt=10,seed=45):
    rng=np.random.default_rng(seed)
    n_species=ns+nt
    source=np.repeat(np.arange(ns),nt)
    target=np.tile(np.arange(ns,n_species),ns)
    class_name=np.asarray(["A"]*n_species)
    order_name=np.asarray([f"o{i%5}" for i in range(n_species)])
    family_name=np.asarray([f"f{i%8}" for i in range(n_species)])
    occurrence_n=np.asarray([35+i%25 for i in range(n_species)])
    return rng,source,target,class_name,order_name,family_name,occurrence_n


def test_ladder_history_addition_collapses_when_history_is_current():
    rng,s,t,c,o,f,n=make_panel(seed=46)
    current=rng.normal(size=len(s))
    history=current+0.01*rng.normal(size=len(s))
    out=ecological_information_ladder(
        class_name=c,order_name=o,family_name=f,occurrence_n=n,
        source_index=s,target_index=t,r_current=current,r_hist=history,
    )
    assert out["layers"]["B1_current_climate"]["total_unique_variance_fraction"]>0
    assert out["layers"]["C1_history_given_current"]["control_unique_variance_fraction_after_fe"]<0.01
    assert out["cross_layer"]["history_fraction_remaining_after_adding_current_to_baseline"]<0.01


def test_ladder_keeps_independent_historical_information():
    rng,s,t,c,o,f,n=make_panel(seed=47)
    current=rng.normal(size=len(s))
    history=rng.normal(size=len(s))
    out=ecological_information_ladder(
        class_name=c,order_name=o,family_name=f,occurrence_n=n,
        source_index=s,target_index=t,r_current=current,r_hist=history,
    )
    assert out["layers"]["C1_history_given_current"]["control_unique_variance_fraction_after_fe"]>0.7
    assert out["cross_layer"]["history_fraction_remaining_after_adding_current_to_baseline"]>0.7
    assert out["interpretation"]["legacy_B_or_C_reopened"] is False
