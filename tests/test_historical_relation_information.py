import numpy as np

from ttf.historical_relation_information import historical_relation_information


def test_history_information_finds_present_convergence_with_distinct_histories():
    ns=12
    nt=10
    species=np.asarray([f"sp{i}" for i in range(ns+nt)])
    source=np.repeat(np.arange(ns),nt)
    target=np.tile(np.arange(ns,ns+nt),ns)
    rng=np.random.default_rng(31)

    r_current=np.clip(0.75+0.08*rng.normal(size=len(source)),0,1)
    latent_source=rng.normal(size=ns)
    latent_target=rng.normal(size=nt)
    r_hist=1/(1+np.exp(-(latent_source[source]-latent_target[target-ns]+0.5*rng.normal(size=len(source)))))

    classes=np.asarray(["Insecta"]*(ns+nt))
    orders=np.asarray([f"o{i%5}" for i in range(ns+nt)])
    families=np.asarray([f"f{i%8}" for i in range(ns+nt)])
    occurrence_n=np.asarray([40+(i%20) for i in range(ns+nt)])

    out=historical_relation_information(
        species_order=species,
        class_name=classes,
        order_name=orders,
        family_name=families,
        occurrence_n=occurrence_n,
        source_index=source,
        target_index=target,
        r_hist=r_hist,
        r_current=r_current,
    )
    eco=out["ecological_geometry"]
    assert eco["high_current_dyads"]>0
    assert eco["high_current_low_history_dyads"]>0
    assert 0<=eco["high_current_low_history_fraction"]<=1
    assert 0<=out["nonredundancy"]["total_unique_variance_fraction"]<=1
    assert out["dyadic_signal_support"]["signal_effective_sources"]>1
    assert out["dyadic_signal_support"]["signal_effective_targets"]>1
    assert out["excluded_at_this_stage"]==["geographic_opportunity","genetic_response"]


def test_history_information_redundancy_collapses_when_hist_tracks_current():
    ns=10
    nt=9
    species=np.asarray([f"sp{i}" for i in range(ns+nt)])
    source=np.repeat(np.arange(ns),nt)
    target=np.tile(np.arange(ns,ns+nt),ns)
    rng=np.random.default_rng(32)
    current=rng.normal(size=len(source))
    hist=current+0.01*rng.normal(size=len(source))

    out=historical_relation_information(
        species_order=species,
        class_name=np.asarray(["A"]*(ns+nt)),
        order_name=np.asarray([f"o{i%4}" for i in range(ns+nt)]),
        family_name=np.asarray([f"f{i%6}" for i in range(ns+nt)]),
        occurrence_n=np.asarray([30+i for i in range(ns+nt)]),
        source_index=source,
        target_index=target,
        r_hist=hist,
        r_current=current,
    )
    assert out["nonredundancy"]["control_unique_variance_fraction_after_fe"]<0.01
