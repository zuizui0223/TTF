from __future__ import annotations

import numpy as np

from ttf.lepidoptera_host_resource_empirical import (
    frozen_pair_surface,
    pair_transfer_congruence,
    score_empirical_host_resource,
)


class Design(dict):
    pass


def toy_design():
    sources=[f"s{i}" for i in range(3)]
    targets=[f"t{i}" for i in range(6)]
    names=tuple(sources+targets)
    target=[]; source=[]; residual=[]; lengths=[]; offsets=[]
    at=[]; ass=[]; cursor=0
    for ti in range(6):
        for si in range(3):
            target.append(3+ti); source.append(si)
            residual.append(float(si-1))
            lengths.append(4); offsets.append(cursor); cursor+=4
            at.extend([0,1,2,3]); ass.extend([0,1,2,3])
    return Design(
        species_order=np.asarray(names,dtype="U"),
        target_index=np.asarray(target,np.int64),
        source_index=np.asarray(source,np.int64),
        host_resource_residual=np.asarray(residual,float),
        alignment_length=np.asarray(lengths,np.int64),
        alignment_offset=np.asarray(offsets,np.int64),
        alignment_target=np.asarray(at,np.int64),
        alignment_source=np.asarray(ass,np.int64),
    )


def toy_response():
    source={
        "s0":np.asarray([-1.,-0.4,.4,1.]),
        "s1":np.asarray([-1.,.2,-.2,1.]),
        "s2":np.asarray([1.,.4,-.4,-1.]),
    }
    target={}
    for i in range(6):
        target[f"t{i}"]=np.asarray([1.,.5,-.5,-1.])+i*.01
    return {**source,**target}


def test_frozen_pair_surface_and_congruence_are_finite():
    z=toy_design()
    names,target,source,residual,pair_id=frozen_pair_surface(z)
    score=pair_transfer_congruence(
        toy_response(),names,target,source,pair_id,
        z["alignment_target"],z["alignment_source"],
    )
    assert len(score)==18
    assert np.isfinite(score).all()
    assert residual.shape==(18,)


def test_empirical_host_resource_uses_three_nuisance_references():
    z=toy_design()
    reference=np.linspace(-1.0,-0.5,499)
    scored=score_empirical_host_resource(
        toy_response(),z,
        private_reference=reference,
        geometry_reference=reference,
        breadth_reference=reference,
        alpha=.05,
    )
    assert np.isfinite(scored.statistic)
    assert 0.0 < scored.p_value <= 1.0
    assert len(scored.per_target)==6
