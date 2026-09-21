from __future__ import annotations

import numpy as np

from ttf.lepidoptera_pair_repeatability import (
    _two_way_demean,
    prepare_pair_repeatability,
    score_pair_repeatability,
)
from ttf.lepidoptera_pair_repeatability_simulate import (
    EVALUATION_TAG,
    REFERENCE_TAG,
    frozen_seed,
    relational_latent_kernel,
)


def toy_design():
    names=tuple([f"s{i}" for i in range(5)]+[f"t{i}" for i in range(6)])
    target=[]; source=[]; offsets=[]; lengths=[]; at=[]; ass=[]
    nuisance=[]; raw=[]; cursor=0
    for ti in range(5,11):
        for si in range(5):
            target.append(ti); source.append(si)
            offsets.append(cursor); lengths.append(20)
            at.extend(range(20)); ass.extend(range(20))
            nuisance.append([ti,si,ti-si,1.0,0.25,1.1,2.0,0.5])
            raw.append((si+1)/(ti+1))
            cursor+=20
    return {
        "species_order":np.asarray(names,dtype="U"),
        "target_index":np.asarray(target,np.int64),
        "source_index":np.asarray(source,np.int64),
        "nuisance_covariates":np.asarray(nuisance,float),
        "raw_host_resource_jaccard":np.asarray(raw,float),
        "alignment_offset":np.asarray(offsets,np.int64),
        "alignment_length":np.asarray(lengths,np.int64),
        "alignment_target":np.asarray(at,np.int64),
        "alignment_source":np.asarray(ass,np.int64),
    }


def test_two_way_demean_removes_group_means():
    source=np.tile(np.arange(5),6)
    target=np.repeat(np.arange(6),5)
    y=2.0*source-3.0*target+np.sin(np.arange(30))
    r=_two_way_demean(y,source,target,5,6)
    for label in range(5):
        assert abs(float(np.mean(r[source==label])))<1e-10
    for label in range(6):
        assert abs(float(np.mean(r[target==label])))<1e-10


def test_pair_split_is_deterministic_and_balanced():
    d=toy_design()
    a=prepare_pair_repeatability(d,np.asarray([32]*11),minimum_rows_per_half=8)
    b=prepare_pair_repeatability(d,np.asarray([32]*11),minimum_rows_per_half=8)
    assert np.array_equal(a.half_a_global_target,b.half_a_global_target)
    assert np.array_equal(a.half_b_global_source,b.half_b_global_source)
    assert np.all(a.half_a_rows==10)
    assert np.all(a.half_b_rows==10)
    assert len(a.pair_index)==30


def test_repeatability_score_is_finite():
    d=toy_design()
    p=prepare_pair_repeatability(d,np.asarray([32]*11),minimum_rows_per_half=8)
    x=np.linspace(-1,1,32)
    response={}
    for i,name in enumerate(d["species_order"]):
        response[str(name)]=np.sin(x*(i+1))+0.05*i*x
    s=score_pair_repeatability(response,tuple(map(str,d["species_order"])),p)
    assert np.isfinite(s.statistic)
    assert np.isfinite(s.covariance)
    assert s.retained_pairs==30


def test_relational_kernel_is_psd_and_seed_namespaces_are_disjoint():
    names=tuple(f"sp{i}" for i in range(12))
    k=relational_latent_kernel(names,dimensions=8)
    assert np.allclose(k,k.T)
    assert np.allclose(np.diag(k),1.0)
    assert float(np.min(np.linalg.eigvalsh(k)))>-1e-10
    seeds={
        frozen_seed(20260921,tag,cell,0)
        for tag in (REFERENCE_TAG,EVALUATION_TAG)
        for cell in (
            "private","geometry_confounded_trap","host_breadth_confounded_trap",
            "host_resource_gradient_trap","relational_latent_positive"
        )
    }
    assert len(seeds)==10
