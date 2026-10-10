"""Mask-survivor qualification tests use only artificial geography and predictors."""
import itertools

import numpy as np
import pytest

from scripts.qualify_historical_host_memory_synthetic import prepare_species
from scripts.qualify_historical_host_memory_survivor_information import EXPECTED_MASK_BINDINGS

from scripts.qualify_historical_host_memory_survivor_synthetic import (
    GROUPS, survivor_beta_worlds, survivor_seed, validate_source_set,
)


def make_toy():
    rng = np.random.default_rng(20261010)
    xyz = rng.normal(size=(7,3))*150.0
    loc = {i:p for i,p in enumerate(xyz)}
    graph = dict(enumerate(itertools.combinations(range(7), 2)))
    predictors = [
        {"edge_index":str(i), "M_host":str(rng.normal()), "M_self":str(rng.normal()),
         "D_host_0BP":str(rng.normal())}
        for i in range(len(graph))
    ]
    pre=prepare_species("Synthetic unrelated species", predictors, loc, graph)
    pre["xyz"]=xyz
    return pre


def test_survivor_seeds_disjoint_across_world_types_and_pre_mask():
    from scripts.qualify_historical_host_memory_synthetic import WORLD_SEEDS, seed
    name="Synthetic unrelated species"
    values=[survivor_seed(group,name) for group in ("reference","evaluation","positive","spatial-basis")]
    assert len(set(values))==4
    assert all(survivor_seed(group,name)!=seed(WORLD_SEEDS[group],name)
               for group in ("reference","evaluation","positive"))


def test_survivor_simulation_reproducible_and_has_no_real_genetic_response():
    toy=make_toy()
    one=survivor_beta_worlds(toy,"evaluation",500)
    two=survivor_beta_worlds(toy,"evaluation",500)
    pos=survivor_beta_worlds(toy,"positive",500)
    np.testing.assert_array_equal(one,two)
    assert np.all(np.isfinite(one)) and np.all(np.isfinite(pos))
    assert not np.array_equal(one,pos)


def test_frozen_world_count_and_group_are_not_tunable():
    toy=make_toy()
    with pytest.raises(RuntimeError,match="frozen survivor synthetic"):
        survivor_beta_worlds(toy,"evaluation",499)
    with pytest.raises(RuntimeError,match="frozen survivor synthetic"):
        survivor_beta_worlds(toy,"rescue",500)


def make_roles(n=231):
    original=[{"species":f"Con{i:03d}","panel":"confirmatory","split_key_sha256":f"k{i}"}
              for i in range(321)]
    original += [{"species":f"Dev{i:03d}","panel":"development","split_key_sha256":f"d{i}"}
                 for i in range(320)]
    # Fixture must preserve the frozen original role set independently of a mutated survivor input.
    survivors=[dict(r) for r in original[:n]]
    mask={"schema":"ttf_historical_host_memory_confirmatory_mask_result_v0.1",
          "source_bindings":{**EXPECTED_MASK_BINDINGS,"mask_rule_sha256":"0"*64},
          "status":("PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION"
                    if n>=200 else "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"),
          "survivor_names":[r["species"] for r in survivors]}
    summary={"schema":"ttf_historical_host_memory_survivor_information_v0.1",
             "decision":"PASS_TO_EXACT_SURVIVOR_SYNTHETIC_REQUALIFICATION",
             "genetic_response_opened":False,"nucleotide_identity_opened":False,
             "post_IBD_turnover_opened":False,"empirical_beta_host_opened":False,
             "survivor_synthetic_qualification_passed":False,
             "mask_result_sha256":"fixture",
             "surviving_species_names":[r["species"] for r in survivors],
             "surviving_frozen_edges":n}
    preds=[{"species":r["species"],"edge_index":"0","M_host":"1",
            "M_self":"0","D_host_0BP":"0"} for r in survivors]
    return summary,survivors,preds,original,mask


def test_exact_survivor_subset_and_roles_authorized():
    summary,roles,predictors,original,mask=make_roles()
    names,by_sp=validate_source_set(summary,roles,predictors,original,mask)
    assert len(names)==231 and set(by_sp)==set(names)


def test_survivor_resplitting_or_species_backfill_fails_closed():
    summary,roles,pred,orig,mask=make_roles()
    roles[0]["panel"]="development"
    with pytest.raises(RuntimeError,match="resplitting"):
        validate_source_set(summary,roles,pred,orig,mask)
    summary,roles,pred,orig,mask=make_roles()
    mask["survivor_names"].append("Dev000")
    with pytest.raises(RuntimeError,match="mask subset"):
        validate_source_set(summary,roles,pred,orig,mask)


def test_survivor_below_floor_and_unqualified_information_stop():
    summary,roles,pred,orig,mask=make_roles(199)
    with pytest.raises(RuntimeError,match="character-mask decision is not PASS"):
        validate_source_set(summary,roles,pred,orig,mask)
    # Even a forged PASS status cannot bypass the predeclared 200-species floor.
    mask["status"]="PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION"
    with pytest.raises(RuntimeError,match="200-species"):
        validate_source_set(summary,roles,pred,orig,mask)
    summary,roles,pred,orig,mask=make_roles()
    summary["decision"]="NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_SURVIVOR_INFORMATION"
    with pytest.raises(RuntimeError,match="did not qualify"):
        validate_source_set(summary,roles,pred,orig,mask)


def test_survivor_simulation_rejects_forged_mask_provenance():
    summary, roles, predictors, original, mask = make_roles()
    mask["source_bindings"]["edge_geometry_sha256"] = "f" * 64
    with pytest.raises(RuntimeError, match="provenance binding mismatch"):
        validate_source_set(summary, roles, predictors, original, mask)
