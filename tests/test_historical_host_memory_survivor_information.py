import pytest

from scripts.qualify_historical_host_memory_survivor_information import qualify_survivors, EXPECTED_MASK_BINDINGS


def fixture(n_survivors=250, unique=0.25):
    confirm=[f"C{i:03d}" for i in range(321)]
    roles=[{"species":sp,"panel":"confirmatory"} for sp in confirm]
    roles += [{"species":f"D{i:03d}","panel":"development"} for i in range(320)]
    d=[{"species":r["species"],"status":"complete",
        "unique_fraction_M_host":str(unique),"predictor_condition_number":"2"} for r in roles]
    d.append({"species":"Excluded","status":"nonfinite_edge_climate",
              "unique_fraction_M_host":"","predictor_condition_number":""})
    ledger=[{"species":sp,"status":"PASS_MASK" if i<n_survivors else "NOT_EVALUABLE_CHARACTER_SUPPORT",
             "survives":i<n_survivors,"invalid_edges":0 if i<n_survivors else 1,
             "valid_edges":10 if i<n_survivors else 9,"frozen_edges":10}
             for i,sp in enumerate(confirm)]
    mask={"schema":"ttf_historical_host_memory_confirmatory_mask_result_v0.1",
          "status":("PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION"
                    if n_survivors>=200 else
                    "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"),
          "source_archive_sha256":"5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5",
          "source_bindings":{**EXPECTED_MASK_BINDINGS,"mask_rule_sha256":"0"*64},
          "confirmatory_before_mask":321,
          "survivor_names":confirm[:n_survivors],
          "surviving_species":n_survivors,
          "ledger":ledger,
          "nucleotide_identity_persisted":False,"genetic_distances_opened":False,
          "empirical_beta_host_opened":False,"no_graph_repair":True,"no_resplitting":True,
          "no_backfill":True,"direct_genetic_opening_authorized":False}
    rule={"post_mask_requalification":{"minimum_confirmatory_survivor_species":200},
          "predictor_information":{
              "minimum_species_fraction_with_unique_fraction_at_least_0.05":0.7,
              "minimum_panel_median_unique_fraction":0.1,
              "maximum_standardized_predictor_condition_number":30,
              "minimum_species_fraction_passing_condition_number":0.9}}
    return mask,roles,d,rule


def test_unchanged_post_mask_information_thresholds():
    r,selected=qualify_survivors(*fixture())
    assert len(selected)==250
    assert r["decision"]=="PASS_TO_EXACT_SURVIVOR_SYNTHETIC_REQUALIFICATION"
    assert r["survivor_synthetic_qualification_passed"] is False
    assert r["genetic_response_opened"] is False


def test_below_200_is_terminal_stopping_not_species_backfill():
    result,selected=qualify_survivors(*fixture(199))
    assert len(selected)==199
    assert result["decision"]=="NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"


def test_low_host_unique_information_stops_without_retune():
    result,_=qualify_survivors(*fixture(250,unique=0.04))
    assert result["decision"]=="NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_SURVIVOR_INFORMATION"
    assert result["checks"]["fraction_nonredundant_ge_0_70"] is False


def test_reject_mask_ledger_drift_and_genetic_unblinding():
    mask,roles,diag,rule=fixture()
    mask["survivor_names"].append("D000")
    with pytest.raises(RuntimeError,match="mask survivor identities"):
        qualify_survivors(mask,roles,diag,rule)
    mask,roles,diag,rule=fixture()
    mask["nucleotide_identity_persisted"]=True
    with pytest.raises(RuntimeError,match="genetic opening"):
        qualify_survivors(mask,roles,diag,rule)


def test_reject_stealth_graph_repair_and_rewiring():
    mask,roles,diag,rule=fixture()
    mask["ledger"][0]["invalid_edges"]=1
    with pytest.raises(RuntimeError,match="invalid frozen edge"):
        qualify_survivors(mask,roles,diag,rule)


def test_mask_status_must_match_exact_survivor_count():
    mask, roles, diagnostics, rule = fixture(250)
    mask["status"] = "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"
    with pytest.raises(RuntimeError, match="mask decision/count mismatch"):
        qualify_survivors(mask, roles, diagnostics, rule)


def test_exact_mask_provenance_cannot_be_substituted():
    mask, roles, diagnostic, rule = fixture()
    mask["source_bindings"]["panel_roles_sha256"] = "0" * 64
    with pytest.raises(RuntimeError, match="provenance bindings"):
        qualify_survivors(mask, roles, diagnostic, rule)
    mask, roles, diagnostic, rule = fixture()
    del mask["source_bindings"]["mask_rule_git_blob_sha1"]
    with pytest.raises(RuntimeError, match="provenance bindings"):
        qualify_survivors(mask, roles, diagnostic, rule)
