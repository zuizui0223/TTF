"""Receipt-only reporting: never opens sequences, fits models or changes inference."""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import tempfile


PREFIX = "ttf_genetic_phylogatr_"
OPEN_FLAGS = (
    "confirmatory_sequence_identity_opened",
    "confirmatory_pairwise_genetic_distances_opened",
    "confirmatory_ttf_statistic_opened",
)
DECISIONS = {
    "TRANSFERABLE_PLACE_COMPONENT": "The qualified test detected out-of-species transfer of post-IBD genetic differentiation within the frozen panel and model family. This supports transferable geographic information, not a causal separation of place from shared demographic history.",
    "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN": "The qualified cross-species test was non-significant, while the separately qualified within-species self test was significant. This is consistent with lineage-conditioned spatial structure in the tested domain; it does not establish zero transfer or a lineage-specific historical cause.",
    "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING": "The qualified cross-species test was non-significant. Qualified empirical self support was insufficient for a lineage-conditioning interpretation. This is not evidence that spatial genetic structure is absent.",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load(path: Path, schema: str | None = None) -> tuple[dict, str]:
    raw = Path(path).read_bytes()
    def reject_constant(value):
        raise ValueError(f"non-finite JSON constant: {value}")
    payload = json.loads(raw, object_pairs_hook=_unique_object, parse_constant=reject_constant)
    _require(isinstance(payload, dict), "receipt must be a JSON object")
    if schema:
        _require(payload.get("schema") == schema, f"unexpected receipt schema: {path}")
    return payload, hashlib.sha256(raw).hexdigest()


def _number(value, name: str, low: float = -1, high: float = 1) -> float:
    _require(type(value) in (int, float), f"{name} must be numeric")
    value = float(value)
    _require(math.isfinite(value) and low <= value <= high, f"{name} out of range")
    return value


def _scores(payload: dict, names: set[str], label: str, *, nullable=False) -> dict:
    _require(isinstance(payload, dict) and set(payload) == names, f"{label} species mismatch")
    return {name: None if value is None and nullable else _number(value, label)
            for name, value in payload.items()}


def _mean_matches(value, scores: dict, label: str) -> None:
    if any(x is None for x in scores.values()):
        _require(value is None, f"{label} must be null when a species is unavailable")
    else:
        actual = _number(value, label)
        _require(math.isclose(actual, math.fsum(scores.values()) / len(scores), abs_tol=1e-12, rel_tol=1e-12),
                 f"{label} does not equal the full species mean")


@dataclass(frozen=True)
class ManuscriptExport:
    markdown: str
    species_csv: str
    manifest: dict


def build_empirical_export(*, result: Path, authorization: Path, qualification: Path,
                           self_qualification: Path, phase4_rule: Path, self_rule: Path) -> ManuscriptExport:
    """Check terminal-receipt consistency and render only validated scalar summaries.

    This is not a rerun of qualification or a provenance authenticity/signature
    service. It checks hashes and consistency of supplied, already-frozen receipts.
    """
    specs = {
        "result": (result, PREFIX + "phase4_empirical_result_v0.1"),
        "authorization": (authorization, PREFIX + "phase4_identity_opening_authorization_v0.1"),
        "qualification": (qualification, PREFIX + "phase3_qualification_v0.1"),
        "self_qualification": (self_qualification, PREFIX + "phase3_self_qualification_v0.1"),
        "phase4_rule": (phase4_rule, PREFIX + "phase4_response_rule_v0.1"),
        "self_rule": (self_rule, PREFIX + "phase3_self_detectability_rule_v0.1"),
    }
    loaded = {name: _load(path, schema) for name, (path, schema) in specs.items()}
    r, a, q, sq, rule, sr = [loaded[name][0] for name in specs]
    hashes = {name: value[1] for name, value in loaded.items()}
    _require(r.get("phase4_authorization_sha256") == hashes["authorization"], "result/authorization hash mismatch")
    for name, key in {
        "qualification": "phase3_qualification_sha256",
        "self_qualification": "phase3_self_qualification_sha256",
        "phase4_rule": "phase4_rule_sha256",
        "self_rule": "phase3_self_rule_sha256",
    }.items():
        _require(a.get(key) == hashes[name], f"authorization/{name} hash mismatch")
    _require(a.get("status") == "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING", "identity opening not authorized")
    _require(r.get("status") == "EMPIRICAL_RESULT_OPENED_UNDER_FROZEN_PHASE4_AUTHORIZATION", "empirical result not completed")
    for key in OPEN_FLAGS:
        _require(r.get(key) is True, f"result is unopened: {key}")
        _require(a.get("outcome_firewall_before_execution", {}).get(key) is False, "authorization firewall drift")
        _require(q.get(key) is False, "qualification contains empirical outcomes")
    for key in OPEN_FLAGS[:2]:
        _require(sq.get(key) is False, "self qualification contains empirical outcomes")
    _require(a.get("outcome_firewall_before_execution", {}).get("decker_empirical_genetic_outcomes_opened") is False,
             "authorization Decker firewall drift")
    for key in ("decker_empirical_genetic_outcomes_opened", "serialized_sequence_identity", "serialized_edge_genetic_distance_vectors"):
        _require(r.get(key) is False, f"result scope drift: {key}")
    fingerprint = a.get("geometry_fingerprint_sha256")
    _require(isinstance(fingerprint, str) and len(fingerprint) == 64, "missing geometry fingerprint")
    for payload in (r, q, sq):
        _require(payload.get("geometry_fingerprint_sha256") == fingerprint, "receipt geometry mismatch")
    _require(q.get("status") == "PASS" and q.get("passed") is True and q.get("phase4_identity_opening_eligible") is True, "Phase-3 qualification did not PASS")
    gates = a.get("phase3_gate_d", {})
    _require(gates.get("status") == "PASS" and all(gates.get(key) is True for key in ("passed", "type1_gate_pass", "power_gate_pass")), "authorization gate mismatch")
    type1, power = q["type1_gate"], q["power_gate"]
    upper = _number(type1["max_observed_wilson95_upper"], "Type-I upper", 0, 1)
    lower = _number(power["observed_wilson95_lower"], "power lower", 0, 1)
    _require(type1.get("pass") is True and type1.get("wilson95_upper_ceiling") == 0.10 and upper <= 0.10, "Type-I gate inconsistent")
    _require(power.get("pass") is True and power.get("wilson95_lower_floor") == 0.80 and lower >= 0.80, "power gate inconsistent")
    split = a["species"]
    training, evaluation = split["train_species"], split["eval_species"]
    _require(isinstance(training, list) and isinstance(evaluation, list) and training and evaluation, "missing species split")
    _require(all(isinstance(name, str) and name.strip() for name in training + evaluation), "invalid species label")
    _require(len(set(training + evaluation)) == len(training + evaluation), "duplicate or overlapping species split")
    names = set(evaluation)
    integrity = r.get("source_integrity", {})
    _require(all(integrity.get(key) is True for key in ("cite_sha256_verified", "genes_sha256_verified", "all_survivor_headers_occurrences_masks_verified")),
             "source integrity flags incomplete")
    source_species = integrity.get("species", {})
    _require(set(source_species) == set(training + evaluation), "source integrity species mismatch")
    for checks in source_species.values():
        _require(all(checks.get(key) is True for key in ("header_sha_verified", "occurrence_sha_verified", "phase2_mask_sha_verified")),
                 "species source integrity flags incomplete")
    primary, self_result = r["primary_place_beyond_ibd"], r["within_species_self_diagnostic"]
    ps = _scores(primary["heldout_species_scores"], names, "primary")
    ss = _scores(self_result["species_scores"], names, "self")
    _mean_matches(primary["statistic"], ps, "primary aggregate")
    _mean_matches(self_result["statistic"], ss, "self aggregate")
    p = _number(primary["profiled_private_p_value"], "primary p", 0, 1)
    sp = _number(self_result["p_value"], "self p", 0, 1)
    _require(primary.get("alpha") == rule["primary_estimand"]["alpha"] == 0.05, "primary alpha drift")
    _require(rule["primary_estimand"].get("name") == "place_beyond_ibd", "primary estimand drift")
    _require(self_result.get("alpha") == sr["inference"]["alpha"] == 0.05, "self alpha drift")
    _require(primary.get("positive") is (p <= 0.05), "primary significance mismatch")
    _require(self_result.get("positive") is (sp <= 0.05), "self significance mismatch")
    selected, components = primary["selected_configurations"], primary["component_p_values"]
    _require(len(selected) == 2 and len(set(selected)) == 2 and set(selected) == set(components), "profile components mismatch")
    component_p = [_number(v, "component p", 0, 1) for v in components.values()]
    _require(p == max(component_p), "profile p is not the maximum component p")
    qualified = sq.get("passed")
    _require(type(qualified) is bool, "self qualification must be final")
    self_status = "PASS" if qualified else "SELF_DETECTABILITY_NOT_QUALIFIED"
    _require(sq.get("status") == self_status, "self qualification status mismatch")
    _require(self_result.get("synthetic_method_qualified") is qualified and self_result.get("synthetic_qualification_status") == self_status, "self result qualification mismatch")
    _require(a.get("fresh_self_detectability", {}).get("passed") is qualified and a["fresh_self_detectability"].get("status") == self_status, "authorization self qualification mismatch")
    decision = ("TRANSFERABLE_PLACE_COMPONENT" if p <= 0.05 else
                "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN" if qualified and sp <= 0.05 else
                "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING")
    _require(r.get("decision") == decision, "decision inconsistent with primary/self evidence")
    total = r["secondary_total_genetic_transfer"]
    _require(total.get("name") == "total_genetic_distance_transfer" and total.get("inferentially_qualified") is False
             and total.get("used_for_primary_decision") is False, "secondary inference boundary drift")
    _require(total.get("biological_ibd_adjustment") is False and total.get("same_primary_geometry_split_kernel") is True
             and total.get("training_geometry_control") == "within_species_edge_length_rank_orthogonalization",
             "secondary computational contract drift")
    _require(not ({"p_value", "positive", "alpha", "confidence_interval"} & set(total)), "secondary contains inferential labels")
    ts = _scores(total["heldout_species_scores"], names, "total", nullable=True)
    _mean_matches(total["statistic"], ts, "total aggregate")
    undefined = sorted(name for name, value in ts.items() if value is None)
    _require(total.get("status") == ("NOT_EVALUABLE_DESCRIPTIVE" if undefined else "DESCRIPTIVE_ONLY"), "secondary status mismatch")
    _require(total.get("undefined_species") == undefined and total.get("frozen_evaluation_species_count") == len(names)
             and total.get("finite_evaluation_species_count") == len(names) - len(undefined), "secondary species counts mismatch")
    rows = io.StringIO(newline="")
    writer = csv.writer(rows)
    writer.writerow(["species", "post_ibd_transfer", "within_species_self", "total_transfer_descriptive", "total_available"])
    for name in sorted(names):
        writer.writerow([name, ps[name], ss[name], "" if ts[name] is None else ts[name], ts[name] is not None])
    total_text = (f"The descriptive total-transfer score was {total['statistic']:.6g}." if not undefined else
                  f"The descriptive total-transfer aggregate was unavailable for the complete panel ({len(undefined)} non-finite species scores); no reduced-panel mean was substituted.")
    markdown = (
        "# Genetic TTF results\n\n"
        f"The frozen panel comprised {len(training)} training and {len(evaluation)} evaluation species. "
        f"Exact-geometry qualification passed: the maximum private-null Wilson 95% upper bound was {upper:.6g} "
        f"(ceiling 0.10), and the shared-A2 power lower bound was {lower:.6g} (floor 0.80).\n\n"
        f"The species-equal post-IBD transfer statistic was T = {primary['statistic']:.6g} "
        f"(profiled-private p = {p:.6g}; alpha = 0.05). {DECISIONS[decision]}\n\n"
        f"Within-species self-detectability was S = {self_result['statistic']:.6g} (p = {sp:.6g}); "
        f"its separate synthetic qualification was {self_status}.\n\n"
        f"{total_text} This is before biological IBD adjustment while retaining the training-only geometric length control. "
        "It has no significance test and cannot rescue the primary result. Differences between these rank correlations are not fractions of genetic variation explained.\n\n"
        "Inference is limited to the frozen phylogatR COI/COX1 panel, graph, split, marker, scale and qualified reference family. "
        "All frozen evaluation species are retained in the accompanying table.\n"
    )
    return ManuscriptExport(markdown, rows.getvalue(), {
        "schema": "ttf_genetic_manuscript_export_v0.1", "mode": "empirical_summary",
        "decision": decision, "source_sha256": hashes, "geometry_fingerprint_sha256": fingerprint,
        "training_species_count": len(training), "evaluation_species_count": len(evaluation),
        "receipt_consistency_checked": True, "statistical_analysis_rerun": False,
        "sequence_identity_read": False, "input_interpretation_text_used": False,
    })


def build_closed_export(receipt: Path) -> ManuscriptExport:
    payload, digest = _load(receipt)
    schema, status = payload.get("schema"), payload.get("status")
    if schema == "ttf_genetic_empirical_opening_state_v0.3":
        _require(status == "EMPIRICAL_OPENING_CLOSED_SINGLE_EXTERNAL_BLOCKER_FRESH_PHYLOGATR_ARCHIVE", "unexpected opening state")
        firewall = payload.get("global_firewall", {})
        _require(firewall and all(v is False for v in firewall.values()), "opening firewall is not closed")
        endpoint = "AWAITING_FRESH_ARCHIVE"
    else:
        allowed = {
            PREFIX + "confirmatory_phase1_geometry_v0.1": "NOT_EVALUABLE_PHASE1_PANEL_TOO_SMALL",
            PREFIX + "confirmatory_phase2_mask_v0.1": "NOT_EVALUABLE_PHASE2_PANEL_TOO_SMALL",
            PREFIX + "phase3_qualification_v0.1": "NOT_EVALUABLE",
        }
        _require(schema in allowed and status == allowed[schema], "not a supported closed terminal receipt")
        for key in OPEN_FLAGS[:2]:
            _require(payload.get(key) is False, "closed receipt contains empirical outcomes")
        _require(payload.get(OPEN_FLAGS[2], False) is False, "closed receipt contains empirical statistic")
        if schema == PREFIX + "phase3_qualification_v0.1":
            _require(payload.get("passed") is False and payload.get("phase4_identity_opening_eligible") is False, "inconsistent failed qualification")
        endpoint = status
    text = ("# Genetic TTF status\n\n"
            f"Recorded endpoint: {endpoint}.\n\n"
            "No empirical genetic transfer result is available from this receipt. "
            "Nucleotide identity remains unopened. This status does not establish absence of transferable genetic structure. "
            "No empirical Results paragraph or species scores have been generated.\n")
    return ManuscriptExport(text, "", {
        "schema": "ttf_genetic_manuscript_export_v0.1", "mode": "closed_status_only", "decision": endpoint,
        "source_sha256": {"receipt": digest}, "sequence_identity_read": False,
        "statistical_analysis_rerun": False, "empirical_results_generated": False,
    })


def write_export(report: ManuscriptExport, output_dir: Path) -> Path:
    """Create a complete new bundle; never overwrite an existing report."""
    destination = Path(output_dir)
    _require(not destination.exists(), "output directory already exists")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".ttf-report-", dir=destination.parent))
    try:
        name = "results.md" if report.manifest["mode"] == "empirical_summary" else "status.md"
        (temporary / name).write_text(report.markdown, encoding="utf-8")
        if report.species_csv:
            (temporary / "species_scores.csv").write_text(report.species_csv, encoding="utf-8")
        manifest = dict(report.manifest)
        manifest["output_sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in temporary.iterdir()}
        (temporary / "export_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        temporary.rename(destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return destination
