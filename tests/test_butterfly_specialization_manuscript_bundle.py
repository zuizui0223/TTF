from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _json(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _words(text: str) -> int:
    return len(re.findall(r"\b[\w–'-]+\b", text))


def test_butterfly_specialization_manuscript_bundle_is_internally_consistent():
    manuscript_path = ROOT / "manuscript/butterfly_specialization_ecology_v0.1.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.1.json")
    readiness = _json(
        "manuscript/butterfly_specialization_submission_readiness_v0.1.json"
    )
    figures = _json("manuscript/butterfly_specialization_figures_v0.1.json")

    dimensionality = _json(
        "benchmarks/exploratory/butterfly_specialization_dimensionality_result_v0.1.json"
    )
    anthropogenic = _json(
        "benchmarks/exploratory/butterfly_anthropogenic_resource_expansion_result_v0.1.json"
    )
    mechanism = _json(
        "benchmarks/exploratory/butterfly_resource_expansion_mechanism_result_v0.1.json"
    )
    hierarchy = _json(
        "benchmarks/exploratory/butterfly_host_specialization_hierarchy_result_v0.1.json"
    )
    climate = _json(
        "benchmarks/exploratory/butterfly_climate_release_postgate_independent_result_v0.1.json"
    )

    assert claim_map["status"] == "WORKING_SUBMISSION_CLAIM_MAP"
    assert readiness["target_journal"]["primary"] == "Global Ecology and Biogeography"
    assert readiness["target_journal"]["article_type"] == "Research Article"
    assert figures["status"] == "RENDERED_AND_VISUALLY_AUDITED"
    assert figures["workflow"]["run_id"] == 36292128188
    assert len(figures["figures"]) == 5

    # Structured GEB abstract should remain below 300 words.
    abstract_start = manuscript.index("## Abstract")
    abstract_end = manuscript.index("---", abstract_start)
    assert _words(manuscript[abstract_start:abstract_end]) <= 300

    # Approximate main-text budget remains under the GEB typical 5,000 words.
    intro_start = manuscript.index("## 1. Introduction")
    legend_start = manuscript.index("## Figure legends")
    assert _words(manuscript[intro_start:legend_start]) <= 5000

    # Submission sections expected in the working bundle.
    for heading in (
        "## 1. Introduction",
        "## 2. Methods",
        "## 3. Results",
        "## 4. Discussion",
        "## 5. Limitations",
        "## 6. Conclusions",
        "## Figure legends",
        "## Data and Code Availability",
        "## References (working)",
    ):
        assert heading in manuscript

    # Core source values used in the paper remain identical to frozen receipts.
    assert (
        dimensionality["host_taxonomy_lower_bound_adequate_subset"][
            "spearman_host_family_vs_geographic_resource_breadth"
        ]
        == 0.398238154462201
    )

    full = anthropogenic["full_resource_eligible_panel"]
    assert full["species"] == 239
    assert full["species_expanded_by_introduced_host_ranges"] == 206
    assert (
        full["spearman"]["host_family_count_vs_log_resource_expansion"]
        == 0.00798541214928876
    )

    adequate_mech = mechanism["host_taxonomy_lower_bound_adequate_subset"]
    assert (
        adequate_mech["spearman"][
            "host_family_count_vs_effective_contributor_number"
        ]
        == 0.4919670615030713
    )
    assert (
        adequate_mech["spearman"][
            "host_family_count_vs_maximum_single_host_fractional_share"
        ]
        == -0.4861161332504768
    )

    one_family = hierarchy["primary_focus"]
    assert one_family["species"] == 82
    assert (
        one_family["spearman"][
            "resolved_host_species_vs_effective_contributor_number"
        ]
        == 0.7338967715221224
    )
    assert (
        one_family["spearman"][
            "resolved_host_species_vs_maximum_single_host_fractional_share"
        ]
        == -0.7065837667974557
    )

    assert climate["climate_crossfit"]["climate_informative_species"] == 24
    assert climate["climate_crossfit"]["species_with_score_above_neutral_0_5"] == 23
    assert climate["primary_test"]["observed_partial_spearman"] == -0.16643611123731772
    assert climate["primary_test"]["one_sided_p_value"] == 0.2237
    assert (
        climate["primary_test"]["decision"]
        == "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"
    )

    # Headline numbers must actually appear in the submitted prose.
    for literal in (
        "206/239",
        "rho = 0.008",
        "rho = 0.734",
        "rho = -0.707",
        "median filtering score = 0.801",
        "partial rho = -0.166",
        "p = 0.2237",
    ):
        assert literal in manuscript


def test_manuscript_claim_map_preserves_inference_boundaries():
    claim_map = _json("manuscript/butterfly_specialization_claim_map_v0.1.json")
    claims = {row["id"]: row for row in claim_map["claims"]}

    assert claims["C1"]["status"] == "exploratory descriptive"
    assert claims["C2"]["status"] == "exploratory descriptive"
    assert claims["C3"]["status"] == "exploratory post-result mechanism analysis"
    assert claims["C4"]["status"] == "exploratory post-result follow-up"
    assert claims["C5"]["status"] == "independent cross-fit descriptive result"
    assert claims["C6"]["status"] == "independent hypothesis test"
    assert (
        claims["C6"]["decision"]
        == "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"
    )

    forbidden = "\n".join(claim_map["forbidden_overclaims"])
    assert "caused butterfly range expansion" in forbidden
    assert "true absence" in forbidden
    assert "confirmatory claims" in forbidden
