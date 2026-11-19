# TTF | 独立移植原研究の重複監査（response-blind preflight）

**2026-11-15 | fresh source-only route | empirical fitness/phenotype unopened | NOT an extension of stopped PR #154**

## Decision

Existing [PR #154](https://github.com/zuizui0223/TTF/pull/154) ends with 44/50 eligible Bontrager plant taxa under its frozen current-climate/history contrast, so **the 44-taxon cohort must not be rescued** by adding later reviews, weakening thresholds, or opening the fitness outcome. The 44 matched taxa, 37 studies, 319 eligible same-garden cells, and 60.12485% *predictor* residual history fraction are design information, **not fitness explained or adaptive evolution**.

This new branch implements only a **bibliographic source-overlap gate** for a genuinely new, separately preregistered panel. A different meta-analysis dataset DOI is not evidence of different primary experiments.

## Verified external source register (publication-level, no raw outcomes)

| Source | Public dataset DOI | Source-level observation | Eligibility / overlap judgement |
|---|---|---|---|
| Bontrager et al. (2021) | [10.5061/dryad.pnvx0k6mc](https://doi.org/10.5061/dryad.pnvx0k6mc) | Original raw study × species × source × garden data; original ZIP mirrored and SHA-verified in PR #154 | Reference cohort; original 44/50 gate already terminal |
| Hargreaves et al. (2020) | [10.5061/dryad.0vt4b8gv1](https://doi.org/10.5061/dryad.0vt4b8gv1) | Public Dryad methods explicitly say much of its literature data was originally collected for Bontrager et al.; geolocation conversion is reported | **Suspected substantial overlap**; do not treat as an independent replication until original study IDs/DOIs crosschecked |
| Hereford (2009) | [10.5061/dryad.315](https://doi.org/10.5061/dryad.315) | Historical meta-analysis in an XLS file; population coordinate support unverified | Original-paper overlap and original-source geography **unknown** |
| Briscoe Runquist et al. (2020) | [10.5061/dryad.2v6wwpzh5](https://doi.org/10.5061/dryad.2v6wwpzh5) | Published effect-size tables on abiotic/biotic local adaptation | Source-population coordinate/pair availability and primary-paper overlap unverified; effect summaries alone cannot be treated as source-level histories |
| Dittmar & Schemske (2023) | [10.5061/dryad.6djh9w10v](https://doi.org/10.5061/dryad.6djh9w10v) | Reciprocal transplant of *Leptosiphon parviflorus*, two habitat populations ~100 m apart | Useful single-species mechanistic candidate but not a many-species independent panel; raster-resolution history contrasts unverified |
| Carley et al. (2025) | [10.5061/dryad.f1vhhmh24](https://doi.org/10.5061/dryad.f1vhhmh24) | Reciprocal transplant and demography across 10 populations of one plant species | One species only; possible mechanistic/fitness replication, not macro species-count rescue |

All source statements above come from publisher/Dryad dataset descriptions. This table is **not a quantitative record-level overlap census**; original primary-study identifiers were not downloaded or matched in this branch.

## Implemented overlap gate

The script [primary_study_overlap_gate.py](./primary_study_overlap_gate.py) accepts **only an independently prepared, metadata-only CSV** with this **exact** header:

~~~csv
source_dataset_doi,source_study_id,primary_paper_doi,primary_first_author,primary_year,primary_title,plant_taxon
~~~

It refuses extra columns (including fitness, effect size, coordinates), malformed lines, missing original-study identifiers, inconsistent study identities, or mixing datasets in one manifest. It checks:

1. Exact original-paper DOI or normalized title/year matching: **bibliographic overlap**.
2. Same first-author/year without exact DOI/title match: **possible overlap; manual review**.
3. No bibliographic linkage: **unlinked, NOT certified independent**.
4. A shared dataset DOI: error, not a second independent cohort.
5. Source study IDs count once across repeated taxa; no pair or garden pseudoreplication is mistaken for independent studies.

Run with sanitized, source-only manifests:

~~~bash
python research/prospective/primary_study_overlap_gate.py \
  --reference bontrager_study_metadata_only.csv \
  --candidate new_source_study_metadata_only.csv \
  --output overlap_audit.json

python -m unittest discover \
  -s research/prospective -p test_primary_study_overlap_gate.py -v
~~~

The output intentionally sets "verified_independent_primary_studies" to **zero**. It is a screening and abstention instrument; **manual adjudication of original experiments** (same paper, accession, cohorts, population IDs, gardens, and year) and independent provenance must precede any certification. A different paper may reuse the same experiment, while one paper may report multiple true experiments.

## Scientific endpoint separation (future, not evaluated)

- **Historical source contrast**: current-source climate matched but local 21ka-to-0BP displacement different. These locations do not establish ancestral occupancy or refugia.
- **Local adaptation**: independently sampled origin-by-garden fitness and appropriate genetic/maternal controls; preregistered response model, study/experiment-level clustering, species-heldout evaluation, false-positive/power calibration and source ascertainment analysis. Fitness contrasts alone do not prove selection by LGM history.
- **Plasticity**: replicated common-garden/reaction-norm measurements across environments for the same genotypes or lineages, preferably with generation or maternal controls. Genetic/phenotypic database labels or transplant fitness alone are insufficient to decompose plasticity.
- **Speciation**: separate lineage divergence/reproductive isolation/timing response data; no inference from short-timescale fitness differences or source climate geometry.

**Next empirical dependency:** Obtain original-study *metadata only* from independent sources and map original paper DOI/experiment IDs for manual overlap resolution. Then freeze a genuinely new response-independent data-qualification contract. Do not reuse or requalify the failed Bontrager or PROCEED 56-plant gates. No biological outcome was opened in this branch.


## 2026-11-15 executed source preflight: public-source lineage yes, record-level overlap unknown

**Original primary-source provenance from Hargreaves, NOT just an inference based on authorship.** Hargreaves's publisher supplementary Appendix B states that the project **leveraged the Bontrager et al. in-preparation transplant experiment database**. That compilation's literature-search inputs explicitly include Hereford (2009), along with Leimu & Fischer (2008), Lee-Yaw et al. (2016) and others. This establishes a linked *source-selection pipeline*, **not a measured count** of identical original experiments or duplicate source-population records. Hargreaves's Dryad record separately states that much of the information was originally collected for Bontrager.

- Hargreaves Dryad dataset: https://doi.org/10.5061/dryad.0vt4b8gv1
- Original publisher supplementary Appendix B, printed p. 9: https://mobot-biodiversity-jc.weebly.com/uploads/1/8/6/0/18603232/appendices_from_hargreaves_et_al__2020__local_adaptation_to_biotic_interactions_59113supplement.pdf
- Bontrager Dryad source: https://doi.org/10.5061/dryad.pnvx0k6mc
- Hereford Dryad source: https://doi.org/10.5061/dryad.315

**Original CSV access test.** On [GitHub Actions #43672366333](https://github.com/zuizui0223/TTF/actions/runs/43672366333), the exact file listed by Dryad as Hargreaves_etal_2020_AmNat_transplant_data.csv (Dryad file 244459) was probed with a byte-limited request: REST API GET yielded **HTTP 401**, public file-stream URL **HTTP 403**. The probe only attempts the first CSV header; it did not read or fit fitness rows, and did not guess underlying original-study IDs. The two source-only programs passed **16/16 synthetic tests** (10 bibliographic overlap, 6 strict source-header access). The workflow succeeded as a correctly detected access gate, **not as a successful data download**.

Source-only decision: **SOURCE_ACCESS_GATED_NO_RECORD_LEVEL_OVERLAP_CLAIM**. A confirmed shared literature-data lineage means these datasets **cannot automatically be counted as independent replication**, but the **number** of duplicate or nonduplicate primary experiments remains unknown (not observed). This is a distinct state from an ecological null or a data census of zero duplicates. The primary-source metadata branch remains Draft, and the old Bontrager 44/50, PROCEED 56-plant and genetic endpoints are unchanged.

Receipt: [hargreaves_bontrager_access_lineage_receipt_v01.json](./hargreaves_bontrager_access_lineage_receipt_v01.json).

**Independent repository-wide CI caveat.** The earlier PR-wide pytest run found **485 passed / 1 failed** for test_ttf_q_submission_scientific_freeze_binds_current_science: frozen-policy string "new outcome opening" versus test expectation "opening a previously closed outcome". This is unrelated to the new source-only code; neither the science freeze JSON nor the global test was modified in this branch.
