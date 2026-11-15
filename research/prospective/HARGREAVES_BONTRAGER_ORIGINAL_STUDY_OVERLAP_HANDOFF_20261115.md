# TTF — Hargreaves × Bontrager original-source overlap (source-only terminal handoff)

**2026-11-15 | Draft PR #155 | new source-only provenance result | no fitness/genetic/phenotypic-rate response opened**

## Decision

**Original-study-label overlap is directly measured, not inferred from different dataset DOIs.** The publicly preserved original Hargreaves CSV contains 148 plant study labels; **102 of those exactly match 102 of the 103 study labels** in the previously verified 107-plant Bontrager reference catalog. Of those 102 identically named source studies, **100 share at least one identical plant taxon text label**, with **128 exact study × taxon keys**. The remaining 46 Hargreaves plant study labels do **not** match this *qualified* Bontrager reference list; they are **not** independently certified experiments, and their apparent novelty must not be used to rescue the Bontrager historical-fitness 44/50-species source gate.

The raw **source-study label** is not proof that the exact same experiment, treatment, source-population pair or fitness measurements were reused. Such finer deduplication requires original experiment IDs, study paper IDs and source/garden identity concordance. Conversely an unmatched label is not proof of an independent original experiment.

## Original published source identities

- Hargreaves et al. 2020, *American Naturalist*, original data DOI https://doi.org/10.5061/dryad.0vt4b8gv1.
- Official Zenodo Dryad preservation record https://zenodo.org/records/4409291 has an **explicit exact original DOI** and lists the actual CSV file `Hargreaves_etal_2020_AmNat_transplant_data.csv`.
- Published original CSV **4,382,308 bytes**, MD5 **`3e8fde16540509341e8f3c1a096c267d`**, downloaded-byte SHA256 **`f00eec59e5388da5a2243d391c4e639dc77eda8d98abbf9a8a6846b14da7bb71`**.
- Source original CSV header is 47 columns, SHA256 **`2ea31cf076d3053096a97374b9b1ff7347d9f759328a3cbcc40ced7c3092df6b`**. It has `study, functional_group, kingdom, taxon, source, site` among its fields but lacks `lat_source, long_source`. It has only destination `lat_site, long_site` and derived distances, insufficient to uniquely reconstruct source locations and origin-specific glacial climates.
- Source CSV contains **11,088 physical rows**, including **10,810 original annual/herbaceous-perennial/woody-perennial plant rows**. These are repeated observations, not independent replicates.
- Source labels: **160 total study labels, 148 plant study labels, 162 plant taxon text labels**; no claim that all 162 are valid distinct biological species.

Bontrager reference: original DOI https://doi.org/10.5061/dryad.pnvx0k6mc, Borealis mirror https://doi.org/10.5683/SP2/DU3ZTQ; verified source catalog contains **103 original study labels, 107 plant binomial species, 711 garden×fitness-type cells**. Source catalog SHA256 `fdfb9a3e825935dc4abc602466142b10bcc8fccbba210b89e65f0575a66fabe9`. Exactly **one Bontrager study label** did not occur among the Hargreaves *plant* study labels: `Smithetal_BiolCons_2005`.

## Original-record source-only overlap

| Measure | Result | What it means |
|---|---:|---|
| Bontrager qualified plant original study labels | 103 | Reference set |
| Hargreaves original plant study labels | 148 | Includes repeated study/population/focal rows |
| Exact **study label** matches | **102** | Direct source metadata overlap; not experimental identity certification |
| Matched studies sharing ≥1 same named plant taxon | **100** | Stronger source/organism evidence but not population-level identity |
| Exact **study × taxon** matched keys | **128** | Labels from shared primary publications; non-independent counts |
| Hargreaves plant study labels unmatched to Bontrager reference | **46** | Unlinked; **not** 46 verified independent experiments |
| Normalized-only study label candidates | 0 | Exact study ID equality carries main source-level match |
| Independently certified additional experiments | **0** | Certification not performed; do not interpret as empirical absence |

All counts are by unique original source ID or study×taxon name, not by arbitrary row or garden-pair enumeration. The Hargreaves source data's "kingdom" includes nonplants; plant-only selection used the three original **precommitted functional_group** categories, not a selection based on the biological outcome.

## Implementation / provenance

1. Strict source search: original Dryad DOI must be present as a *structured exact identifier* of a Zenodo record (not title similarity); record **4409291** located by [Actions #43694395913](https://github.com/zuizui0223/TTF/actions/runs/43694395913).
2. Original 47-column first-line-only source schema recovered with zero biological data rows read by [Actions #43694961822](https://github.com/zuizui0223/TTF/actions/runs/43694961822).
3. Frozen source-only 47-column contract: `hargreaves_zenodo_study_overlap_contract_v01.json`; full original bytes checked against published MD5 before row parsing. Strict UTF-8 was invalid, so the immutable original was parsed using **strict CP1252**, under separately logged [technical encoding amendment](./hargreaves_original_encoding_amendment_v01.json). No imputation, no replacement characters.
4. Only `study,functional_group,kingdom,taxon,source,site` fields are selected. Biological response, treatment-effect, model coefficient and fitness columns were **not accessed or interpreted**. The entire raw file was byte-hashed; do not confuse this with opening biological fitness outcomes.
5. Actual source-only study-ID overlap **PASS**: [Actions #43695463005](https://github.com/zuizui0223/TTF/actions/runs/43695463005). Five local parser tests, eight original header tests and eight index tests pass within this workflow.
6. Exact SHA-pinned 7,670-byte JSON result receipt saved permanently at [source_registry/hargreaves_bontrager_original_study_overlap_v01.json](./source_registry/hargreaves_bontrager_original_study_overlap_v01.json); original receipt SHA256 **`c279c042b1386ebb35af38045e80cc7047c4cd7820aba19ce332d6f980f121f6`**. It contains matched source labels and unmatched labels, no raw population geocoordinates, source IDs, observed fitness values, or row-level effects.

## Scientific interpretation and stopping boundary

**Supported:** Hargreaves and Bontrager share almost all of Bontrager's original qualified plant study labels. Different dataset DOIs cannot be counted as independent sources, and the Hargreaves file is missing the original source coordinates required for the proposed independent-source 21ka-versus-0BP climate-history comparison.

**Not supported:** genetic adaptation, historical climate effects on fitness, phenotypic plasticity, reproductive isolation, speciation, or an effect of LGM on local adaptation. These biological responses were never fitted in this route.

**Frozen and untouched:** PR #154's terminal Bontrager **44/50** source qualification; PROCEED 56-plant synthetic null calibration failure; older genetic He, Fst, COI, Jensen and TTF-Q endpoints. PR #155 remains DRAFT and must **not** be merged as a biological fitness prediction. The historical predictor's 60.12% residual source variation is not 60.12% fitness explained.

**Future independent program (not part of this exhausted source-screen family):** Identify genuinely different original field transplant experiments with original *source-population* lat/lon and shared-garden fitness comparisons. Before opening biological fitness, preregister primary study/experiment deduplication, sampling/missingness rules, species/study-disjoint heldout evaluation and synthetic spatial false-positive/power calibration. Genetic versus plastic effects require independent common-garden or replicated reaction-norm genotype evidence; speciation requires divergent lineages/reproductive isolation response data.

**Terminal source-only verdict for this pair:** `HARGREAVES_BONTRAGER_SHARED_SOURCE_STUDIES_NO_INDEPENDENT_HISTORICAL_FITNESS_RESCUE`.
