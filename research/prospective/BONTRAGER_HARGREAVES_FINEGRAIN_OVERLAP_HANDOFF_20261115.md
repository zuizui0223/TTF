# TTF — Bontrager × Hargreaves: independent data or duplicated source designs?

**2026-11-15 | Draft PR #155 | response-blind original-source comparison | biological fitness, adaptation, plasticity, speciation outcomes NOT fitted**

## Main conclusion

The two publicly released plant transplant meta-analysis datasets **are not independent empirical panels at the primary-study/source-population/garden design level**. Different dataset DOIs do not establish independent experimental replication.

Earlier, a comparison to the **qualified** 103-study Bontrager reference found exact Hargreaves study labels for **102/103**, leaving **46** Hargreaves plant study labels unlinked to that *qualified subset*. The new first-shot, SHA-locked, source-only audit of **both full original plant files**, rather than only the Bontrager eligible subset, changes the source-overlap accounting:

| Source-only unit | Bontrager full plant original | Hargreaves full plant original | Exact shared identifier keys |
|---|---:|---:|---:|
| Study label | 114 | 148 | **113** |
| Study × binomial taxon | 143 | 174 | **137** |
| Study × taxon × original source ID | 1,554 | 1,670 | **1,511** |
| Study × taxon × transplant garden ID | 538 | 615 | **502** |
| Study × taxon × source × garden | 5,100 | 5,373 | **4,947** |
| Above × original performance **type** | 6,751 | 7,333 | **6,536** |

The full Bontrager original has **11 plant study labels excluded from its earlier qualified 103-label reference**. All 11 also appear in Hargreaves, so previously "unlinked" Hargreaves labels decrease from 46 to **35**. Importantly, **35 is NOT the count of confirmed independent new experiments**: it is the count of Hargreaves study labels not found in any Bontrager full original plant record. These labels and their primary-paper/experiment IDs would need additional verification.

Of **6,751** unique Bontrager source×garden×fitness-type metadata keys, **6,536 (96.8%)** occur verbatim in Hargreaves. These are labels of experimental design and measurement category, **not duplicate fitness numeric values, raw biological replicates, test observations, or species-independent predictive performance**.

## Geographic cross-check

Within the **502** shared study×plant taxon×garden-ID keys, original **transplant garden** WGS84 coordinate agreement (unique finite point on both sides) was:

- **488** within 1 km (**97.2%** of the 502 shared site keys)
- **4** farther than 1 km but within 10 km
- **9** farther than 10 km
- **1** missing or ambiguous original site coordinate

The 13 >1 km differences cannot be ignored; they may reflect geocoding revisions, different localizations or a mistaken equivalence of garden IDs. No raw coordinates or site names were exported. The two datasets therefore provide **strong metadata evidence of source-design reuse**, not universal positional identity. Hargreaves's original 47-field CSV does **not** contain the original source-population latitude/longitude fields required to reconstruct independent LGM exposure, even where source IDs match.

## Source provenance and why this does not reopen adaptation responses

1. Original Bontrager source: DOI https://doi.org/10.5061/dryad.pnvx0k6mc; official Zenodo archived original at https://zenodo.org/records/4666308. Original dryad.zip 811,204 bytes, MD5 `2de98c1bb2cd2ddffcb090eb65a74d9f`, SHA256 `2995d0f165b27ccf2f38f3d0b286845c107762a0c2b02df2290e9a5df6b396f0`.
2. Original Hargreaves source: DOI https://doi.org/10.5061/dryad.0vt4b8gv1; official Zenodo archived original at https://zenodo.org/records/4409291. Hargreaves original CSV 4,382,308 bytes, MD5 `3e8fde16540509341e8f3c1a096c267d`, SHA256 `f00eec59e5388da5a2243d391c4e639dc77eda8d98abbf9a8a6846b14da7bb71`.
3. Pre-data comparison rule: [frozen contract](./bontrager_hargreaves_finegrain_overlap_contract_v01.json), exact string equality after trimming `study/source/site/fitness_type`, strict binomial normalization of taxon; no fuzzy author matching, response-derived tuning or changed qualification cutoffs.
4. Two original file-byte digests were checked **before reading** allowed metadata fields. Only `study,functional_group,taxon,source,site,fitness_type,lat_site,long_site` fields were selected from the respective 47-column CSVs. Hargreaves required strict CP1252 decoding; Bontrager original was UTF-8-sig. No original numeric fitness values, effect sizes, genotype or climate-fitness relationships were selected or modeled. Unselected response fields are necessarily tokenized during CSV parsing; they were not inspected or used in qualification.
5. Original source-only crosswalk [GitHub Actions #43696872468](https://github.com/zuizui0223/TTF/actions/runs/43696872468) **SUCCESS**: 9 synthetic known-truth tests passed and source metadata join completed. The frozen aggregate-only receipt from this original run has SHA256 `8268387cdaec3a835086418f7ce91f91777cbbc1026e148c6e6e0399647106c3`. The receipt contains no raw source or garden coordinates and must remain the authority, not any retuned rerun.

## Scientific boundary and next distinct program

- **Supported:** shared source study labels, many shared source-population/garden/metric-type identifiers and mostly concordant original transplant-garden coordinates. Substantial provenance overlap precludes claiming two separate meta-analysis releases provide independent large comparative adaptation evidence.
- **Not yet supported:** identical original experiments in every shared key, independent experimental replication for 35 unmatched study labels, genetic adaptation to LGM, plasticity decomposition, speciation timing or lineage divergence. This audit contains **no phenotype/fitness regression**.
- **Already stopped and NOT reopened:** PR #154 historical source-match Bontrager **44 of 50 taxa** macro gate; PROCEED 56-plant synthetic calibration gate; previously frozen He/Fst/COI/Jensen/TTF-Q outcome tests.
- **Only defensible next empirical target:** an entirely new prospective panel with authenticated **independent original primary experiments**, at least two original **geographically identified source populations per same garden**, original source lat/lon and source-type outcome available only after response-independent qualification, primary experiment-ID adjudication, species/study-disjoint heldout evaluation and synthetic null/power calibration. **Do not append 35 Hargreaves study labels or 11 Bontrager formerly ineligible labels to the exhausted prior cohort as a rescue.**

**Source-only terminal status**: `SHARED_SOURCE_GARDEN_DESIGNS_CONFIRMED_NO_INDEPENDENT_HISTORY_FITNESS_RESCUE`. PR #155 remains DRAFT, with no biological-response or main-branch merge.
