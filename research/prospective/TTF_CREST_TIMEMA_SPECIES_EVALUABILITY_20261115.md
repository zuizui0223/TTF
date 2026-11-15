# Timema 2024 × TTF-Q: species-disjoint coverage and stopping decision (2026-11-15)

**Status: SOURCE GEOMETRY PASS, BIOLOGICAL RESPONSE NOT OPENED; CONFIRMATORY CROSS-SPECIES TTF-RI HOLD.** This file records an outcome-blind manuscript triage, not a genetic or reproductive-isolation estimate.

## Verified authority

Original article: Nosil, Gompert & Funk (2024), *Nature Communications* 15:2273, DOI https://doi.org/10.1038/s41467-024-46294-9.

Original publisher source spreadsheets (SHA256 checked):
- [Source Data / Figure3 XLSX](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-46294-9/MediaObjects/41467_2024_46294_MOESM6_ESM.xlsx): `c309402e1b1fa58ab0fb0763086a155aff028cdaa04b804c2c91796c5bed61e2`.
- [Supplementary Data 1–5 XLSX](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-46294-9/MediaObjects/41467_2024_46294_MOESM4_ESM.xlsx): `7e2b6210c182859d7b4d77ad8b0ca0858a1e437c0106a9c42e6c3787c24a33db`.

Both files independently agree on all 42 original comparison identifiers, species-label strings and population-pair cells. No host-choice outcome, sexual-isolation value, or GBS/Fst value was decoded during the source-only qualification. The official Dryad v3 source metadata are pinned, but its raw CSV download returned HTTP 401; the publisher originals are the separately verified public authority for this check.

## What the comparison graph actually contains

- **42 dyad rows**: **32 within-label** and **10 between-label** comparisons.
- **10 distinct literal source species-name labels**, not 10 independently verified current taxonomic species. The manuscript describes *nine sexually reproducing species*. This discrepancy is unresolved and **must not be silently reconciled** by dropping or merging an outcome-dependent taxon.
- The source comparison graph separates into **five disconnected components**:
  1. `bartmani–chumash–podura`: 19 dyad rows, 4 between-label comparisons.
  2. `cristinae`: 11 dyads, all within-label.
  3. `californicum–poppensis`: 7 dyads, 4 between-label.
  4. `knulli–landelsensis–petita`: 3 dyads, 2 between-label.
  5. `curi`: 2 dyads, all within-label.
- Thus, there are only **three disjoint source components containing between-label comparisons**. The 10 interspecific edges are not 10 independent replicated mechanisms of speciation.

## Heldout-geometry results

The exact code **first** enumerated 490 species-label-disjoint candidate partitions with nonempty train and evaluation dyads (fixed anchor removes complementary duplicates).

An additional *source-only candidate design gate* required at least:
- 3 distinct source species labels per arm;
- 3 within-label comparisons per arm;
- 2 between-label comparisons per arm.

Among 456 partitions meeting the arm-label-size requirement, **103** retained both within- and between-label comparisons in both arms. A balanced maximum-minimum comparison-count candidate is:

| Source-only measure | Train | Unseen-label evaluation |
|---|---:|---:|
| Source species labels | 4 | 6 |
| Total dyad rows | 21 | 21 |
| Within-label dyads | 17 | 15 |
| Between-label dyads | **4** | **6** |
| Effective diversity of species endpoint *appearances* (inverse Simpson, NOT effective sample size) | **2.1565** | **2.8544** |
| Highest endpoint share | **64.29%** | **52.38%** |

Training labels: `bartmani, chumash, curi, podura`; evaluation labels: `californicum, cristinae, knulli, landelsensis, petita, poppensis`. Zero dyads cross these source components. The **21/21 balance obscures the 4/6 between-species data bottleneck and dominant-source composition**. Partition selection is based solely on the source geometry; a model using the same source as an all-taxon response-trained model would not be an independent validation.

## Experimental response limits remain unresolved

The article reports host preference assays for all **42** dyad comparisons (**8,749 trials**) and sexual-isolation assays for only **30** comparisons (**2,436 trials**). It also reports genetic data for only a subset (e.g., 36 comparisons in Fig. 4), including some estimates inferred from an alternative marker rather than directly genotyped GBS.

The publisher Source Data's Figure3 rows serialize 42 cells under the "sexual isolation" column. **A serialized XLSX cell is not proof of a measured sexual-isolation outcome.** Those cells can include missing-value strings/formulas and must be checked at least for validated numerical/experiment support before opening the biological response. The extra storage-type-only audit does not read any effect sizes. The true overlap of sexual-isolation **observations** and observed-not-imputed GBS across species-disjoint arms remains unknown.

No claim of 42 sexual outcomes, 10 independent species, 103 statistically powered valid designs, or positive/negative TTF-RI biological transfer is warranted.

## TTF × CREST conclusion and next disposition

**Geometry: PASS. Empirical cross-species RI transport: HOLD.** Next source-only gate: authenticate response availability, strict taxonomic units, GBS direct-versus-inferred provenance, and dyad connectivity under both arms. Only then freeze a single constrained heldout model and calibrate its error/power before opening a reproductive-isolation predictor.

If the 10 source labels cannot be reconciled with the nine-sexual-species authority, or neither arm retains sufficient independent between-label comparisons with actually measured sexual RI and observed genomic differentiation, terminate as `NOT_EVALUABLE`, not an ecological null. The claim that habitat isolation behaves differently within and between species is **already an original result of Nosil et al. 2024**, not a new discovery here.

## Reproduction

- [Qualified source geometry machine receipt](timema_species_qualified_source_geometry_receipt_20261115.json)
- [Original publisher source graph code](ttf_crest_timema_publisher_species_holdout_geometry_20261115.py)
- [Species- and response-availability metadata-only qualifier](ttf_crest_timema_species_qualification_20261115.py)
- [Successful independent Actions run #43696875531](https://github.com/zuizui0223/TTF/actions/runs/43696875531)
- Artifact ID `14155335792`, ZIP SHA256 `2e767eb60230c8b369af3b345db3975514f6a82e88f03bfb42feb1e859cbf718`.

No frozen TTF-Q, COI, He, Fst, break-location response or prior Senecio P0/P1 result was changed.
