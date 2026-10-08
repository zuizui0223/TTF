# Historical host connectivity: a fresh test of biotic memory in Lepidoptera

Status: **prospective; no fresh genetic response opened**

## 2026-10-08 prequery integrity correction — primary panel v0.2

**Current state: HOLD external execution pending a new v0.2 authority. No genetic response has been opened.**

An independent replay of the exact frozen candidate CSV and WCVP v13/HOSTS sidecar exposed two prequery issues.

1. The committed v0.1 panel CSV contained `36964` for the *Fraxinus excelsior* host of *Craniophora ligustri*. The exact frozen sidecar uses `369664`. Rebuilding the v0.1 panel from frozen inputs reproduced **140 species** and the **original receipt SHA-256 `67e8d6ad...` byte-for-byte**; the corrupted GitHub CSV was restored to those originally frozen bytes.
2. The original exclusion universe included the 211 response-opened Phase-4 survivors, S1/S2 and fresh-1000, but omitted non-surviving members of the **original 250-species Phase-1 design** and the **conditional 500-species design**. The original 250 were rederived response-blind and contained all 211 survivors. The complete strict design exclusion union has **2,492 species**. Five v0.1 insects overlapped the previously designed conditional cohort: *Apoda y-inversum*, *Chionodes electella*, *Herpetogramma aeglealis*, *Nisoniades rubescens*, and *Stigmella hybnerella*.

The corrected v0.2 panel removes exactly these five insects, **without backfill**, yielding **135 insects and 222 accepted host IDs**. The retained insects' host IDs, host names and relative deterministic order are unchanged. Receipt: `benchmarks/frozen/genetic_historical_host_prequery_integrity_audit_v0.2.json`. The old 140-species external execution authority **must not** be applied to the 135-species panel.

The two previous authoritative v0.1 replacement workflows stopped before Neotoma or GBIF external queries. Subsequent runs triggered by restoring the old CSV, if any, are not authoritative for v0.2. No Neotoma coverage, host occurrence outcome, LGM predictor, or genetic response was consulted to select the five excluded species.

**An additional biological gate is needed before any LGM host-resistance analysis.** Exact present-day HOSTS associations do not establish that the same insect–host interactions existed at 21 ka. The 135-species panel includes cultivated host taxa (e.g., *Zea mays*, *Malus domestica*, *Prunus domestica*, *Brassica oleracea*). These cannot simply be hindcast as if modern domesticated taxa and their present consumer associations existed during the LGM. An independently sourced temporal-eligibility protocol must be frozen and applied response-blind before the historic-interaction prediction is interpreted.

The licensed scientific question remains historical **host-resource opportunity**, not a directly observed ancient interaction network. Neotoma pollen at genus level can indicate past host-lineage occurrence but cannot alone establish the exact host species or its use by the insect.

## Independent cross-design freshness audit (2026-10-08)

The corrected 135-species panel was checked against the previously **response-blind, genetically unopened** Study C historical-design list (1,000 species). Four insects are shared: *Erebia aethiops*, *Glena cribrataria*, *Hypena laceratalis*, and *Teleiopsis diffinis*.

These are **design-only overlaps**, not genetic-response reuse: Study C closed as `NOT_EVALUABLE` before nucleotide identity or source-target genetic transfer was opened. The present study remains genetically fresh, but it is not disjoint from every historical response-blind design universe. The four species remain in the frozen 135-species panel; none may be removed, substituted or backfilled based on paleo-coverage or future genetics. Machine-readable audit: `benchmarks/frozen/genetic_historical_host_cross_program_overlap_audit_v0.1.json`.

The host-genetic positive prediction also requires the independent **insect-abiotic LGM** counterfactual frozen in `docs/supporting/genetic_historical_host_biotic_abiotic_contrast_v0.2.json`. A host effect without that control cannot establish biotic historical memory. Existing Neotoma runs are external, descriptive fossil-evidence censuses, not empirical genetic results.

## Question

> Do present-day genetic discontinuities in narrow-host Lepidoptera remember where their larval host plants were fragmented during the Last Glacial Maximum?

The primary prediction is deliberately within species rather than another source-target similarity test.

For every fresh insect species, genetic edges are defined before nucleotide identity. Larval host plants are obtained from a fixed HOSTS snapshot and resolved to accepted WCVP species. Only insects with 1-5 accepted species-level hosts enter the focal domain.

Each host plant receives a response-blind palaeoclimatic suitability model. Host suitabilities are combined into an insect-specific resource surface at present and at 21 ka BP. The same frozen insect graph edge is then assigned present and LGM host-mediated resistance.

The primary within-species model asks whether **LGM host resistance remains positively associated with post-IBD genetic turnover after present host resistance is controlled**.

## Why this is biologically different from the failed host-resource relation

The completed host-resource study asked whether two Lepidoptera species with similar **current host-resource footprints** transferred spatial genetic information better. It returned a qualified null.

The present study does not reuse those genetic outcomes and does not ask whether species are ecologically similar.

Instead it asks whether a dependent consumer's own population structure carries a lagged spatial imprint of where its interaction partners could persist during the last glacial cycle.

That is a historical biotic-constraint hypothesis.

## Freshness

The pre-census excludes, by species identity, the complete 211-species opened Phase-4 survivor set, the complete S1 design, the complete S2 design, and the entire fresh-1000 design universe used by later relational/codistributed work.

Even under that deliberately strict exclusion, the exact phylogatR archive supplies the first **600** geometry-qualified fresh Lepidoptera after checking only 1,084 ranked candidates.

HOSTS then gives species-level larval hosts for **230/600** species. **146** have 1-5 species-level hosts before WCVP synonym resolution, so a >100-species narrow-host focal panel is realistically available without touching genetic response.

## Primary biological contrast

The key comparison is not host versus no host and not specialist versus generalist.

It is:

- **present host resistance** — where the larval resource landscape is discontinuous now;
- **LGM host resistance** — where that same interaction-constrained landscape was discontinuous at 21 ka.

If only current host geography matters, the LGM term should disappear after the present term is included.

If historical interaction structure leaves a genetic memory, the LGM term should remain positive.

## Palaeo validation

Neotoma pollen and plant-macrofossil records can provide a secondary, direct palaeoecological check for host genera where dated fossil coverage exists. They do not select the confirmatory result and do not rescue species after genetic opening.

## Claim if positive

> **The spatial genetic structure of dependent herbivores retains a memory of the past geography of their interaction partners.**

A stronger causal claim about an exact refugium or fossil population would require additional evidence.

## Current evidence state

Everything reported so far remains response-blind.

The strict fresh precensus excludes the complete opened 211-species Phase-4 survivor set, the full S1 and S2 design universes, and the entire fresh-1000 universe. From the remaining archive, 600 geometry-qualified Lepidoptera were frozen. Exact HOSTS × WCVP resolution then produced a **140-species narrow-host panel** with 1-5 accepted species-level larval hosts per insect.

The primary panel contains **230 distinct frozen accepted host IDs**. Panel membership, host breadth and host IDs are fixed; there is no backfill.

The first external-data workflows exposed a transport/identity inconsistency before any Neotoma or GBIF query: downstream workflows attempted to reconstruct the frozen host IDs from WCVP tables rather than replay the exact WCVP/HOSTS sidecar bytes that defined the panel. The scientific panel did not fail. A response-blind sidecar-replay repair is now frozen and binds the exact original sidecar artifact and SHA-256 values. Replacement Neotoma and GBIF runs were designated before their results were inspected.

The remaining opening sequence is now frozen:

1. exact native-range host occurrence gate;
2. current/LGM host-resistance construction;
3. predictor-information gate;
4. deterministic development/confirmatory split;
5. confirmatory canonical/noncanonical character-mask gate;
6. exact survivor synthetic Type-I/power qualification;
7. one-shot confirmatory nucleotide-identity opening only after every prior gate passes.

The exact response-blind geometry materialization rule for the 140-species panel is also frozen so that locality and edge geometry can be extracted once from the authenticated phylogatR archive without serializing sequence identity.

Still unopened:

- fresh nucleotide identity;
- fresh pairwise genetic distance;
- fresh post-IBD turnover;
- fresh historical-host coefficient;
- any host-history/genetics association.

Thus no biological result yet exists for historical biotic memory. A later gate failure is `NOT_EVALUABLE`, not a biological null.
