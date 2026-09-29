# Palearctic insect LGM comparative phylogeography v0.1

## Biological question

> Do Palearctic terrestrial insects that shared Last Glacial Maximum climatic
> refugial opportunity show genetic discontinuities in the same geographic
> places?

This is intentionally narrower than the 211-species genetic TTF manuscript.
The global result asked whether one spatial field transfers broadly across
unseen lineages. The present study asks whether a specific historical
mechanism explains why transfer should be stronger for some species pairs than
for others.

## Why LGM SDMs enter before genetics

For each eligible insect species, a present-day climatic SDM is fitted from
external GBIF occurrences and projected to 21 ka using the same
CHELSA-TraCE21k variables. The pairwise focal relation is LGM suitability
overlap. Present-day suitability overlap is mandatory nuisance information.

The design therefore distinguishes:

`shared present niche` from `shared LGM climatic opportunity`.

TTF-Q is run on that predictor geometry before any species-level genetic
subpanel response is inspected.

## Response

The genetic response is not the old global PASS/FAIL transfer statistic.

Each species receives a spatial break-intensity surface from positive
post-IBD genetic residuals on the already frozen genetic graph, using the
inherited 500-km kernel. For a source-target species pair, response is Schoener
D between those two break surfaces.

The primary model asks whether pairwise break-surface concordance increases
with LGM suitability overlap after current suitability overlap, geographic
co-opportunity, lineage, sampling imbalance, and exact endpoint identity are
removed.

## Analysis order

1. Start from the exact 211 Phase-2/Phase-4 survivor identities.
2. Use taxonomy and genetic-panel coordinates only to define the Palearctic
   terrestrial-insect panel.
3. Freeze a balanced deterministic source/target split.
4. Acquire fresh external GBIF occurrence data.
5. Fit and validate present SDMs; project the same models to 21 ka.
6. Freeze `R_LGM`, `R_current`, geography, lineage and sampling controls.
7. Run TTF-Q.
8. If the study-specific authorization rule fails, stop.
9. Only after PASS, construct genetic break surfaces and open the pairwise
   response once.

## Status

The response-blind source, eligibility rule, LGM-SDM relation definition and
genetic response definition are frozen. Authoritative WWF realm assignment,
SDM construction and TTF-Q qualification remain unopened.
