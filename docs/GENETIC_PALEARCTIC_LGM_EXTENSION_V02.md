# Palearctic non-Lepidoptera LGM extension v0.2

## Question

> Do Palearctic terrestrial holometabolous insects outside Lepidoptera that occupied similar LGM climatic refugial space show post-IBD genetic differentiation in more similar places?

This supersedes the 44-species v0.1 panel **before any v0.1 LGM relation, TTF-Q result, or subgroup genetic response existed**.

## Why v0.1 was stopped

The formal response-blind v0.1 census retained 44 species, but 26 were Lepidoptera. That would make the paper substantially overlap the independent butterfly ecology programme.

The v0.2 scope therefore excludes Lepidoptera by design. This is a manuscript-domain decision made without species-level genetic response, pairwise subgroup T_st, beta_LGM, or any LGM-relation result.

The earlier v0.1 census remains in the repository as provenance and is not rewritten.

## v0.2 panel

The source universe remains the exact 211 Phase-2 survivors. The new panel is:

- Insecta;
- Lepidoptera excluded;
- terrestrial Holometabola restricted to Coleoptera, Hymenoptera and a conservative whitelist of terrestrial Diptera families;
- aquatic Coleoptera families excluded;
- at least **12 exact selected-COI localities inside the RESOLVE 2017 Palearctic realm**.

Unlike v0.1, Palearctic *fraction* is descriptive only. A Holarctic species can enter if its Palearctic subgeometry itself contains at least 12 localities. This matches the biological question: we need enough within-Palearctic geography to evaluate an LGM-refugial hypothesis; we do not need the species to be Palearctic-endemic.

The response-blind census gives **23 species**:

- Coleoptera: 14
- Hymenoptera: 8
- Diptera: 1

No backfill is allowed.

## Dyadic design

Every retained species appears as both source and target. Self-dyads are excluded.

With all 23 species retained, the ecological design contains **506 directed dyads** and 23 source plus 23 target endpoint clusters. This is preferable to an arbitrary 11/12 role split for this new pair-specific ecological question.

External GBIF/climate missingness may only remove species. At least 20 species must survive; otherwise the route is NOT_EVALUABLE.

## LGM model

The climate representation remains the already frozen response-blind v0.1 design:

- independent Palearctic GBIF occurrences;
- CHELSA V2.1 current bio1/bio7/bio12/bio15;
- CHELSA-TraCE21k 21 ka BP corresponding variables;
- pooled four-dimensional whitened climate space;
- species-specific Gaussian KDE climatic envelopes;
- current and LGM suitability evaluated on the same fixed Palearctic grid;
- primary pair relation R_LGM = Schoener D between LGM suitability surfaces;
- R_present = Schoener D between current suitability surfaces and is mandatory nuisance.

This is a deterministic presence-only climatic-envelope SDM/hindcast. Thresholded refuge polygons are not allowed to replace the continuous primary relation after seeing results.

## TTF-Q before genetics

The original v0.1 TTF-Q opening thresholds are inherited **without relaxation**:

- total unique R_LGM fraction >= 0.15;
- effective source and target signal counts >= 10;
- maximum source and target signal share <= 0.15;
- null qualification passes at A=0,1,2;
- evaluable A=2 grid MDE <= 0.10.

If any fails, the study stops without recomputing v0.2 pairwise genetic response.

## What would be opened only after PASS

The response would be directed source-to-target post-IBD transfer T_st recomputed on the Palearctic-clipped v0.2 genetic geometry. The primary coefficient is beta_LGM > 0 after present-climate similarity, current geographic co-opportunity, same-order, same-family, occurrence-count imbalance, and source/target identity are controlled.

Because the parent 211-species archive has already been opened historically, this remains a pre-specified **secondary ecological reanalysis**, not a fresh independent confirmation.
