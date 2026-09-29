# Palearctic-core Lepidoptera × LGM-SDM re-analysis v0.1

## Biological question

Do terrestrial insects that shared similar Last Glacial Maximum suitable geography also show more similar spatial genetic differentiation after ordinary isolation by distance is removed?

The working empirical system is the response-blind subset of the frozen 211-species phylogatR panel whose selected COI sampling geometry is both taxonomically homogeneous (Lepidoptera) and confined to a conservative temperate Palearctic core.

## Why this is separate from the 211-species flagship claim

The flagship asks whether one geographic field transfers across a broad heterogeneous set of unseen lineages.

This study asks a conditional comparative-phylogeographic question:

> does historical geographic similarity explain which lineage pairs share spatial genetic structure?

The predictor is pair-specific and historical. The intended response, if later authorized, is pair-specific source-to-target post-IBD genetic transfer/concordance.

## Opening order

1. **P0 — response-blind panel census.**
   Use only frozen survivor names, taxonomy, aligned FASTA headers and occurrence coordinates. Do not inspect nucleotide identity, genetic distance, Phase-4 species scores or subpanel transfer values.

2. **P1 — LGM SDM construction.**
   Build present climate SDMs from external GBIF occurrences, project each to CHELSA-TraCE21k 21 ka, and construct the pairwise LGM suitability-overlap relation.

3. **P2 — TTF-Q qualification.**
   Characterize endpoint survival, unique LGM information beyond present climate/geography/taxonomy, residual signal concentration, null qualification and conditional detectability.

4. **STOP if not evaluable.**
   No genetic subpanel result is opened.

5. **P3 — one frozen genetic re-analysis if evaluable.**
   Construct directed source-to-target post-IBD genetic cross-prediction/concordance and test its relation to frozen R_LGM.

## Important provenance limitation

The 211-species full-panel Phase-4 genetic response was opened historically before this new question was conceived. Therefore a later P3 test would be a prospectively specified **re-analysis conditional on a known global aggregate result**, not a fully untouched confirmatory experiment.

The protection against response-driven subgroup fishing is narrower but real: P0-P2 and the LGM predictor are frozen without reading species-level or pairwise subpanel genetic responses.

## Current P0 result

A strict Palearctic-core Lepidoptera panel exists and exceeds the predeclared continuation minimum. The frozen receipt is:

`benchmarks/frozen/palearctic_lepidoptera_response_blind_screen_v0.1.json`.

The next permitted scientific action is P1 LGM predictor construction only.
