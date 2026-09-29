# Palearctic terrestrial-insect LGM study v0.1

Status: **response-blind design frozen; genetic subpanel response remains closed.**

## Biological question

The 211-species Phase-4 result asked whether post-IBD spatial genetic structure was transferable across a taxonomically broad animal panel. The present study asks a narrower comparative-phylogeographic question:

> Among Palearctic terrestrial insects, do species that shared more similar climatic geography during the Last Glacial Maximum show more similar present-day spatial genetic discontinuities?

This is not a rescue of the global transfer result and not a reopening of relational v0.3. It is a new, prospectively specified secondary analysis using an existing genetic source.

## Why LGM SDMs change the biological question

The predictor is no longer generic present-day environmental similarity. Independent GBIF occurrences are used to fit a fixed climatic SDM for each response-blind eligible insect species. The model is projected onto CHELSA-TraCE21k climate at 21 ka BP.

The primary pairwise predictor, `R_LGM`, is continuous Schoener overlap between the normalized LGM suitability surfaces of two species. Present-day SDM overlap is mandatory nuisance structure.

Therefore a positive result would mean:

> species whose independently reconstructed LGM climatic suitability occupied more similar parts of the Palearctic also transfer more post-IBD spatial genetic structure between them.

The manuscript must call this **shared LGM climatic-refugial opportunity**, not direct proof of realized refugia.

## Response-blind feasibility

From the exact frozen 211-species authorization list and the exact source archive metadata:

- 183/211 survivors are Insecta;
- a deliberately coarse temperate-Eurasia coordinate screen gives 51 insects with >=80% of exact source occurrence localities inside the proxy;
- restricting that proxy to the predeclared terrestrial-core order set leaves 47 species.

Those numbers are feasibility only. The formal panel will be determined independently using GBIF occurrences and the TNC/WWF terrestrial ecoregion realm field.

## Formal sequence

1. Freeze the present contract before formal Palearctic realm census.
2. Build the response-blind eligible species panel from the exact 211 names, taxonomy, independent GBIF occurrences and the external Palearctic realm polygon.
3. Assign source/target roles by frozen hash.
4. Fit the fixed current-climate SDMs and project them to 21 ka BP.
5. Build `R_LGM`, `R_current`, geographic co-opportunity and declared controls.
6. Bootstrap occurrence inputs to quantify relation repeatability.
7. Run TTF-Q: endpoint-identity survival, nuisance-control survival, signal concentration, null qualification and detectability envelope.
8. **Only if the frozen response-opening domain passes**, compute the new subpanel genetic transfer response and the single one-sided `beta_LGM` test.
9. Otherwise close as `NOT_EVALUABLE_PALEARCTIC_INSECT_LGM` without computing the subpanel response.

## Important boundary

The overall 211-species Phase-4 result is already known. Accordingly, this study cannot be described as an independent fresh-data confirmatory replication. Its protection is narrower: species-level transfer scores and all subpanel response quantities are forbidden design inputs, and the new subpanel, SDM relation and qualification rule are frozen before those subpanel outcomes are computed or inspected.
