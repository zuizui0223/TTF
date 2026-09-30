# Palearctic terrestrial-insect LGM extension v0.1

## Biological question

> Do Palearctic terrestrial insects that shared similar climatic refugial opportunity at the Last Glacial Maximum show genetic differentiation in more similar places?

This is intentionally narrower than the 211-species flagship question. The original flagship asks whether a general spatial genetic field transfers across unseen lineages. This extension asks whether transfer becomes biologically predictable after restricting the system to a coherent region and life-history domain and explicitly representing late-Quaternary climatic history.

## Why this is not a rescue of the 211-species result

The overall 211-species Phase-4 result is already known. Therefore this same-archive extension cannot be described as an independent confirmatory experiment.

The valid contribution is narrower: freeze a subgroup and an LGM predictor **without inspecting subgroup-specific species scores or pair-specific transfer values**, characterize that design with TTF-Q, and compute the subgroup genetic response only if the predeclared qualification rule passes.

A failure of the panel, ecological data, predictor geometry, null qualification or detectability closes the route as NOT_EVALUABLE.

## Response-blind source

The authoritative species universe is the compact 211-species Phase-2 survivor list frozen from the pre-opening authorization. Only taxonomy, COI-family FASTA headers, linked occurrence coordinates and external ecological data may be used before this route is authorized to open a genetic response.

A coarse non-authoritative screen performed before this contract found 183 Insecta among the 211 survivors and roughly 48 species under a deliberately conservative terrestrial + approximate Palearctic screen. That count has no membership authority. The formal filter is fixed below and may still fail.

## Formal panel

The panel is deliberately conservative.

- Class must be Insecta.
- Allowed terrestrial-core orders: Lepidoptera, Hymenoptera, Coleoptera, Hemiptera, Orthoptera and Psocodea.
- Diptera and orders with characteristically aquatic immature stages are excluded rather than adjudicated after seeing results.
- A small fixed blacklist removes aquatic Coleoptera families.
- At least 80% of exact selected-COI localities must intersect RESOLVE 2017 polygons labelled Palearctic.
- At least 40 species must remain, yielding at least 20 source and 20 target clusters after a fresh hash-based role assignment.

There is no backfill.

## LGM model

The LGM layer is a deterministic presence-only climatic-envelope SDM rather than a result-selected MaxEnt tuning exercise.

Independent GBIF occurrences are filtered and thinned before climate extraction. Current CHELSA bio1/bio7/bio12/bio15 define a common four-dimensional whitened climate space. Each species receives a Gaussian KDE in that space. The same frozen KDE is evaluated on current and 21-ka CHELSA climate over a fixed Palearctic land grid.

The normalized current and LGM suitability maps then provide two continuous dyadic relations:

- R_present: Schoener D between current suitability maps;
- R_LGM: Schoener D between LGM suitability maps.

R_LGM is primary. R_present is mandatory nuisance.

This avoids thresholding a refugial map into an arbitrary binary refuge before qualification. Thresholded connected refugia can be descriptive later, not a replacement primary predictor.

## TTF-Q comes before genetics

Before any pair-specific genetic transfer response is computed for this route, TTF-Q must answer:

1. how much R_LGM survives source/target identity;
2. how much survives present climate, geography, taxonomy and sampling controls;
3. whether residual signal is endpoint-concentrated;
4. where the intended crossed-dyad test satisfies the declared Type-I ceiling;
5. the MDE envelope over nuisance heterogeneity.

The genetic response is opened only if the frozen utility rule in the machine-readable contract passes.

## If it passes

The empirical test is:

T_st ~ source FE + target FE + R_LGM + R_present + geography + same order + same family + occurrence-count imbalance.

The directional prediction is beta_LGM > 0.

A positive result would turn the genetic TTF story into a concrete comparative-phylogeographic mechanism: spatial structure is not universally transferable, but lineages sharing similar LGM refugial opportunity reuse geography more similarly.

A qualified null would still be biologically interpretable: even within a coherent Palearctic terrestrial-insect domain, LGM climatic refugial similarity does not explain transfer heterogeneity under the frozen representation.
