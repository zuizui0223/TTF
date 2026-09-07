# Independent biogeographic benchmarks

TTF must not be validated only on the flower-colour system that motivated its development. Empirical qualification therefore uses boundaries and traits defined independently of RGFCA and of TTF outcomes.

## Anti-circularity rule

A named biogeographic boundary is **not** an input to the primary TTF field.

For every positive benchmark:

1. freeze the source dataset, eligible species, sampling units, graph rule, genetic/trait distance, train/evaluation rule, and external boundary geometry before looking at TTF outcomes;
2. fit the boundary field using training species only and without a named-boundary indicator;
3. score entirely held-out species at the edge level;
4. run sharedness inference on held-out species scores;
5. only after that result is fixed, compare transition support with the independently defined biogeographic boundary.

This separates:

- **TTF sharedness:** does transition structure learned from some species transfer to unseen species?
- **biogeographic concordance:** does transferable structure coincide with an independently defined biogeographic break?

The second cannot be used to tune the first.

## Primary genetic candidate: Amazonian 66-species UCE panel

Johnson et al. (2023; Molecular Ecology 32:2186–2205; doi:10.1111/mec.16886) archived a comparative Amazonian genomic dataset on Dryad (doi:10.5061/dryad.rxwdbrvc1). The source dataset contains UCE data from 587 samples of 66 Amazonian bird species, approximately divided among upland forest, floodplain forest, and riverine-island habitats. The Dryad deposit includes a small sampling metadata file, a trait table, a README, and a 32.48 MB VCF archive in addition to larger alignment files.

This panel is now the primary candidate for the **species-level statistical benchmark** because it supplies many more independent species units than the Iranian nine-species panel.

### Metadata-first contract

Before any VCF-based TTF outcome is opened:

- retrieve and freeze the Dryad file/version inventory;
- inspect only the README, variable key, coordinate-bearing sampling table, and other non-genetic-outcome metadata needed for eligibility;
- reconcile specimen IDs and species names against the public source-code repository;
- determine per-species sample counts and usable coordinate support;
- freeze a geometry-only eligibility rule;
- run shared/private semi-synthetic calibration on the retained real specimen geometry.

Species cannot be included or removed because their FST, Dxy, clustering, published population structure, or eventual TTF score is convenient.

Only if actual-geometry calibration passes are VCFs acquired for TTF genetic-distance construction.

### Why this panel is not yet a named-boundary positive

The Johnson dataset was assembled to compare genetic structure across Amazonian habitat types, not to provide 66 independent replicates of one pre-labelled river barrier. TTF therefore treats it first as a **boundary-blind transferable-transition panel**. Any later concordance with the Amazon/Solimões, Madeira, Negro, or another biogeographic feature must use independently frozen geometry after sharedness inference has been fixed.

A smaller white-sand-bird dataset is retained separately as a biological concordance panel for the Amazon River.

## Iranian Plateau RADseq panel: small-panel stress benchmark

Noroozi et al. (2026; Molecular Ecology 35:e70355; doi:10.1111/mec.70355; PMCID PMC13126619) sampled nine endemic mountain plant species across 113 populations. The published material provides Table S1 coordinates and nine species-specific RADseq SNP datasets. TTF acquired and audited these inputs before opening empirical sharedness outcomes.

Population counts in Table S1 and Appendix S1 NEXUS individual identifiers were reconciled across all nine species. The source geometry was then used with **synthetic trait worlds only**.

### What the real-geometry stress test showed

The idealized TTF v0.2 inference had already passed high-precision synthetic qualification, but the frozen Iranian 3-train/6-evaluation geometry did not.

The raw v0.2 actual-geometry pilot gave:

- maximum zero-shared rejection: **0.15**;
- fully shared amplitude-2 power: **0.48**.

Several prospective attempts were then compared without opening empirical SNP sharedness outcomes:

| Variant | max type-I | A=2 full-shared power |
| --- | ---: | ---: |
| 3-fold cross-fit, raw | 0.32 | 0.68 |
| 3-fold cross-fit + IBD nuisance | 0.25 | 0.74 |
| single 3/6 split + IBD nuisance | 0.08 | 0.50 |
| paired opportunity subtraction, raw | 0.05 | 0.24 |
| paired opportunity subtraction + IBD | 0.02 | 0.11 |
| predictor orthogonalization to opportunity | 0.10 | 0.45 |
| predictor orthogonalization to opportunity + edge length | 0.07 | 0.33 |

The repeated pattern matters: methods that control false positives recover only about 0.3–0.5 power at the prospectively defined moderate fully shared signal. The panel is therefore not being tuned until positive.

### Interpretation

The Iranian dataset is retained as an **applicability/detection-floor stress test**. It demonstrates that a biologically attractive comparative panel can still be too small for a species-level superpopulation claim. Empirical SNP sharedness remains sealed for TTF inference.

The useful result is methodological: TTF must distinguish “no shared structure” from “insufficient number of independent species units.”

## Genetic distance contract

For any genetic benchmark, TTF stores pairwise distances as a fixed symmetric matrix. **Distance-matrix cells are never permuted.** Graph edges are selected solely from geography and then read their corresponding genetic dissimilarities.

Where population allele frequencies are the observation unit, Rogers distance is the prospective primary metric and Nei distance is a robustness metric unless the source-data structure requires a different contract to be frozen before outcomes are inspected.

Raw genetic distance can contain isolation by distance. A positive raw genetic TTF result establishes transferable genetic turnover, not automatically a discrete barrier mechanism. Barrier attribution requires an separately qualified distance nuisance or paired predictor comparison.

## External biological positive: Amazonian white-sand birds

The seven-species UCE panel of Amazonian white-sand birds is reserved for external biological concordance. The source study reports species-specific variation in genetic structure while identifying the Amazon River as the barrier shared across all seven species. This closely matches the biological idea behind TTF, but seven species are too few to serve as the sole primary species-level inferential panel after the Iranian stress-test experience.

It can instead answer a later question: does a transferable field discovered in a larger Amazonian panel align with a known cross-species river break?

## Other trait domains

### Phenology

A Pacific Northwest herbarium panel provides a route to a temporal trait. Raw day-of-year is not interpreted directly as a boundary trait. TTF should first remove a prospectively specified within-species expectation from year, temperature, and elevation, then analyze residual phenological turnover or fitted breakpoint parameters. The target is whether a maritime–interior regime change transfers among species.

### Body size and morphometrics

The Great Plains suture zone has strong morphometric/genetic hybrid-zone examples. Many published replicates are taxon-pair hybrid systems rather than populations within a single nominal species, so this remains a secondary system-level benchmark until a compatible replicate-unit adapter is qualified.

### Song dialects

Song is valuable because its generating mechanism is cultural rather than genomic or pigmentary. Location-blind acoustic embeddings can define pairwise song distance without hand-assigning dialect labels. Current Andean examples are strongest as single-species demonstrations, so they remain mechanism prototypes until a multi-species audio panel is frozen.

## Negative controls

### Point Conception

Point Conception is an important marine biogeographic transition but published comparative phylogeography has found that intraspecific genetic breaks do not consistently coincide with it. It is therefore an informative negative control: an important community-level boundary must not automatically force a positive TTF result.

### Wallace Line

The Wallace Line is a scope negative. It strongly separates faunas, but species replacement alone is not the TTF estimand. A TTF positive is expected only if enough individual species themselves span the line and exhibit transferable within-species transitions.

## Claim ladder

1. **Idealized synthetic qualification:** distinguish private strong structure from shared structure.
2. **Actual-geometry qualification:** demonstrate adequate type-I control and detection power with the intended species counts and coordinates.
3. **Independent trait benchmark:** recover transferable structure in a non-flower trait without using a named boundary in training.
4. **External boundary concordance:** compare the frozen transferable field with an independently defined biogeographic boundary.
5. **Mechanistic attribution:** distinguish smooth distance, geography/barriers, environment, and resistance using held-out predictor increments.
6. **Cross-domain replication:** repeat in a second biological trait domain.

An empirical TTF sharedness claim cannot skip level 2.
