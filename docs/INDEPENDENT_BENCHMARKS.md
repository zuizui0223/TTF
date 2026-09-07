# Independent biogeographic benchmarks

TTF must not be validated only on the flower-colour system that motivated its development. Its empirical qualification therefore uses boundaries and traits that were defined independently of RGFCA and of TTF outcomes.

## Anti-circularity rule

A named biogeographic boundary is **not** an input to the primary TTF field.

For every positive benchmark:

1. freeze the source dataset, eligible species, sampling units, graph rule, genetic/trait distance, train/evaluation split rule, and external boundary geometry before looking at TTF outcomes;
2. fit the boundary field using training species only and without a boundary indicator;
3. score entirely held-out species at the edge level;
4. run the sharedness inference on held-out species scores;
5. only after that inferential result is fixed, compare transition support with the independently defined biogeographic boundary.

This separates two questions:

- **TTF sharedness:** does transition structure learned from some species transfer to unseen species?
- **biogeographic concordance:** does that transferable structure coincide with a previously defined biogeographic break?

The second cannot be used to tune the first.

## Primary executable positive: Iranian Plateau RADseq panel

Noroozi et al. (2026; Molecular Ecology, doi:10.1111/mec.70355; PMCID PMC13126619) sampled nine endemic mountain plant species across the Iranian Plateau. Each species spans at least three major mountain regions and was sampled at 8–15 localities. The published analysis identified phylogeographic break zones supported by several species, notably between Alborz and Zagros and between the Azerbaijan Plateau and Zagros. Those mountain regions were already defined as areas of endemism from floristic data, so their boundaries are independent of the focal genetic distances.

TTF will use:

- Table S1 for population coordinates and population identities;
- Appendix S1/S2 for RADseq SNP genotypes;
- population-level genetic distances as the trait dissimilarity object;
- a within-species geographic graph among populations;
- species-disjoint training and evaluation;
- held-out-species inference as the primary sharedness test.

### Genetic distance contract

The first benchmark will prospectively calculate two population genetic distances where the supplied genotype format permits it:

- Nei distance;
- Rogers distance.

One will be designated primary before TTF outcomes are inspected; the second is a robustness metric. If reconstruction from the archived SNP format cannot exactly reproduce the source definitions, the discrepancy is recorded rather than silently switching metrics.

TTF stores pairwise distances as a fixed symmetric matrix. **Distance-matrix cells are never permuted.** Graph edges are selected solely from geography, then read their corresponding genetic distances from the matrix.

### Isolation-by-distance caveat

A positive raw genetic-distance TTF result establishes transferable genetic turnover, not automatically a discrete barrier mechanism. Geographic distance can itself create genetic differentiation. Therefore:

- raw genetic-distance transfer is the benchmark for the sharedness estimator;
- attribution to named mountain breaks requires a separately qualified IBD nuisance layer or paired predictor comparison (`D` versus `DG`/barrier features);
- nuisance fits must be cross-fitted so that evaluation species do not tune their own expected IBD relationship.

## Independent replication: Amazonian white-sand birds

The seven-species UCE panel of Amazonian white-sand birds is reserved as an independent positive replication. The source study reports species-specific variation in genetic structure while identifying the Amazon River as the only barrier shared across all seven species. This is conceptually close to the TTF estimand—heterogeneous within-species amplitude with a transferable transition location—but it is kept separate from the primary Iranian Plateau development benchmark to avoid tuning on multiple known positives.

## Broader discovery panel: Amazonian birds

The 66-species Amazonian UCE dataset is a discovery/robustness panel rather than a prespecified single-boundary positive. Species eligibility must be determined from coordinates and sampling support before any TTF score is inspected. Species may not be added or dropped because they improve the sharedness statistic.

## Other trait domains

### Phenology

A Pacific Northwest herbarium panel provides a route to a temporal trait. Raw day-of-year should not be interpreted as a boundary trait. TTF should first remove a prospectively specified within-species expectation from year, temperature and elevation, then analyze residual phenological turnover or fitted breakpoint parameters. The key target is whether a maritime–interior regime change transfers among species.

### Body size and morphometrics

The Great Plains suture zone is biologically important and has strong morphometric/genetic hybrid-zone examples. However, many published replicates are taxon-pair hybrid systems rather than independent populations within a single nominal species. It is therefore a secondary system-level benchmark until a compatible replicate-unit adapter is qualified.

### Song dialects

Song is especially valuable because its generating mechanism is cultural rather than genomic or pigmentary. Location-blind acoustic embeddings can define pairwise song distance without assigning dialect labels by hand. The best-known Andean examples are currently strongest as single-species demonstrations, so they are treated as mechanism prototypes until a multi-species panel is frozen.

## Negative controls

### Point Conception

Point Conception is an important marine biogeographic transition but published comparative phylogeography has found that intraspecific genetic breaks do not consistently coincide with it. It is therefore an informative negative control: an important community-level boundary must not automatically force a positive TTF result.

### Wallace Line

The Wallace Line is a scope negative. It strongly separates faunas, but species replacement alone is not the TTF estimand. A TTF positive is expected only if enough individual species span the line and exhibit transferable within-species transitions there.

## Claim ladder

1. **Synthetic qualification:** distinguish private strong structure from shared structure.
2. **Independent trait benchmark:** recover transferable structure in a non-flower trait without using a named boundary in training.
3. **External boundary concordance:** show that the learned transferable field aligns with a prospectively frozen biogeographic boundary.
4. **Mechanistic attribution:** distinguish smooth distance, geography/barriers, environment and resistance using held-out predictor increments.
5. **Cross-domain replication:** repeat in a second biological trait domain.

Only levels 1–2 are required to establish that TTF is a working sharedness estimator. Levels 3–5 support progressively stronger biological interpretations.
