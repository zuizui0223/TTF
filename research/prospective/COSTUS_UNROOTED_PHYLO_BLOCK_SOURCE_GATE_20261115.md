# Costus TTF source-only amendment: unrooted phylogenetic distance blocks (2026-11-15)

**Status:** registered before opening any 2025 floral scent response values. Retained source SHA256 are unchanged. Prior *rooted clade* qualification workflow #43700817708 ended in `NOT_EVALUABLE_ORIGINAL_CHRONOGRAM_NOT_ROOTED`. An apparent Newick tree layout is **not** sufficient evidence for a biologically rooted phylogeny; do not force a root or claim independent origins of hummingbird pollination.

## Scientific purpose

A defensible **out-of-phylogenetic-neighborhood** TTF-P pilot, not a held-out independent transition / true adaptation mechanism test. This independent source-derived grouping can be used later for a species-level heldout scent model only if source geometry is adequate and a single new outcome evaluation is separately preregistered.

## Sources

- 2022 Kay & Grossenbacher chronogram `chronogram_v8_2.tre`, byte-identical 2025 OSF mirror, SHA256 `24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1`.
- 2025 Darragh et al. OSF `Floral Scent Analyses/phylo_info.csv`, SHA256 `9ce8f35b736652d3f0bbd2970707ce1aebebc0c4786962a23a983da7e81a8226`.
- Source species-level response `means_species.csv`, SHA256 `056b91173cd14907347557c267c802fa644c710830d039370b1b416e755f1af1` has not had scent magnitudes read. Its unnamed taxon-key column has exactly 27 keys, and all are present among 69 original tree tips and the 27 source pollination labels.
- Publisher source pollination labels: 14 `Bee`, 13 `Hummingbird`, **not independently certified actual visitation records**, because some 2022 pollination states were inferred from morphology.

## Exact geometry-only algorithm, prospectively frozen before any scent responses

1. Parse the original tree with DendroPy preserving the original unrooted flag. **Never call reroot, root-on-outgroup or force rooted.** Read only the 27 matching source taxon tips and their *patristic* pairwise distances (without interpreting the tree as rooted).
2. Initialize 27 singleton taxon groups. Use deterministic **average-linkage agglomerative clustering** on patristic distances. For each step merge the two groups with the smallest arithmetic mean cross-group tip distance, breaking numeric ties after 12-decimal rounding lexicographically on original normalized tip labels.
3. Examine the finite *predeclared* ordered group-count family `K = 4, 3, 2`. Use the **first** K with **every heldout group** containing at least three original source tips, at least one original source `Bee` label and at least one `Hummingbird` label, and every complementary training set containing at least four tips of each class. If none passes, report `SOURCE_GEOMETRY_NOT_EVALUABLE` and stop.
4. The prospective evaluation unit would be a **source-tip / taxon**, not a scented individual, pollen deposition event, genetic sequence, or independent evolutionary transition. This source file includes multiple accessions or infraspecific forms; a strict species-level identity gate must ensure any repeated species do not cross train and test before calling it species-heldout.
5. Never use scent compound richness, total amount, gene expression, or any response-derived grouping to select K or neighborhood assignments.

## Boundaries

- These are *unrooted phylogenetic neighborhood blocks*, **not** uniquely reconstructed independent pollinator-shift lineages or molecular phylogenetic root states.
- All 2025 scent data have been analyzed already by the original authors; a new fixed blocked prediction estimator will be a **methodological retrospective reanalysis**, not independent discovery of the relationship between hummingbird pollination and scent loss.
- Comparison across three source papers still requires taxon/voucher ID reconciliation and the direct-observed-versus-morphology-inferred pollinator label distinction.
- This amendment is an explicitly separate response-blind source qualification, not a post-outcome rescue, because no floral scent or gene-expression values were opened in this sequence.
