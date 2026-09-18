# Molecular Ecology submission adapter — genetic TTF

Status: **candidate first-shot adapter; journal choice remains a human editorial decision**.

This file adapts the frozen genetic TTF manuscript to the current Molecular Ecology submission surface without changing the empirical analysis, estimand, null family, species set, figures, result or biological interpretation.

## Why this journal is a natural fit

Molecular Ecology explicitly includes population structure, phylogeography and analytical methods development within scope. The genetic TTF paper uses molecular genetic data to ask a general comparative-phylogeographic question and couples a predictive analytical framework to a prospectively qualified empirical test.

Current author-guideline constraints checked for this adapter:

- Original Articles: up to 8,000 words for Abstract + Introduction + Materials and Methods + Results + Discussion.
- Abstract: no more than 250 words.
- Keywords: four to six.
- Cover letter: include a statement explaining fit to the journal's aims and scope.
- Data accessibility and author-contribution statements are required at submission/production stages.

The canonical empirical manuscript remains `manuscript/genetic_ttf_flagship_v0.2.md`; this adapter is formatting only.

## Candidate title

**Within-species genetic structure is detectable but not detectably transferable across unseen lineages**

## Candidate abstract — 213 words

Comparative phylogeography seeks spatial genetic patterns shared across species, but retrospective concordance does not show whether geography learned from some lineages predicts differentiation in unseen lineages. We define shared geographic structure as out-of-species transferability. For each species, exact sampling localities were connected by a density-scaled graph, genetic differentiation was summarized on graph edges, and a leave-two-localities-out nuisance fit removed ordinary isolation by distance before a spatial turnover field was learned from training species. A fresh phylogatR COI/COX1 analysis was prospectively gated so that geometry, character admissibility and exact-geometry synthetic qualification were fixed before nucleotide identities were opened. Phase 1 froze 250 species; 211 survived the character-support gate and retained a 103-training/108-evaluation split. The survivor geometry passed predeclared cross-species calibration and separate within-species self-detectability qualification. In the single authorized empirical opening, post-IBD transfer to unseen species was not detected (T = 0.0326, profiled-private p = 0.7393), whereas the calibrated within-species self statistic was significant relative to its frozen structural-null reference (S = -0.2261, upper-tail p = 0.0090). The pre-IBD total-transfer score was 0.0155 and remained descriptive only. Thus, within this COI/COX1 domain, spatial genetic structure was detectable within species but not detectably reusable across unseen lineages under the tested geographic field, consistent with lineage-conditioned spatial structure rather than a universally transferable map of differentiation.

## Candidate keywords

- comparative phylogeography
- isolation by distance
- spatial population genetics
- predictive transfer
- COI/COX1
- cross-species generalization

## Cover-letter scope statement

This manuscript fits Molecular Ecology because it uses molecular genetic data to address a general question in population structure and comparative phylogeography: whether the geography of within-species genetic differentiation learned from some species predicts differentiation in entirely unseen species beyond isolation by distance. The study combines a new predictive analytical framework with a prospectively qualified empirical COI/COX1 test across 211 species, while explicitly separating method qualification, within-species detectability and cross-species transfer. The result is broadly relevant to how comparative phylogeography distinguishes reusable properties of place from lineage-conditioned spatial structure.

## Remaining human inputs for a Molecular Ecology submission

- authors, order, affiliations and corresponding author;
- funding and acknowledgements;
- competing-interest statement;
- author-contribution statement;
- final data-accessibility wording for the authenticated source archive;
- confirmation of the exact journal file designations and figure/export format at submission time.

## Frozen boundary

The adapter must not turn the qualified non-significant cross-species result into an equivalence claim, proof of zero transfer, or evidence for a specific historical mechanism. It also must not treat the negative raw self-statistic sign against zero; that statistic is interpreted relative to its prospectively frozen structural-null distribution.
