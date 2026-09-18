# Genetic TTF submission readiness

Status: **science complete; human metadata and journal-specific formatting pending**.

This ledger describes the submission state after the single authorized fresh phylogatR Phase-4 opening. It is downstream of the frozen empirical result and cannot authorize any rerun, retuning, new subset, new marker, or stronger biological claim.

## Frozen scientific state

- Terminal decision: `LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN`.
- Primary `place_beyond_ibd`: `T = 0.03256548471362711`, profiled-private `p = 0.7392607392607392`.
- Qualified within-species self statistic: `S = -0.226097652897688`, upper-tail `p = 0.008991008991008992`; the raw sign is interpreted only relative to the frozen structural-null distribution, not against zero.
- Descriptive pre-IBD total-transfer score: `0.015501562599030518`; no significance test and no rescue authority.
- Confirmatory survivor panel: 211 species, split 103 training / 108 evaluation.
- Empirical result SHA-256: `f5d19fa50c7c18cbf2a110c0afb9c70dcd543c44b77521c015938013843e72fb`.
- Exact-manifest Phase-4 authorization SHA-256: `a66edba0a6378f2fa2e24eba73aa5d48f6a58905420e13b53be802f603c3a0e1`.

## Manuscript package

- Main text: `manuscript/genetic_ttf_flagship_v0.2.md`.
- Current title: “Within-species genetic structure is detectable but not detectably transferable across unseen lineages”.
- Current manuscript length: approximately 5,517 whitespace-delimited words including title, headings and references.
- Current abstract: 260 words.
- Frozen figures:
  - `manuscript/figures/genetic_phase4_v0.1/figure1_qualification_margins.svg`
  - `manuscript/figures/genetic_phase4_v0.1/figure2_post_ibd_species_scores.svg`
- Frozen supplementary reporting:
  - `manuscript/generated/genetic_ttf_phase4_v0.1/results.md`
  - `manuscript/generated/genetic_ttf_phase4_v0.1/species_scores.csv` — Table S6, all 108 evaluation species.
  - `manuscript/generated/genetic_ttf_phase4_v0.1/export_manifest.json`
- Terminal empirical receipt: `benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json`.

All reporting artifacts are receipt-derived. Sequence identities and edge-level genetic-distance vectors are not serialized in the manuscript package.

## What remains before a journal submission

These are human/editorial inputs, not scientific-analysis blockers:

1. Final author list, order, affiliations and corresponding-author details.
2. Target journal selection and journal-specific title/abstract/section/reference formatting.
3. Funding, competing-interest, author-contribution and acknowledgement statements.
4. Final data-availability wording for the authenticated source archive under its provenance/redistribution conditions.
5. Journal-specific supplementary-file naming and figure export requirements.

The current 260-word abstract may need a small trim if the selected journal uses a 250-word limit. That is a formatting edit only; it must not change the frozen estimand, result or interpretation.

## Claim boundary for submission

The licensed biological statement is that, within the frozen phylogatR COI/COX1 domain and tested spatial scale, the qualified analysis did not detect cross-species transfer of post-IBD differentiation while the separately qualified within-species self statistic was significant relative to its frozen null. This is **consistent with lineage-conditioned spatial structure within the tested domain**.

The manuscript must not convert this into:

- proof of zero transfer;
- a universal statement across markers, clades or spatial scales;
- identification of a specific historical, environmental or barrier mechanism;
- a causal partition of geography from shared demographic history.

Any mechanistic explanation, independent replication domain, predictor competition, or new marker family is a new prospective study rather than a revision of this one-shot result.
