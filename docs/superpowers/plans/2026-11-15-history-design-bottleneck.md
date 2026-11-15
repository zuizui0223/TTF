# History design bottleneck Implementation Plan

**Goal:** Locate loss of historical-climate contrast in the original 74-row evaluation geometry and define transparent, idealized targets for a new independent design. No trait outcomes.
**Architecture:** Add deterministic aggregate-only span/cell diagnostics; reuse unchanged exact-source reconstruction and CPU-path identity selection; never alter old terminal decision or fit a new biological model.
**Spec:** This file. Continue the user's instruction to advance the design; no new response opening is authorized.
**Tech Stack:** Python/NumPy/SciPy, existing PR #154. Local isolated scratch checkout; connector is the only repository network access (local clone unavailable).

## Fixed design
- Original source commit 85241f69c0c7fb2ddc7d7c7bbbdd70f4aa7d64bd; exact X/H SHA authority from the terminal receipt.
- Training 91 pairs/29 species; evaluation 74 pairs/27 species/45 studies. Scale from training species-equally only.
- Report all prespecified nuisance spans, not a best model: intercept, environment11, baseline21, species, study, additive species+study, groups+environment11, groups+baseline21, study×species cells, cells+baseline21.
- Rank tolerance 1e-10, column normalization as in history_absorbed_design; verify final rank/energy against the prior function.
- Counts of singleton/repeated cells, repeated-cell species/studies, and within-cell history variation; no identifiers, coordinates or pointwise weights exported.
- Decompose rank/energy along groups -> +environment11 -> +method10, with explicit order dependence. Descriptive only, not causal attribution or response R².
- New-design ideal paired contrasts: two independent error observations per within-species/study block, perfectly matched baseline covariates, standardized history difference d in {0.5,1,2}. I=m*d²/2, residual df=m-1. Fix gamma=.5, sigma=.5, two-sided alpha=1/60, desired power=.8; calculate minimal m via noncentral t. Not a claim that real data satisfy assumptions or a guaranteed sample-size recommendation.
- Existing Haldane v1 remains STOP; no random worlds; no source outcome values; no altered hashes.

## Tasks
- [x] Write tests, confirm explicit missing-implementation failure (core and source runner).
- [x] Implement pure diagnostics and idealized planning math; 17 local tests pass, including independent group demeaning and explicit paired-block design.
- [x] Add source runner/workflow; integration must reproduce old final rank and energy and retain input locks.
- [ ] Execute on GitHub, inspect aggregate artifact, quantify design constraints, record results/citations and no biological inference.

## Review focus
Nonfinite/misaligned inputs; rank-deficient duplicated columns; saturation/no historical variation; singleton vs repeated cells; assumptions and denominator differences in proposed sample counts. Each has a deterministic test.
