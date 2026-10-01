# TTF-Q submission scientific freeze v1

Status: **scientific content frozen; author input only remains**

Freeze date: 2026-10-01

Source main commit: `944d6a06311360107f85f4c62c39ab5fd549ff40`

Machine-readable receipt:

`benchmarks/frozen/ttf_q_submission_scientific_freeze_v1.json`

## What is now frozen

The TTF-Q scientific package is complete and validated:

- exact response-blind information-decomposition theory;
- known-truth information benchmark;
- known-truth calibration/detectability benchmark;
- 653-species response-blind current/historical climate application;
- prospective 21-species Palearctic stopping application;
- five manuscript figures rendered only from frozen inputs;
- blinded and unblinded manuscripts;
- anonymous peer-review bundle machinery;
- submission-readiness checker;
- current-manuscript DOCX renderer smoke test.

The prospective Palearctic case remains closed as
`NOT_EVALUABLE_NONLEPIDOPTERA_PALEARCTIC_LGM_RELATION`.

No subgroup nucleotide identity, pairwise `T_st`, or `beta_LGM` was opened.

## Validation at freeze

The source-main state passed:

- **485 pytest tests**;
- five-figure rendering;
- submission-readiness structural audit;
- current blinded-manuscript and title-template DOCX rendering;
- continuous line-number and page-field checks.

The readiness audit reports:

- `scientific_package_ready = true`;
- structural failures = 0;
- abstract = 330 words;
- blinded manuscript = 6,404 words;
- expected figures = 5.

## What remains

Only author-controlled submission inputs remain.

### Hard blockers

1. Choose an open-source software license.
2. Complete the title page:
   - final author list/order;
   - affiliations and full addresses;
   - corresponding-author address and email;
   - acknowledgements/funding;
   - author contributions;
   - conflict-of-interest declaration.

### Author action, not a scientific blocker

Send or decline the prepared MEE presubmission scope enquiry.

### Non-blocking until public release

- final `CITATION.cff`;
- permanent archive DOI.

## Freeze rule

The following do **not** require scientific unfreezing:

- adding the explicitly chosen license;
- filling title-page identity/declaration fields;
- sending editorial correspondence;
- final DOCX/PDF rendering;
- technical packaging repairs with no scientific effect;
- release metadata after publication/archive decisions exist.

The following **do** require an explicit scientific-unfreeze decision:

- another analysis;
- another subgroup;
- another predictor;
- changing a threshold;
- opening a previously closed outcome;
- adding a new substantive ecological/genetic claim;
- materially reinterpreting the current conclusions.

The default next action is therefore **submission completion, not more science**.
