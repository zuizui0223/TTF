# Methods-only repository boundary

This document governs the split of ecological paper material from TTF.

## Keep in TTF

Material is retained when its primary purpose is to define, qualify, audit or demonstrate the transferability method:

- estimands, nulls and held-out-unit inference;
- simulation families and calibration;
- geometry and observation-support gates;
- genetic Phase 1–4 frozen receipts;
- finite-family and stopping rules;
- explicit `NOT_EVALUABLE` / abstention examples needed by the methods manuscript;
- tests required to verify those claims.

## Move to paper-specific repositories

Material is removed from TTF when its primary scientific contribution is a biological result rather than a methods result. The current migration target is the butterfly specialization ecology paper in `zuizui0223/chocho`.

The removal set includes paper-specific:

- manuscript and submission files;
- ecological synthesis documents;
- acquisition, mapping, analysis and rendering scripts;
- exploratory/frozen result receipts used only by that ecological manuscript;
- CI workflows dedicated to the ecological paper;
- source modules and tests with no remaining methods dependency.

## Migration gate

Deletion from TTF is authorized only after the destination repository satisfies all of the following:

1. source commit provenance is recorded;
2. paper manuscript and claim-map files are present;
3. required code and frozen receipts are present;
4. a clean environment installs declared dependencies;
5. paper-specific tests pass without network access;
6. previously skipped geometry tests run with declared geospatial dependencies.

TTF history is not rewritten. The split is represented by ordinary commits so provenance remains inspectable.
