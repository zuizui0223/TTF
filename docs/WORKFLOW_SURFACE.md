# Workflow surface

TTF separates currently supported automation from historical execution provenance.

## Current Actions surface

After the methods-only split and workflow cleanup, the active workflow directory contains **35** current or deliberately rerunnable definitions.

The repository-wide CI gate is:

- `.github/workflows/tests.yml` — installs `.[test]` and runs the full offline pytest suite.

Targeted smoke/contract workflows and selected manual qualification/reproduction entrypoints remain active where they still provide a useful supported interface to the methods repository.

The former `methods-ci.yml` workflow is archived because it duplicated `tests.yml` and its push trigger targeted only the now-merged refactor branch.

## Historical provenance

**63** completed, superseded, diagnostic, transport-repair, and one-shot workflow definitions are stored under `provenance/workflows/`.

The move preserves each original Git blob SHA. Frozen authorization/result files may therefore continue to record the historical path `.github/workflows/<name>.yml`; those strings describe where the workflow lived when the execution was frozen and are intentionally not rewritten.

The complete classification and blob mapping is:

`provenance/workflow_surface_classification_v0.1.json`

## Admission rule

A workflow belongs in `.github/workflows/` only when it is one of:

1. a current CI or interface smoke gate;
2. a deliberately supported reproducibility/qualification entrypoint whose rerun semantics remain meaningful.

A completed one-shot diagnostic, superseded generation, recovery workflow, or transport repair belongs in provenance after closure.

This is repository-surface management only. Moving a workflow does not change any estimator, threshold, frozen rule, candidate family, empirical result, or prior decision.
