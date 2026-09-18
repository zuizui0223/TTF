# Genetic TTF result-to-manuscript export

Current terminal state: the exporter has been run on the frozen fresh Phase-4 receipt chain (103 training / 108 evaluation species; Branch B), and the validated scalar prose is incorporated into `manuscript/genetic_ttf_flagship_v0.2.md`.

The exporter converts already-completed scalar receipts into a manuscript Results passage and a complete evaluation-species table. It does not read FASTA/occurrence files, fit a model, recompute a test, qualify a geometry, or authorize identity opening.

## Empirical report

After the authorized Phase-4 runner has completed:

```bash
PYTHONPATH=src python scripts/export_genetic_ttf_manuscript.py \
  --result /path/to/phase4_empirical_result.json \
  --authorization /path/to/phase4_authorization.json \
  --qualification /path/to/phase3_qualification.json \
  --self-qualification /path/to/phase3_self_qualification.json \
  --output-dir /path/to/new_report
```

The default rule paths are the checked-in Phase-4 response and fresh self-detectability rules. To read an archived version, supply `--phase4-rule` and `--self-rule` explicitly; their hashes must match the authorization.

Outputs:

| File | Content |
| --- | --- |
| `results.md` | Qualification bounds, full species counts, primary estimate and p-value, the appropriate interpretation branch, self support, and a separately labeled total-transfer descriptor. |
| `species_scores.csv` | Every frozen evaluation species with primary, self and descriptive total scores. Unavailable total scores are blank and explicitly flagged. |
| `export_manifest.json` | Input receipt SHA-256 values, geometry fingerprint, decision, counts and output hashes. |

The output directory must be new. Validation happens before any report files are written; files are staged together and then moved into place. Inputs are read-only.

## Closed or non-evaluable route

Historical or failed routes can still be reported without empirical data:

```bash
PYTHONPATH=src python scripts/export_genetic_ttf_manuscript.py \
  --closed-receipt benchmarks/frozen/genetic_empirical_opening_state_v0.3.json \
  --output-dir /path/to/new_status_report
```

This writes **only** `status.md` and `export_manifest.json`. It produces no empirical Results section or species scores. Closed mode also accepts the terminal panel-too-small Phase-1/2 schemas and a failed fresh Phase-3 qualification with the empirical firewall closed. It does not interpret an arbitrary pipeline-progress JSON or synthetic PASS as an empirical result.

## Validation boundary

The report checks consistency of supplied frozen receipts, not their external authenticity and not the complete statistical analysis from raw sources.

- The runner now includes `phase4_authorization_sha256` in its result. The exporter requires this field; legacy results without this link are unsupported.
- The supplied qualification, self qualification and two rule files must match the authorization hashes.
- Geometry fingerprints, species split and source-integrity flags must agree. Every frozen evaluation species must appear in all three summaries.
- The cross-species gate must be PASS with consistent Type-I/power flags and numerical bounds under the frozen 0.10/0.80 limits.
- Primary and self aggregates must equal the means of their complete species tables. p-values, alpha, component maximum, boolean significance flags, self qualification and the final decision must agree.
- The descriptive total has no p-value, significance label or qualification. Its calculation metadata must match the declared pre-IBD response with retained geometric length control. Any unavailable species makes the full total aggregate unavailable, with no reduced-panel average.
- Duplicate JSON keys, non-finite numbers, missing species and contradictory receipts are rejected.
- The Results text is selected from validated numeric/boolean fields. Free-form `interpretation` text in the input is not copied.

A significant transfer result supports predictive geographic information within the qualified domain. It does not causally distinguish place from shared demographic history. A qualified non-significant result with qualified significant self support is described as **consistent with** lineage conditioning, not proof of zero transfer. The self statistic uses a frozen upper-tail structural-null reference; its raw numerical sign is not interpreted against zero.

The runner remains covered by existing Phase-4 authorization code hashes. The exporter is downstream reporting software: it cannot reopen or rerun the empirical analysis. This change does not revise any frozen gate, null family or source-access restriction.

## Development evidence

The tests use fabricated numeric reporting fixtures plus frozen receipts. They exercise all decision branches, alpha-boundary behavior, complete species tables, unavailable descriptors, mismatched receipt hashes, incorrect species means, forged pass labels with failing bounds, source-integrity drift, JSON ambiguity, CLI behavior, input immutability and overwrite refusal. The terminal empirical export was generated only after the one-shot Phase-4 result had been frozen.
