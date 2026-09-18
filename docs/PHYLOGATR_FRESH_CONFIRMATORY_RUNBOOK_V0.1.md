# Fresh phylogatR confirmatory runbook v0.1

Status: **response-blind execution handoff**. This runbook changes no estimator, threshold, species rule, split, graph rule, bandwidth, reference family, or biological interpretation. The canonical empirical opening state remains `benchmarks/frozen/genetic_empirical_opening_state_v0.3.json`.

## Firewall

Until a valid Phase-4 identity-opening authorization exists: do not inspect fresh A/C/G/T identity, compute fresh pairwise genetic distances, compute the fresh empirical TTF statistic, use Decker empirical outcomes as a tuning signal, change the 0.10 Type-I Wilson ceiling or 0.80 shared-A2 Wilson floor, or drop species/edges/localities to rescue a failed gate. `NOT_EVALUABLE` is never a biological negative.

## 0. Obtain exactly one untouched portal archive

Use the frozen portal query: no location constraint / global, Kingdom = Animalia only, and no manual lower-taxon include/exclude choices. If the global Animalia request exceeds the portal UI limit, STOP rather than selecting an arbitrary smaller subset. Do not inspect nucleotide identities after download.

Keep the untouched `.zip`, `.tar.gz`, or `.tgz`; this archive is the byte-level provenance object for both Phase 1 and Phase 2.

## 1. Phase 1 — response-blind geometry

```bash
python scripts/run_phylogatr_confirmatory_phase1_intake.py \
  --archive /path/to/phylogatr_results.zip \
  --output-dir frozen/fresh_phylogatr_phase1
```

Proceed only if `phase1_manifest.json` has status `FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY` and `phase1_intake_receipt.json` has status `PHASE1_FROZEN`. If the panel is below the frozen minimum, STOP. The Phase-1 receipt freezes the archive SHA-256 without interpreting nucleotide identity; never regenerate or recompress the archive afterward.

## 2. Phase 2 — validity mask only

Pass the **same untouched archive** directly; do not manually extract it.

```bash
python scripts/run_phylogatr_confirmatory_phase2_intake.py \
  --archive /path/to/phylogatr_results.zip \
  --phase1-dir frozen/fresh_phylogatr_phase1 \
  --output-dir frozen/fresh_phylogatr_phase2
```

The wrapper verifies the Phase-1 hashes, archive SHA-256/format/root, `cite.txt`, and `genes.txt`, re-extracts safely to a temporary directory, runs only the frozen canonical-valid/noncanonical mask freezer, writes `phase2_intake_receipt.json`, and deletes the extraction. Nucleotide identity is not persisted and pairwise differences are not computed.

Proceed only when `phase2_manifest.json` status is `PASS_TO_SYNTHETIC_GATE`; otherwise STOP with identity closed.

## 3. Formal Phase 3 Gate-D — sole cross-species qualification

```bash
python scripts/authorize_phylogatr_phase3_gate_d.py \
  --geometry frozen/fresh_phylogatr_phase2/phase2_geometry.csv \
  --phase1-manifest frozen/fresh_phylogatr_phase1/phase1_manifest.json \
  --phase2-manifest frozen/fresh_phylogatr_phase2/phase2_manifest.json \
  --phase3-rule docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json \
  --repo-root . \
  --output frozen/fresh_phylogatr_phase3/authorization.json
```

Then execute only the shard ranges produced by `plan_phylogatr_phase3_shards.py`, followed by all planned `run_phylogatr_phase3_reference_shard.py` jobs, `aggregate_phylogatr_phase3_references.py`, all planned `run_phylogatr_phase3_observed_shard.py` jobs, and `aggregate_phylogatr_phase3_qualification.py`.

The resulting `ttf_genetic_phylogatr_phase3_qualification_v0.1` receipt is the **only** formal Gate-D PASS/NOT_EVALUABLE decision. Retention=1.0 is never rerun by the fragility diagnostic. If formal Phase 3 is NOT_EVALUABLE, identity remains closed.

## 4. Fresh within-species self-detectability

On the same Phase-2 survivor geometry and inherited split, use the frozen rule `docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json` and the scripts `run_phylogatr_phase3_self_reference_shard.py`, `aggregate_phylogatr_phase3_self_references.py`, `run_phylogatr_phase3_self_evaluation_shard.py`, and `aggregate_phylogatr_phase3_self_qualification.py`. This qualification does not alter formal Gate-D and does not open identity.

## 5. Response-blind fragility diagnostic

First create the deterministic geometry plan with `plan_phylogatr_gate_d_fragility.py`. Retention fractions are frozen at `1.0, 0.875, 0.75, 0.625, 0.5`.

Then execute `plan_phylogatr_gate_d_fragility_synthetic.py`, planned diagnostic reference shards, `aggregate_phylogatr_gate_d_fragility_references.py`, planned diagnostic observed shards, and `aggregate_phylogatr_gate_d_fragility_curve.py`.

Retention=1.0 imports the sole formal Phase-3 receipt and is not rerun. Identical geometry fingerprints reuse the prior summary; structural-support failures are recorded, not repaired. No thinned level has formal Gate-D or Phase-4 authority, and diagnostic values cannot rescue NOT_EVALUABLE or change thresholds. Completion status is `DIAGNOSTIC_FRAGILITY_CURVE_COMPLETE`.

## 6. Phase 4 identity-opening authorization

Attempt only after formal Phase-3 PASS and completion of the self/fragility provenance chain. A recommended filename convention is:

```bash
python scripts/authorize_phylogatr_phase4_identity_opening.py \
  --geometry frozen/fresh_phylogatr_phase2/phase2_geometry.csv \
  --phase1-manifest frozen/fresh_phylogatr_phase1/phase1_manifest.json \
  --phase2-manifest frozen/fresh_phylogatr_phase2/phase2_manifest.json \
  --phase3-rule docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json \
  --phase3-authorization frozen/fresh_phylogatr_phase3/authorization.json \
  --references frozen/fresh_phylogatr_phase3/references.json \
  --qualification frozen/fresh_phylogatr_phase3/qualification.json \
  --self-rule docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json \
  --self-references frozen/fresh_phylogatr_phase3/self_references.json \
  --self-qualification frozen/fresh_phylogatr_phase3/self_qualification.json \
  --fragility-execution-rule docs/supporting/genetic_phylogatr_gate_d_fragility_execution_v0.1.json \
  --fragility-curve frozen/fresh_phylogatr_phase3/fragility_curve.json \
  --phase4-rule docs/supporting/genetic_phylogatr_phase4_response_rule_v0.1.json \
  --opening-state benchmarks/frozen/genetic_empirical_opening_state_v0.1.json \
  --repo-root . \
  --output frozen/fresh_phylogatr_phase4/identity_opening_authorization.json
```

The v0.1 opening-state file is intentionally the immutable Phase-4 authorization baseline required by the frozen runtime contract; v0.3 remains the current source-access governance state. The fragility curve is a procedural prerequisite only: its values do not enter the scientific PASS/FAIL decision.

Do not inspect nucleotide identity until the authorization receipt status is exactly `AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING`.

### Frozen recovery for the executed 2026-09-16 lineage

The executed projected Phase-1 manifest was generated response-blindly and its exact SHA-256 was frozen into the Phase-3 authorization before Phase 3 ran, but the ephemeral manifest bytes are no longer retained. Do **not** reconstruct that file from partial summaries. For this one exact lineage, use the repository-frozen recovery record `docs/supporting/genetic_phylogatr_phase1_provenance_recovery_v0.1.json`.

Replace only the `--phase1-manifest ...` argument above with:

```bash
  --use-frozen-phase1-provenance-recovery \
```

The recovery path is accepted only when the frozen Phase-1 manifest SHA-256, exact Phase-2 manifest SHA-256, exact Phase-3 authorization SHA-256, Phase-1 summary, source-integrity witness, and outcome firewalls all match. It cannot reconstruct the missing manifest, alter Gate-D, change thresholds/estimands, or authorize Phase 4 by itself.

## 7. One empirical test

After valid authorization only, run `scripts/run_phylogatr_phase4_empirical_test.py` with the exact authorized files. Do not change species, marker, graph, split, bandwidth, reference family, or thresholds. The primary estimand is `place_beyond_ibd`; total genetic transfer is descriptive and cannot rescue the primary result.

## Terminal interpretation

- qualified cross-species positive: evidence for a transferable place component within the frozen domain;
- qualified cross-species non-significant + separately qualified empirical self positive: evidence consistent with lineage-conditioned spatial structure within the tested domain;
- cross-species non-significant without qualified self support: `NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING`;
- any failed geometry, mask, Type-I, power, source-integrity, or provenance gate: `NOT_EVALUABLE`, never a biological null.

## Current executed lineage

The authenticated broader portal archive has already been projected response-blindly to the frozen Animalia universe. Phase 1 and Phase 2 completed; formal fresh Phase 3 passed; fresh self-detectability completed; and the response-blind fragility curve completed. The remaining pre-identity step is the exact Phase-4 authorization above. No fresh nucleotide identity, pairwise genetic distance, or empirical TTF statistic may be opened before that authorization receipt exists.
