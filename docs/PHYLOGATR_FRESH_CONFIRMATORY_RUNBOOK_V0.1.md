# Fresh phylogatR confirmatory runbook v0.1

This is the execution handoff for the fresh genetic TTF confirmatory arm. It does not change any estimator, threshold, species rule, split, graph rule, bandwidth, reference family, or biological interpretation. It exists to remove analyst-choice surfaces between the already-frozen stages.

## Non-negotiable firewall

Until a valid Phase-4 identity-opening authorization exists:

- do not inspect fresh A/C/G/T nucleotide identity;
- do not compute fresh pairwise genetic distances;
- do not compute the fresh empirical TTF statistic;
- do not use Decker empirical outcomes as a tuning signal;
- do not change the 0.10 Type-I Wilson ceiling or 0.80 shared-A2 Wilson floor;
- do not drop species, edges, localities, or alter the inherited split to rescue a failed gate.

`NOT_EVALUABLE` is never a biological negative.

## 0. Portal acquisition

Use the already-frozen portal query only:

- location: none / global;
- kingdom: Animalia;
- do not manually include or exclude lower taxa;
- do not inspect downloaded nucleotide identities.

If the global Animalia query exceeds the portal UI download limit, STOP. Do not choose a smaller taxon or geographic subset ad hoc. A response-blind partition/union rule must be frozen first.

Keep the untouched `.zip`, `.tar.gz`, or `.tgz`. The archive itself is the byte-level provenance object across Phases 1 and 2.

## 1. Phase 1 — response-blind geometry only

Preferred archive-direct command:

```bash
python scripts/run_phylogatr_confirmatory_phase1_intake.py \
  --archive /path/to/phylogatr_results.zip \
  --output-dir frozen/fresh_phylogatr_phase1
```

Expected files:

- `phase1_geometry.csv`
- `phase1_manifest.json`
- `phase1_intake_receipt.json`

Proceed only when both conditions hold:

- manifest status = `FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY`;
- intake receipt status = `PHASE1_FROZEN`.

If the panel is below the frozen minimum, STOP. Do not open Phase-2 masks.

For archive input, the Phase-1 receipt freezes the full archive SHA-256 without interpreting nucleotide identity. Do not replace, recompress, rename-and-repack, or otherwise regenerate the archive after this point.

## 2. Phase 2 — character-validity mask only

Use the **same untouched archive** from Phase 1. Do not manually extract it.

```bash
python scripts/run_phylogatr_confirmatory_phase2_intake.py \
  --archive /path/to/phylogatr_results.zip \
  --phase1-dir frozen/fresh_phylogatr_phase1 \
  --output-dir frozen/fresh_phylogatr_phase2
```

The wrapper must verify the Phase-1 manifest/geometry hashes, archive SHA-256, archive format, discovered archive root, `cite.txt`, and `genes.txt`, then safely re-extract to a temporary directory. The temporary extraction is deleted after Phase 2.

Expected files:

- `phase2_geometry.csv`
- `phase2_manifest.json`
- `phase2_intake_receipt.json`

Allowed character-level information in Phase 2 is only canonical-valid versus noncanonical/missing. Nucleotide identity is not persisted and pairwise differences are not computed.

Proceed only when manifest status = `PASS_TO_SYNTHETIC_GATE`. Otherwise STOP with nucleotide identity closed.

## 3. Formal Phase 3 Gate-D — the only cross-species qualification decision

Authorize formal Phase 3 on the exact Phase-2 survivor geometry:

```bash
python scripts/authorize_phylogatr_phase3_gate_d.py \
  --geometry frozen/fresh_phylogatr_phase2/phase2_geometry.csv \
  --phase1-manifest frozen/fresh_phylogatr_phase1/phase1_manifest.json \
  --phase2-manifest frozen/fresh_phylogatr_phase2/phase2_manifest.json \
  --phase3-rule docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json \
  --repo-root . \
  --output frozen/fresh_phylogatr_phase3/authorization.json
```

Then use the frozen Phase-3 execution scripts in this order:

1. `plan_phylogatr_phase3_shards.py`
2. all planned `run_phylogatr_phase3_reference_shard.py` jobs
3. `aggregate_phylogatr_phase3_references.py`
4. all planned `run_phylogatr_phase3_observed_shard.py` jobs
5. `aggregate_phylogatr_phase3_qualification.py`

The planner receipt is authoritative for shard ranges. Do not invent, omit, or duplicate replicate ranges manually.

The resulting `ttf_genetic_phylogatr_phase3_qualification_v0.1` receipt is the **only** formal Gate-D PASS / NOT_EVALUABLE decision. Retention=1.0 is never rerun by the fragility diagnostic.

If formal Phase 3 is NOT_EVALUABLE, fresh nucleotide identity remains closed. The response-blind fragility curve may still be completed for method diagnosis, but it cannot rescue the formal result.

## 4. Fresh within-species self-detectability qualification

Before empirical interpretation can distinguish a cross-species negative from lineage-conditioned spatial structure, complete the separately frozen self-detectability qualification on the same Phase-2 survivor geometry and inherited split.

Use:

- `docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json`;
- `run_phylogatr_phase3_self_reference_shard.py`;
- `aggregate_phylogatr_phase3_self_references.py`;
- `run_phylogatr_phase3_self_evaluation_shard.py`;
- `aggregate_phylogatr_phase3_self_qualification.py`.

This self qualification does not alter formal Gate-D and does not open nucleotide identity.

## 5. Response-blind Gate-D fragility diagnostic

The diagnostic is mandatory to complete before any Phase-4 identity opening, but its **values are not a second scientific gate**.

First create the deterministic nested geometry plan from the exact Phase-2 survivor geometry using `plan_phylogatr_gate_d_fragility.py`. Retention levels are frozen at:

`1.0, 0.875, 0.75, 0.625, 0.5`.

Then create and execute the synthetic diagnostic plan with:

1. `plan_phylogatr_gate_d_fragility_synthetic.py`
2. planned `run_phylogatr_gate_d_fragility_reference_shard.py` jobs
3. `aggregate_phylogatr_gate_d_fragility_references.py`
4. planned `run_phylogatr_gate_d_fragility_observed_shard.py` jobs
5. `aggregate_phylogatr_gate_d_fragility_curve.py`

Rules:

- retention=1.0 imports the one formal Phase-3 qualification receipt and is not rerun;
- identical geometry fingerprints reuse the already-computed summary rather than adding Monte Carlo noise;
- a thinned level below formal structural support is recorded as `STRUCTURAL_SUPPORT_BELOW_FORMAL_MINIMUM`, not repaired;
- no thinned level has formal Gate-D decision authority;
- no thinned level has Phase-4 opening authority;
- no diagnostic value can rescue NOT_EVALUABLE or change thresholds.

A completed curve has status `DIAGNOSTIC_FRAGILITY_CURVE_COMPLETE`.

## 6. Phase 4 identity-opening authorization

Only attempt this stage when formal Phase 3 has PASS status. The authorization script additionally requires the completed fragility curve and completed self-detectability qualification/provenance chain.

Use one exact provenance chain. A recommended file convention is shown below; the chosen Phase-3/self/fragility output filenames must be the exact files produced and frozen by the preceding stages.

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

`genetic_empirical_opening_state_v0.1.json` is intentionally the immutable Phase-4 authorization baseline used by the frozen runtime contract. The current source-access governance state remains v0.3; do not substitute v0.3 for the Phase-4 baseline unless the authorization/runtime contract is prospectively revised and requalified before any identity opening.

The fragility curve is a **procedural prerequisite only**. Its diagnostic values do not enter the scientific PASS/FAIL decision; only formal Phase-3 PASS does.

Do not inspect nucleotide identity until the authorization receipt status is exactly:

`AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING`.

## 7. Phase 4 empirical result

After valid authorization only, run `scripts/run_phylogatr_phase4_empirical_test.py` with the exact authorized files. Do not change species, marker, graph, split, bandwidth, reference family, or any threshold.

The primary estimand is `place_beyond_ibd`. The total-genetic-distance transfer surface remains descriptive and cannot rescue the primary result.

## Terminal interpretation branches

- formal cross-species positive: evidence for a transferable place component of within-species genetic differentiation within the frozen domain;
- formal cross-species negative + separately qualified empirical self positive: evidence consistent with lineage-conditioned spatial structure within the tested domain;
- formal cross-species negative without qualified self support: `NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING`;
- any failed geometry, mask, Type-I, power, source-integrity, or provenance gate: `NOT_EVALUABLE`, never a biological null.

## Current external blocker

Until an untouched authenticated phylogatR portal archive is supplied, stop before Step 1. No further same-purpose source hunting, public substitute search, threshold tuning, or empirical opening is authorized.
