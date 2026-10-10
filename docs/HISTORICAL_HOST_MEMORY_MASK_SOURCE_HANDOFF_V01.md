# Historical-host memory: exact mask-only source handoff (2026-10-10)

## 2026-10-11 response-blind attrition robustness and provenance update

The post-mask information gate was **not weakened or waived**. Instead, a
separate, source-locked, deterministic extremal audit established that for
**every possible subset of 200–321** species from the original 321-species
confirmatory panel, the original information criteria already pass if every
retained species keeps all its frozen graph edges. At the maximally adverse
200-species subset, the minimum possible median host-specific unique fraction
is **0.147075945318** (required >=0.10), the minimum fraction with unique
fraction >=0.05 is **0.92** (required >=0.70), and the minimum fraction
satisfying condition number <=30 is **1.00** (required >=0.90).

Reason: only **16/321** original confirmatory species have host-specific
unique fraction <0.05 and **0/321** fail the condition-number limit. Taking
the lowest 200 species is the exact worst-case median for a 200-species
subset. The bounds improve as the subset cardinality increases.

The sharp combinatorial result is fully reproducible using
`scripts/audit_historical_host_memory_information_attrition_bound.py` and
`benchmarks/frozen/historical_host_memory_attrition_robustness_result_v0.1.json`.
CI run [#38062731257](https://github.com/zuizui0223/TTF/actions/runs/38062731257)
passed five unit tests and the exact frozen input audit.

The canonical-mask source receipt now cryptographically carries forward the
frozen candidate, original roles, localities, edges, pre-mask synthetic PASS,
and exact frozen mask-rule Git blob. Both downstream survivor qualification
scripts reject substituted source bindings (CI mask #38062525294 13 tests;
survivor #38062510294 29 tests).

**Still unresolved:** at least 200 masks must actually pass, with every
original edge; the new independent post-mask synthetic test must also PASS.
This proof certifies **no empirical mask survival rate** and **no empirical
genetic host-memory effect**. The exact untouched phylogatR archive remains
absent from accessible connected sources.


## Completed and frozen

The original 642-species historical host-memory predictor has 641 externally complete
species (320 development and 321 confirmatory). The exact 321-species synthetic
qualification passed under its predeclared simulation family: 26/500 null
rejections (5.2%; Wilson upper 7.51%) and 500/500 positive-world rejections
(100%; Wilson lower 99.24%). The exact immutable result is
`benchmarks/frozen/historical_host_memory_synthetic_qualification_result_v0.1.json`.

**This is not evidence for an empirical genetic host-memory effect.**

The mask-only source and survivor information programs have passed synthetic
fixture tests. The independent post-mask survivor simulation has also passed
synthetic fixture tests. These tests do not count real mask survivors.

## Why the next computation has not run

The only authorized source is the **untouched original** phylogatR ZIP:
- size: **274,988,692 bytes**
- SHA-256: `5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5`

The currently connected Library and Google Drive searches yielded only
`phylogatr_results_exact_20260918.sha256.txt`, not the ZIP itself.
A SHA text receipt cannot substitute for raw bytes. The geometry-only archive
already in GitHub contains no aligned sequence characters and **cannot**
produce canonical-site masks.

No nucleotide identities, empirical genetic distances, post-IBD turnover or
host-memory coefficients have been opened for this prospective 321-species
confirmatory panel.

## Exact sequence for authorized mask-only work

1. Recover an untouched local ZIP and check its SHA-256 and byte size.
2. Determine whether the original source may legally be redistributed or
   uploaded to a public repository. **Do not upload restricted or personal
   files to a public GitHub Actions artifact merely for convenience.**
   If redistribution is uncertain, keep the source local in a protected
   environment and use the Python commands there.
3. Confirm the frozen mask-only rules:
   `docs/supporting/historical_host_memory_character_mask_rule_v0.1.json`.
4. Verify ZIP members, extract into a temporary private directory, and open
   only Boolean A/C/G/T-validity masks for **321 confirmatory species**. The
   development genetic response remains permanently closed.
5. Retain species only when **every pre-existing graph edge** has at least one
   valid cross-locality sequence pair with jointly canonical columns
   >=ceil(0.5 × frozen alignment length). No backfill or rewiring. If fewer
   than **200** confirmatory species survive, stop.
6. Re-evaluate the original information thresholds on **exact survivors**.
7. Only if qualified, repeat **1,999 reference + 500 null + 500 positive**
   synthetic worlds with the new, already-frozen independent survivor
   namespace; keep all biological/effect/nuisance rules unchanged.
8. Regardless of PASS, **stop before nucleotide identity**. Further empirical
   opening requires a separate, explicit, exact-source authorization and
   honest interpretation of historical host plausibility.

## Execution paths

**Local/private:** Use
`scripts/verify_and_extract_historical_host_mask_source.py`,
`scripts/freeze_historical_host_memory_character_masks.py`,
`scripts/qualify_historical_host_memory_survivor_information.py`, and
`scripts/qualify_historical_host_memory_survivor_synthetic.py`
with the pinned source CSVs and prior frozen artifact inputs. The manual
Actions workflow contains these precise command arguments and is a reference
for the steps even when run privately. Never persist extracted FASTA files in
Git or report artifacts.

**GitHub Actions (only if redistribution rights and access are appropriate):**
`.github/workflows/historical-host-memory-manual-mask-intake-v01.yml` is
**workflow_dispatch only**, requiring a provider artifact ID that contains
exactly one untouched original phylogatR ZIP. It does not run on pushes,
never auto-opens a nucleotide identity, and refuses source size/SHA drift.
It uploads only source provenance and mask/survivor/synthetic results, never
the extracted original FASTA. Be aware that the input artifact itself may be
visible to people with access to the repository.

## Current status

- Source-identity repair and 73-species parallel-panel overlap: checked
- Current/LGM predictor information qualification: PASS
- Initial independent-null/positive synthetic qualification: PASS
- Mask-readiness and survivor synthetic fixture tests: PASS
- Actual confirmatory mask survivors: **UNKNOWN**
- Exact survivor information and survivor synthetic qualification: **NOT RUN**
- Empirical host-memory hypothesis: **NOT TESTED**
