# Historical-host parallel panels: source identity and overlap firewall (2026-10-08)

## What was measured — not a genetic outcome

Two **response-blind, frozen** candidate universes were compared using original, persistent Library CSV bytes, not a GitHub spreadsheet/text-preview reconstruction:

| Protocol | Candidate species | Original SHA-256 |
| --- | ---: | --- |
| historical-host-memory, within-species CHELSA climate analog (642 candidates) | 642 | `c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b` |
| historical-host-connectivity, host resistance and LGM resources (panel v0.2) | 135 | `a402483497bc7bdf9b7427bfc859ce9dcc9183eebf5964d4591f352847130484` |

Both CSVs contain exactly the declared number of unique `species` strings. The **species-set intersection has 73 species** (54.1% of the connectivity panel; 11.4% of the memory panel). The species are listed in `benchmarks/frozen/historical_host_parallel_panel_input_overlap_audit_v0.1.json`, which is a descriptive source-identity receipt only. There are 704 distinct species across the union, not 777 independent species.

**Inferential implication:** These are not two fully independent species-disjoint replications. Even if each prospective estimator later passes its gates, reporting them as two independent genetic samples without dependence accounting would overstate replication. The biological estimands differ (historical host-resource climatic compatibility versus spatial host resistance) and must not be pooled by label alone. The original panel memberships, future qualifications, and no-backfill rules remain unchanged.

## Unresolved candidate CSV transport

During a GitHub Actions overlap preflight, `benchmarks/frozen/historical_host_memory_candidate_642_v0.1.csv` **as checked out from GitHub** hashed to:

`81fe699448774ebd9a31de0076dfc5d9b2970dc8c09ab40efa62554bdb83e1b1`

rather than the original frozen `c36c...`. The first overlap workflow correctly failed closed. Neither an apparently matching set of species names nor row counts can resolve a byte-identical scientific input requirement. The separate debug workflow has not yet established the exact transport transformation.

**HOLD:** No historical-host-memory predictor should run using the GitHub copy until its exact candidate-source hash is restored to the frozen original, or a fully audited, deterministic transport repair restores those exact original bytes. Never quietly rewrite the frozen expected hash to the repository copy. Identical rule applies to the localities and edges sources.

The new `verify_frozen_input_hashes` and `verify_host_pair_identity` checks in `scripts/build_historical_host_memory_predictor.py` enforce these pre-response source conditions. Scoped unit tests passed (14/14 in GitHub Actions run 37783074401).

## Biological next gate

The confirmed cross-lineage codistribution null and historical-host designs do **not** establish that past host distributions caused present genetic breaks. The testable mechanism remains: does host-specific LGM resource opportunity explain present post-IBD differentiation after the insect's own abiotic climate history and current host mismatch are controlled?

- Preserve both distinct frozen questions and all stopping/qualification rules.
- Confirm and restore exact external candidate/geometry/host-pair identities before computing historical predictor values.
- Require temporal host plausibility checks: a modern cultivated taxon/interaction cannot be projected as an actual 21 ka host interaction by assumption.
- Do not interpret the Neotoma technical timeout as palaeoecological absence.
- No confirmation genetic nucleotide identities, pairwise distances, turnover, or historical-host coefficients have been opened for either prospective historical-host panel.

This note changes **neither** the 211-species main analysis nor the 326-species codistribution result; both remain frozen.
