# Direct INSDC confirmatory source plan v0.1

Status: **source preflight only**. No nucleotide identity or empirical genetic outcome is authorized by this document.

## Why this route exists

OSC authentication currently blocks the untouched phylogatR portal archive. Rather than reopening mirror hunting, this arm goes to an upstream public source: a fixed GenBank/INSDC release. The preferred route is GenBank-only because source features can carry `lat_lon`, organism/taxonomy, accession and COI-family annotation in the same archival record. GBIF accession-coordinate linkage is a contingency, not the default.

## Fixed source

Use **GenBank flat-file release 273.0** (2026-08-15; cutoff 2026-08-02) and do not mix in daily updates. This makes the source universe reproducible. GenBank is an INSDC member; ENA/DDBJ exchange the same archival sequence data and are not independent validation arms.

## Phase A — metadata-only enumeration

Parse only fields needed to establish eligibility:

- accession/version
- organism and taxonomy needed to establish Animalia and exclude Aves/Chiroptera
- source `lat_lon`
- COI-family feature annotation (`COI`, `CO1`, `COX1` aliases)
- release/provenance metadata

Do **not** parse or persist ORIGIN sequence letters during this phase. Exclude the frozen Decker 221 species before any identity opening. Collapse exact duplicate localities and apply the frozen locality/geometry criteria. If at least 150 species survive, freeze the panel, source digest, species-name-hash cap and split before any identity access.

## Phase B — mask-only admissibility

Only after Phase A is frozen, read sequence characters into an in-memory canonical-valid/noncanonical mask. Persist masks/hashes, never nucleotide letters. Apply the same frozen >=50% comparable-site edge admissibility semantics. No edge rewiring or marker switching.

## Phase C — source-specific synthetic qualification

Run the already-defined TTF genetic inference family on the exact direct-INSDC geometry, but issue a **new source-specific authorization and Gate-D receipt**. A PASS from the phylogatR portal arm cannot be inherited because ascertainment differs.

## Phase D — one empirical opening

Only a source-specific Phase-4 authorization may open exact nucleotide identity and compute the frozen p-distance response. The primary estimand remains `place_beyond_ibd`; total genetic transfer remains descriptive.

## GBIF contingency

GBIF ordinary occurrence search is capped at 100,000 records per query. Never use the first 100,000 records as a panel. Bulk occurrence downloads require a GBIF account. Therefore GBIF is used only if GenBank-only metadata geometry is insufficient and a deterministic, non-truncated accession-coordinate retrieval contract is frozen first.

## Stop rules

- GenBank-only >=150 eligible species: implement the direct-INSDC arm; do not add GBIF merely to increase sample size.
- GenBank-only <150: evaluate the frozen GBIF-linkage contingency without row-level outcome inspection.
- No deterministic >=150-species route: `NOT_EVALUABLE_SOURCE_FEASIBILITY` and stop.
- Never tune taxon partitions, locality thresholds, locus aliases or source versions after seeing genetic outcomes.
