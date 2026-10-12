# TTF historical host-memory: geographic sampling context for unstable locality deletions

**2026-10-12 — exploratory and response-exposed.** The one-shot 259-species primary remains `NOT_EVALUABLE_EMPIRICAL_RESPONSE` because four original COI species have exactly zero standardized response variance. No panel alteration, new p-value, genetic response opening, or LGM host-occurrence inference is made by this audit.

## Frozen source integration

The previously computed locality-deletion signs (`6,374` attempts) and per-species original slope estimates were joined to the original *response-blind* 17,202-locality and 71,728-edge genetic-sampling geometry by species and exact locality indices.

All four inputs must pass their original byte-exact SHA-256 checks (see `scripts/audit_locality_geography_drivers.py`). For each of the 259 original mask-surviving species, the program recomputes nearest-other-sample chord distance, distance to its sampled ECEF centroid, median distance to other sampled localities, and degree in the original frozen graph. These physical/topological values are converted to within-species average-tie rank / n (not compared across species in kilometers).

The script produces only aggregate values. It never opens FASTA sequences, source nucleotide identities or per-edge pairwise genetic distances. The original 259 species and all their original edges remain unchanged.

## Descriptive results

Among 255 nonzero-baseline species there were `6,308` defined one-locality deletions and `420` **strictly opposite-sign** slopes. Three deletions yielded exactly zero and are not counted as sign reversals. The four originally COI-invariant species remain outside the sign-comparison denominator without being dropped from the study.

| Relative locality position, within its species | Mean rank for 420 flip-causing deletions | Mean rank for 5,888 other defined deletions | Species with flip sites ranked higher vs other sites |
|---|---:|---:|---:|
| Nearest-neighbor distance | 0.5722 | 0.5165 | **82 / 129** |
| Median distance to all other sites | 0.5584 | 0.5175 | 71 / 129 (2 ties) |
| Distance to sample centroid | 0.5395 | 0.5189 | 70 / 129 (1 tie) |
| Graph-edge degree | 0.5202 | 0.5201 | 60 / 129 (1 tie) |

The within-species paired median rank difference for nearest-neighbor separation is **+0.1062** in those 129 species with at least one sign reversal. There is no corresponding overall tendency for flip-causing sites to have high graph degree.

Nearest-neighbor ECEF *chord* distance is monotonically related to great-circle distance for any two terrestrial sites, but the ranking measures **observed sample isolation**, not a dated dispersal barrier or paleoclimate refuge. For other distance metrics, ECEF geometry is only a descriptive approximation.

The differences are mild and observational, may be driven by graph sampling or incomplete geographic coverage, and are not an independent causal test. Each species contributes multiple correlated deletions; do not treat 6,308 deletion scenarios as independent biological replicates. No p-values, confidence intervals, or selection rules have been fitted.

## Independent verification

A separate implementation using pandas `read_csv(float_precision="round_trip")` and scipy `rankdata` reproduced the four summary means, paired medians and species counts. Default pandas floating conversion changes one near-tied nearest-neighbor rank; the round-trip float reader matches Python's original CSV floats and the source-locked primary audit.

The public result is `benchmarks/frozen/historical_host_memory_postoutcome_geographic_locality_drivers_v0.1.json`. The script and tests are shared without species-level response data. Owner-only derived deletion rows remain in the earlier private ZIP.

## Relation to the ecological research

Historical larval host-resource opportunity is a climate-suitability proxy, not demonstrated ancient host occupancy. The observed tendency of sign-flipping sites to be geographically more isolated is compatible with **sampling-edge leverage** in mitochondrial population-structure measurements. It cannot separate dispersal barriers, range histories, sequence sampling, or IBD correction as causal mechanisms. Nor does it rescue the original not-evaluable 259-species preregistered result.

A future genuine confirmatory ecological test needs a fresh, independent panel, frozen response-definability handling, and a separate evaluation of sampling geometry and host-resource history **before** any genetic effect is seen.
