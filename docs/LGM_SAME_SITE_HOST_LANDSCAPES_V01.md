# Same recorded places, different host landscapes: LGM follow-up v0.1

Date: 2026-10-19. Base: `7319fc9b555aabb2bd6e543dec90ee632588e33c` (PR #143).

## Central question

**Does the LGM-to-0BP change in potential host-climate connectivity alter which
geographic connections are relatively costly, rather than merely their average
cost, when insect species are compared on exactly the same recorded endpoints?**

The eventual biological hypothesis is that a place becomes an effective barrier
through a species' ecological resource requirements and history, not geographic
position alone. This analysis does not test that genetic mechanism. It develops
a controlled ecological comparison before acquiring genuinely independent
population-genetic validation.

This is a post-exposure exploratory successor to Issue #141, not a fresh
confirmation, rescue of the closed 259-species test, or reopening of the
Palearctic non-Lepidoptera study. The study region is the previous compact
Bungalotis pilot region, not a new global or Palearctic sample.

## What was actually executed

Three source files were recovered byte-for-byte and verified against the live
repository's bindings: the 642-species candidate table, 17,202 locality rows, and
the original 64-edge Bungalotis CHELSA LCP output. No nucleotide file, genetic
distance, genetic coefficient, or new climate raster outcome was accessed.

1. Reanalyse the existing pilot descriptively, separating **magnitude** from
   **relative edge ordering**. Keep the same-pixel zero in original aggregate
   summaries and report it separately for positive-path rank/ratio summaries.
2. Screen every one of the 642 candidates, using only recorded coordinates and
   already-recorded host identities, in the original unrounded coordinate
   rectangle plus 0.12 degrees. Bounds: west -85.63, south 10.335, east -85.1604,
   north 11.17. This rectangle lies inside, but is not identical to, the
   outward-snapped 101 x 57 native raster window.
3. Retain every species with at least 12 distinct recorded coordinates inside
   that rectangle. Enumerate all unordered species pairs, including pairs with
   no coordinate overlap. Compare coordinates numerically without rounding or
   a geographic tolerance, and normalize accepted plant IDs exactly.
4. Export matched endpoint queries for every species pair sharing at least 12
   recorded coordinates. This is an explicitly exploratory engineering tier,
   not a statistical qualification threshold; counts at 1, 3, 6 and 12 are all
   reported and all 666 pairs are retained. Each pair of shared coordinates
   retains both insects' original locality indices.

## Executed pilot diagnostic: a mean gap is not a different barrier map

On the 63 positive-path edges (one of the original 64 pairs shares a raster
pixel), the Spearman correlation of host versus insect-self LGM-minus-0BP
cost changes is **0.957417**. After dividing both changes by the original
geographic chord distance, it is still **0.777170**.

The median LGM/0BP cost ratio is **1.170306** for the host analogue and
**3.256513** for the insect-self analogue. Thus the reference landscapes differ
in magnitude while retaining substantial common ordering of expensive
connections. Neither a mean gap nor 64 rarefaction-draw means establishes that
historical barriers occur in different places.

These correlations are descriptive. There are no p-values or confidence
intervals, edges share nodes, and dividing by distance does not remove every
spatial confound. They use the original **665-vs-26** reference clouds. The
published 26-point rarefaction receipt retains draw-level means, not edgewise
cost matrices, so equal-count **spatial-rank** robustness has not been checked.
This is not evidence that the host-history hypothesis is false.

## Executed geometry census

| Design quantity | Result |
|---|---:|
| Candidate species screened | 642 |
| Locality rows screened | 17,202 |
| Species with >=12 coordinates in the rectangle | 37 |
| Families, using the frozen taxonomy | 10 |
| Regional species-coordinate records | 957 |
| Distinct regional recorded coordinate pairs | 377 |
| All unordered regional species pairs | 666 |
| Pairs sharing >=1 / >=3 / >=6 / >=12 coordinates | 546 / 266 / 78 / 4 |
| Unique species in the >=12 shared-coordinate tier | 6 |
| Newly enumerated matched coordinate-pair queries in that tier | 343 |

The four high-overlap pairs all have disjoint **recorded accepted-host-species**
sets, not necessarily disjoint host families or disjoint climatic niches:

| Insect A | Insect B | Shared recorded coordinates | Matched coordinate-pair queries |
|---|---|---:|---:|
| Cautethia spuria | Crinodes besckei | 16 | 120 |
| Crinodes besckei | Xylophanes pluto | 14 | 91 |
| Isognathus rimosa | Xylophanes pluto | 12 | 66 |
| Quinta cannae | Xenophanes tryxus | 12 | 66 |

The six species represent three frozen families. The 343 query rows are not
343 independent biological replicates or 343 validated positive-length raster
paths. They are a new exploratory all-pairs query geometry, not a modification
of any original frozen genetic k-neighbor graph.

Two species pairs have identical recorded host sets and shared coordinates:
Cautethia spuria / Xylophanes pluto (10 shared coordinates, Chiococca alba) and
Conchylodes ovulalis / Sabulodes arge (3, Parthenium hysterophorus). Their full
host-reference cost surfaces must be identical under the same source, mask,
grid and algorithm. That is an algebraic implementation control, **not evidence
of a biological host effect**. Their 48 matched queries are exported separately.

## Next ecological computation, specified but NOT executed here

For the matched endpoints, compute host and insect-self climate cost at both
LGM and 0BP on the same native grid. Preserve reference taxonomy and retain a
shared, explicit comparison domain. Separate the regional record cloud from a
species' global climatic niche. Recheck finite climate, source equivalence,
land/ice assumptions, positive path lengths, and ROI boundary sensitivity.
Do not snap invalid endpoints to a nearby favorable pixel.

The necessary outputs are **per species x matched edge x time x reference
scenario**, including all equal-count sensitivity draws, rather than only
species means. This permits separate comparisons of:

- average historical cost change (magnitude);
- the ordering and spatial allocation of historical cost change (shape);
- which part of host change remains beyond insect-self climate and common
  route geometry.

A large magnitude difference with nearly unchanged edge ordering cannot be
promoted to a claim of different barrier locations. Conversely, robust shape
changes would motivate, but would not themselves establish, a genetic
mechanism. Use joint spatial/locality-aware validation for any later inference;
ordinary edge-wise regression p-values would pseudoreplicate sites.

For identical-host implementation controls, use the full identical reference
or the same deterministic subset, not independent random subsets whose
sampling variation would manufacture a difference.

CHELSA-TraCE21k provides climatic reconstruction, not verified ancient host
occupancy. Its original model description is Karger et al. (2023),
https://doi.org/10.5194/cp-19-439-2023. Two time slices do not establish continuous
connectivity through the intervening 21,000 years.

## Boundaries that remain binding

Coordinate equality can reflect rounding or shared locality labels; coordinate
uncertainty, collection timing and microhabitat equivalence have not been
validated. Neither host lists nor present native-range unions establish actual
host use at these particular sites or historical occupancy. The 37-species
subset is geographically selected after previous ecological outcome exposure;
it is not a source-disjoint genetic panel. The original 259-species state
remains `NOT_EVALUABLE_EMPIRICAL_RESPONSE` with no replacement p-value.

## Reproduction

Run from a TTF checkout containing the new script and with the three exact
source CSVs staged in a directory:

```bash
python scripts/analyze_historical_host_same_site_design.py \
  --sources-directory /path/to/source-csvs --output-directory /path/to/results
python -m pytest -q tests/test_same_site_lgm.py
```

Outputs: the complete 642-species screen, the 37-species panel, all 666 species
pairs, 343 matched high-overlap queries, 48 identical-host control queries, and
`same_site_design_summary.json`. SHA-256 checks prevent substituting a different
candidate universe, coordinate table or pilot cost table.

Verification in this session: **12 new focused tests passed**, exact-source
replay completed, and separate pandas joins/ranks checked the main counts,
matched coordinate identities and correlations. The full historical TTF test
suite and new multispecies raster runs were not executed.
