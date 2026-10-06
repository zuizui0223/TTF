# Historical biotic memory v0.1

Status: **fresh response-blind successor; genetic response closed**

## Biological question

> **Did past fragmentation of larval host resources leave a genetic memory in the insects that depended on them?**

The study does not ask whether two Lepidoptera species use similar hosts today. That prospective relation was already tested and was null.

The new hypothesis is within species:

> **Present-day post-IBD spatial genetic differentiation should be stronger across locality pairs whose larval host resources were historically less connected.**

The strongest version compares biotic history against abiotic history:

> historical host-resource disconnection should predict present genetic differentiation beyond geographic distance, present host availability, and the insect's own palaeoclimatic opportunity.

## Why this is a new family

Previous TTF ecological tests asked whether a source species transfers better to a target species when the pair shares traits, host-resource geography, environmental niche, historical climate, or sampled geography.

This study changes the response structure. It asks whether the **historical geography of an interacting partner** predicts the focal species' own spatial genetic structure.

No coefficient or subgroup from an already-opened genetic panel is reused.

## Fresh panel

The exact phylogatR archive is reused only as a source container; the genetic response panel is species-disjoint from all prior design universes.

Response-blind exclusions include the original 250-species phylogatR Phase-1 panel, the 500-species conditional programme, the 1000-species relational fresh universe, the 339-species S1 trait-gradient design, and the 500-species S2 host-resource design.

Their exact union contains 2,492 species, including 1,518 Lepidoptera.

After those exclusions:

- 25,505 unused Lepidoptera COI species remain in archive metadata;
- 8,922 have at least 12 aligned records;
- 475 also have a valid LepTraits host-family annotation;
- 220 pass the frozen response-blind locality/graph requirements;
- **192 use exactly one valid LepTraits host family and form the primary feasibility panel.**

No nucleotide identity or genetic distance from these candidates has been opened.

## External interaction stage

Species-level larval host identities are acquired before any genetic response from a frozen capture of Global Biotic Interactions (GloBI) and, where available, the underlying cited interaction datasets.

The GloBI response is treated as an external ecological input, not as outcome data. Raw records and citations must be frozen before downstream filtering.

The initial 192-species panel is not guaranteed to survive species-level host resolution. The next action is a response-blind feasibility census only.

## Historical host stage

Two evidence layers are deliberately separated.

1. **Primary scalable reconstruction:** resolve host plants, estimate their present climatic niche from external occurrence data, and project the same host niche onto frozen palaeoclimate time slices. This produces a host-availability history independent of focal insect genetics.
2. **Palaeoecological validation where available:** use Neotoma pollen, plant macrofossil, or other plant records as an external check on reconstructed host presence/history. Fossil coverage cannot be used to cherry-pick a genetic-positive subgroup after response opening.

The historical metric is not frozen yet. It will be frozen only after response-blind host-identity and palaeo-data feasibility are known.

## Candidate genetic estimand

For focal species s and frozen genetic graph edge e=(i,j), define a response-blind historical host-connectivity score H_s,e from the host landscape.

The future primary genetic response would be endpoint-safe post-IBD edge turnover.

The intended biological test is:

`post-IBD turnover_s,e ~ historical_host_disconnection_s,e + present_host_connectivity + insect_palaeoclimate_connectivity + frozen sampling controls + species effects`.

Direction: **greater historical host disconnection predicts greater present genetic differentiation.**

This model is only a target. No metric, threshold, time slice, resistance transform, or genetic opening is authorized until the response-blind external-data stages are frozen.

## Claim boundary

A future positive result could support the statement that historical resource continuity helped shape the present phylogeography of herbivorous insects.

It would not prove direct co-migration of a particular host and insect population, nor that the current interaction persisted unchanged through the entire late Quaternary.

A null after a qualified test would mean the frozen historical host reconstruction does not explain present spatial genetic structure in the tested panel.

## Immediate next gate

Acquire and freeze species-level host records for the 192 fresh single-host-family candidates.

Pass/fail at this gate concerns ecological data availability only. Genetics remains closed.
