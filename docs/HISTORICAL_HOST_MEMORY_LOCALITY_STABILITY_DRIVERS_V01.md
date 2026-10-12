# TTF historical host-memory: drivers of locality-deletion sign instability

**Analytical scope (2026-10-12):** post-outcome, exploratory, aggregate-only; the original one-shot 259-species confirmatory analysis remains `NOT_EVALUABLE_EMPIRICAL_RESPONSE`. Four COI-invariant species are **not** reclassified or removed. Nothing here is a replacement primary test or a basis for a new p-value.

## Reproducible source

Input A: `private/species_stability.csv` (259 original species), SHA-256 `b447b360edc3d73ae62fe4c668f2cadf9c13f94a0292d8aee5e96828c6724f23`.

Input B: `private/locality_deletion_sensitivity.csv` (6,374 original locality-deletion scenarios), SHA-256 `ba6f0f5f1ec9b9da4e697a0a271f0156ac2ad80d385da0a793bba8e200f97686`.

Both are derived from the previously opened 259-species one-shot data and saved privately in `TTF_hhm259_locality_stability_20261012.zip`. No original FASTA identities, sequences, or per-edge genetic distances are opened here.

Run:
```bash
python scripts/audit_locality_stability_drivers.py \
  --species-stability PRIVATE/species_stability.csv \
  --site-deletions PRIVATE/locality_deletion_sensitivity.csv \
  --output results/locality_drivers_aggregate_v0.1.json
python -m pytest -q tests/test_audit_drivers.py
```

## Core descriptive results

- Frozen original species: 259; of those, 255 had a nonzero empirical slope, four remained strictly zero, with no standardized slope.
- 6,374 attempted original locality deletions: 6,373 finite and one undefined due to insufficient endpoint-disjoint IBD training. For the 255 originally nonzero species, 6,308 site deletions were finite.
- **420 / 6,308 = 6.66%** of finite deletions strictly reversed coefficient sign (negative×positive). **129/255 = 50.6%** of those species had at least one sign-reversing deletion. Three nonzero baselines became exactly zero under one deletion and are *not* counted as strict reversals.
- Descriptive Spearman correlation with each species' *fraction* of deletions reversing sign: `abs(standardized beta)`: −0.614; number of original localities: −0.351; maximum leave-two-endpoint-out IBD prediction leverage: +0.314; max residual host predictor leverage: +0.271; frozen host-history unique information: +0.031. **No p-values or confidence intervals are calculated.**
- Baseline standardized effect size quartiles (smallest→largest): strict-reversing species 57/64, 40/64, 19/63, 13/64. Quartile counts are fixed by sorting observed 255 defined estimates, not an outcome-independent sampling design.
- Within the weakest half of abs-standardized baselines, 58 species in the *high IBD leverage* half had 210/1025 (20.49%) site-deletion sign reversals, versus 69 species in the *low leverage* half with 166/2353 (7.05%). The corresponding figures in the strongest half were 34/1184 (2.87%) and 10/1746 (0.57%). This is an **unadjusted descriptive four-cell comparison**, not an independent conditional estimate.
- Replacing exactly one species' baseline slope with any one of its well-defined deletion slopes changes the 259-species arithmetic mean from +0.001260429 to a range of **+0.000745261 to +0.002017705**. Thus the *descriptive* mean remains positive under any **single** original locality deletion, but this does not correct the earlier finding that the mean is unusually sensitive to the five highest positive species, nor does it rescue the confirmatory one-shot test.

## Interpretation

There are two distinct forms of fragility:

1. **Within-species sign fragility** is concentrated where the already estimated slope is close to zero, and it is descriptively greater where the endpoint-disjoint IBD fit has high prediction leverage or few localities.
2. **Across-species arithmetic mean concentration** is dominated by a few large positive unstandardized slopes, a different issue diagnosed by the earlier tail-mean analysis.

Host-history predictor *uniqueness* is an external-data identifiability measure, not a guarantee that empirical response regression has a stable sign. The 259-species source also has known post-mask taxon and graph selection biases. None of these results identifies an ancient insect–host association, ecological causation, a new confirmation, or fitness/adaptation.

### Scientific boundary

The original `NOT_EVALUABLE_EMPIRICAL_RESPONSE` is **unchanged**. Because empirical nucleotide identities were already opened, do not relabel a 255-species subset or a nonstandardized-response model as a preregistered confirmatory replication. A future genuinely confirmatory test needs a **separate fresh panel** and a response-definability gate frozen before any sequence outcome is observed.

### Reproduction and privacy

`results/locality_drivers_aggregate_v0.1.json` and `results/independent_verification.json` include only aggregate quantities. `private/` retains derived per-species coefficients and site-deletion records for the owner's reproducibility; never commit the private directory to a public GitHub repository. `scripts/` and `tests/` can be shared publicly with the aggregate receipt.
