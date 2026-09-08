# v0.8 validity and deployment-adequacy outcome

## Frozen interpretation

v0.8 separates two questions that the earlier conjunctive qualification mixed:

1. **Gate V:** does the profiled-private inference control false positives on fresh external sampling geometries?
2. **Gate D:** is one exact intended deployment geometry sufficiently informative to detect the predeclared moderate shared signal while retaining that false-positive control?

The numerical criteria were frozen before the v0.8 external validity panels were declared. No v0.2--v0.7 failure is reclassified by the v0.8 outcome.

## Gate V — PASS

Result: `results/v08_external_validity_v0.1.json`.
Workflow run: `34178680246`.

Only zero-shared private cells were opened. No shared positive-control world was generated in Gate V.

### European birds

| private amplitude | rejection rate | Wilson 95% upper |
| ---: | ---: | ---: |
| 0.5 | 0.046 | 0.06808 |
| 1.0 | 0.054 | 0.07743 |
| 2.0 | 0.032 | 0.05134 |
| 3.0 | 0.050 | 0.07277 |

Maximum Wilson upper = **0.07743 <= 0.10**.

### North-American waterbirds / shorebirds

| private amplitude | rejection rate | Wilson 95% upper |
| ---: | ---: | ---: |
| 0.5 | 0.058 | 0.08206 |
| 1.0 | 0.044 | 0.06572 |
| 2.0 | 0.034 | 0.05377 |
| 3.0 | 0.028 | 0.04645 |

Maximum Wilson upper = **0.08206 <= 0.10**.

Therefore `gate_v_pass = true`. This supports the narrow method claim that the frozen profiled-private inference controls the predeclared strong-private false-positive stress across both fresh external geometries. It is not a power claim.

## Gate D — FAIL

Result: `results/v08_gate_d_rgfca_reserve_v0.1.json`.
Workflow run: `34179030658`.

The exact frozen RGFCA reserve contains 250 species x 20 records, with a fixed 125/125 species split, k=3 and bandwidth=500 km. Empirical colour outcomes remained unopened throughout Gate D.

### Private validity on the deployment geometry

| private amplitude | rejection rate | Wilson 95% upper |
| ---: | ---: | ---: |
| 0.5 | 0.030 | 0.04890 |
| 1.0 | 0.026 | 0.04397 |
| 2.0 | 0.034 | 0.05377 |
| 3.0 | 0.046 | 0.06808 |

Maximum Wilson upper = **0.06808 <= 0.10**: deployment type-I criterion PASS.

### Predeclared moderate shared signal

Shared A=2 rejection/power = **0.666** with Wilson 95% interval **[0.62353, 0.70594]**.

The required lower bound is 0.80, therefore deployment power criterion FAIL and `gate_d_pass = false`.

## Consequence

The result is not evidence that empirical flower-colour sharedness is absent. The empirical colour outcome was not opened. It says that the exact 250 x 20 reserve design is **non-qualifying for the predeclared moderate-sharedness estimand** under the now externally validated false-positive procedure.

This panel is frozen as a failed deployment design and must not be retuned or reused as confirmatory evidence for a successor design. Any successor must be a newly frozen, outcome-unopened design layer with a separate authorization and lineage.
