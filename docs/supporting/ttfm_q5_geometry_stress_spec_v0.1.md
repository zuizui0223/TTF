# TTF-M / TTF-C Q5 heterogeneous-geometry stress v0.1

Status: **prospective design only; not yet executed**.

The idealized qualification established that TTF-M and TTF-C can control false positives and detect their intended shared fronts under a balanced 40-system × 60-record geometry. Q5 asks whether those operating properties survive heterogeneous sampling geometries closer to opportunistic ecological data.

This is a new qualification layer. It does not retroactively modify the idealized PASS.

## Frozen stress dimensions

Each synthetic world draws system-specific observation geometry before states are generated. The geometry draw is independent of the response states.

1. **Unequal sample size**
   - `n_s` sampled from `{30, 45, 60, 90, 120}` with all five levels represented in train and evaluation systems.

2. **Spatial clustering**
   - balanced uniform sampling;
   - moderate clustered sampling;
   - strong clustered sampling with sparse boundary coverage.

3. **Domain truncation**
   - full domain;
   - one-sided domain loss;
   - interior gap removing observations near part of the transition front.

4. **System-specific missingness**
   - 0%, 10%, and 25% response-blind observation loss applied after geometry generation.

5. **Train/evaluation geometry shift**
   - matched geometry distribution;
   - evaluation systems shifted toward lower density / stronger clustering than training systems.

## Graph scaling

The graph rule remains density-scaled rather than fixed-k:

`k_s = max(2, round(0.15 * n_s))`

subject to `k_s < n_s`.

No k, bandwidth, or geometry-specific tuning is permitted after Q5 outcomes are opened.

## Mandatory cells

Q5-NULL-COMPONENT: shared component turnover but no mismatch/relation front, under heterogeneous geometry.

Q5-NULL-PRIVATE-M: strong system-private mismatch fronts with heterogeneous geometry.

Q5-NULL-PRIVATE-C: strong system-private relation rotations with heterogeneous geometry.

Q5-POWER-M: shared mismatch-magnitude front under heterogeneous geometry.

Q5-POWER-C: shared relation front at constant mismatch magnitude under heterogeneous geometry.

Q5-SHIFT-M: Q5-POWER-M with train/evaluation geometry shift.

Q5-SHIFT-C: Q5-POWER-C with train/evaluation geometry shift.

## Formal qualification thresholds

When Q5 is authorized, use the same operating criteria as the idealized gate unless a different rule is frozen before any Q5 result exists:

- alpha = 0.05;
- each type-I cell: Wilson 95% upper rejection bound `<= 0.10`;
- each power cell: Wilson 95% lower rejection bound `>= 0.80`.

World count, bootstrap count, exact geometry mixture weights, random seeds, bandwidth, and batching must be written to a machine-readable rule and exact blobs authorized before execution.

## Interpretation

- Idealized PASS + Q5 FAIL means the estimand/inference works in balanced geometry but is not yet qualified for heterogeneous opportunistic data.
- Q5 PASS supports robustness to the declared synthetic geometry family only; it does not establish validity on any particular empirical dataset.
- Real-data use still requires an observation-support/classifiability gate defined independently of the empirical outcome.

## Claim ceiling

Q5 is not a universal robustness theorem. It qualifies a finite predeclared family of sampling geometries. It must not be used to claim robustness to untested spatial bias, taxonomic bias, temporal bias, measurement error, or arbitrary missing-not-at-random processes.