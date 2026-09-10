# TTF-C Q4.1 strict-private directionality qualification — prospective spec v0.1

Status: **prospective design; no Q4.1 outcomes opened**.

Q4 v0.1 remains immutable FAIL. Q4.1 is a new protocol motivated by the estimand clarification in `docs/Q4_DIRECTIONALITY_FAILURE_DIAGNOSIS.md`; it is not a rerun, threshold adjustment, or rescue of Q4 v0.1.

## 1. Target

Qualify a directional `breakdown` / `recoupling` label only when:

1. direction-free TTF-C relation turnover transfers to held-out systems; and
2. mismatch change across a separately predeclared orientation coordinate has a consistent held-out sign.

TTF-C itself remains a transferable **turnover-intensity field** estimator, not an exact-boundary identity test.

## 2. Strict-private null requirement

A private-front type-I control must have strong relation/mismatch transitions within every system but no non-constant population turnover-intensity field in the absolute spatial coordinates used by TTF-C.

Q4.1 therefore uses a rotation-invariant construction:

- each system receives an anchor phase `phi_s ~ Uniform(0, 2*pi)` before outcomes;
- records receive a local orientation coordinate `g_s in [-1,1]`;
- absolute transfer coordinates are an arc on the unit circle centered at `phi_s`;
- the system's low-to-high mismatch transition occurs at `g_s=0`;
- `phi_s` is independent across systems and uniform over the full circle.

The local direction is therefore common, but the absolute front phase is strict-private. Rotational symmetry makes the population expected absolute turnover-intensity field constant in phase.

This is not the same world as Q4 v0.1 PRIVATE, whose independent fronts were restricted to a common central subdomain.

## 3. Shared-front-zone semantic control

The old bounded-center construction is retained under a new explicit interpretation:

`shared_front_zone_directional_change`

Exact front locations differ across systems, but all are drawn from a common restricted spatial zone. This is **not a type-I null**. It is a characterization world showing that TTF may detect transferable front-density / turnover-intensity structure without exact front identity.

Its outcome must not enter Q4.1 PASS/FAIL.

## 4. Mandatory qualification cells

The formal Q4.1 gate will prospectively freeze fresh independent seeds and include:

- **BREAK:** shared spatial relation front; mismatch increases along positive `g`; relation and breakdown power required, recoupling type-I controlled.
- **RECOUPLE:** same shared spatial front; mismatch decreases along positive `g`; relation and recoupling power required, breakdown type-I controlled.
- **ROTATE:** shared relation rotation at constant mismatch magnitude; relation power required, both directional labels type-I controlled.
- **STRICT_PRIVATE:** rotation-invariant absolute front phase with common local low-to-high orientation; relation, breakdown, and recoupling type-I controlled.
- **FLIP:** shared spatial relation front but 50:50 breakdown/recoupling signs; relation power required, both directional labels type-I controlled.
- **ZONE:** bounded independent exact fronts in a shared spatial zone; characterization only, excluded from gate logic.

## 5. Frozen-before-outcome quantities required

Before Q4.1 execution, a machine-readable rule must freeze:

- systems, records, train/evaluation split;
- spatial arc width and full-circle anchor distribution for STRICT_PRIVATE;
- graph density scaling, kernel bandwidth and priors;
- local front neighbourhood width;
- minimum records per side;
- sign-consistency requirement;
- replicate counts, bootstrap counts and seeds;
- alpha and Wilson type-I/power thresholds;
- exact role of every cell, including `ZONE = characterization_only`;
- exact code/workflow blob SHAs in a one-run authorization.

No value may be chosen by inspecting Q4.1 outcomes.

## 6. Pass logic

Every declared qualification role must pass. The intended convention remains the existing TTF qualification convention unless a different value is frozen before outcomes:

- type-I: Wilson 95% upper <= 0.10;
- power: Wilson 95% lower >= 0.80;
- alpha = 0.05.

`ZONE` cannot rescue or fail the gate.

## 7. Claim ceiling

A Q4.1 PASS would establish directional operating behavior only on the declared synthetic family. It would not establish causality, fitness loss, mechanism failure, heterogeneous/opportunistic geometry robustness, or empirical ecological validity. Q5 remains a separate required gate. No flower-colour, pollinator, island, or other empirical mismatch dataset is authorized here.