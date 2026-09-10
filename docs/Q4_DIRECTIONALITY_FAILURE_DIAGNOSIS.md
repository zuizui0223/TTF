# Q4 directionality v0.1 — immutable failure diagnosis

Status: **post-outcome diagnosis of an immutable failed gate**.

This document does not change, rerun, reinterpret as PASS, or delete `results/ttfm_q4_directionality_v0.1.json`. Q4 v0.1 remains formally failed. Q5 and empirical applications remain unauthorized.

## 1. Frozen outcome

The predeclared Q4 v0.1 gate failed because `Q4_PRIVATE_private_directional_change` exceeded both relevant type-I ceilings:

- direction-free TTF-C relation detection: 73/500 = 0.146; Wilson 95% upper = 0.1796492662;
- breakdown label: 58/500 = 0.116; Wilson 95% upper = 0.1470418370;
- recoupling label: 0/500.

The other mandatory controls behaved as intended: shared breakdown and recoupling were directionally separated, neutral constant-mismatch relation rotation produced no directional label, and a 50:50 sign-conflict front retained direction-free relation transfer but produced neither shared directional label.

The legal gate result is therefore `Q4_pass = false`.

## 2. Why the failed PRIVATE cell is not a clean TTF null

In the frozen Q4 simulator, every PRIVATE system receives an independently sampled front center

```text
Z_s ~ Uniform(-0.75, 0.75)
```

on an observation domain approximately `[-1, 1]`. Exact front locations are independent across systems, but their population distribution is not translation-invariant over the observed domain: transition fronts are concentrated in the central band and excluded near the outer domain.

TTF does not test the proposition

> all systems contain exactly the same deterministic front.

It tests whether a turnover-intensity field learned from training systems predicts where turnover is high within unseen systems. If a system-specific front has local turnover profile `K_w(x-Z_s)` and front location density `p(z)`, then schematically

\[
E\{u_s(x)\} \propto (K_w * p)(x).
\]

A non-uniform `p(z)` therefore induces a non-constant population turnover-intensity field even when `Z_s` values are independent. That field can legitimately transfer to held-out systems.

Thus

```text
independent exact front locations
    !=
absence of transferable turnover-intensity structure.
```

The Q4 PRIVATE construction mixed these two concepts. Its elevated rejection rate cannot, by itself, be assigned wholly to anti-conservative inference.

## 3. Revised semantic hierarchy

The method must distinguish at least three nested notions:

1. **exact-front sharing** — the same deterministic front location occurs across systems;
2. **front-density / turnover-intensity sharing** — exact fronts vary, but their spatial distribution is non-uniform in a way that predicts unseen-system turnover;
3. **strict private structure** — system-specific transitions contain no non-constant population turnover-intensity field on the comparison domain.

Core TTF, TTF-M, and TTF-C target the second/general notion: held-out transferability of turnover intensity. Exact-front sharing is a special case, not the estimand definition.

## 4. Consequence for a valid strict-private null

A replacement null must make the expected turnover-intensity field constant under the geometry used for transfer. One clean construction is a translation/rotation-invariant domain:

- generate a system-specific front anchor independently and uniformly over a full periodic spatial domain;
- express absolute transfer coordinates on that periodic domain;
- define a separate predeclared local orientation coordinate relative to the anchor, generated before state outcomes;
- generate the same low-to-high mismatch direction in local orientation coordinates.

The local directional sign can therefore be common while the absolute spatial front is genuinely private. Because the anchor phase is uniform over the full periodic domain, the population turnover-intensity field is rotationally invariant rather than concentrated in a shared zone.

This is conceptually different from the failed Q4 PRIVATE world and must be treated as a **new prospective gate**, not a rerun or rescue of Q4 v0.1.

## 5. A useful alternative, not a null

The original bounded-center construction should be retained as a distinct `shared_front_zone` alternative. It represents systems whose exact transition locations differ but whose transitions preferentially occur in a common region. Detecting such a world is compatible with the TTF estimand and helps define its spatial-scale semantics.

It must not be used as a type-I null after this diagnosis.

## 6. No tuning allowed

The following are forbidden as responses to the Q4 v0.1 failure:

- changing alpha, Wilson ceilings/floors, sign-consistency fraction, or bootstrap counts and rerunning the same gate;
- shrinking or widening the PRIVATE front range to obtain a preferred rejection rate;
- changing the front window or locator after inspecting Q4 performance and calling the result Q4 v0.1;
- relabelling the frozen failed result as passed.

A new protocol must use fresh independent seeds, an explicitly strict-private construction, its own frozen rule and exact-blob authorization.

## 7. Current claim ceiling

The idealized direction-free TTF-M/TTF-C gate remains passed. Q4 v0.1 remains failed. Therefore the method is currently qualified only for direction-free mismatch/coupling transition transfer on the declared idealized family. It is **not yet qualified for transferable breakdown/recoupling labels**, heterogeneous geometry (Q5), or any empirical mismatch application.