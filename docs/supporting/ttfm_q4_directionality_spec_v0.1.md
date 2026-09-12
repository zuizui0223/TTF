# TTF-C Q4 directionality gate v0.1

Status: **prospective design only; not yet executed**.

This document extends the already-qualified direction-free TTF-C relation-front estimand with an explicitly oriented deterioration/recoupling label. It does not alter TTF-M, TTF-C, or any frozen v0.1-v0.12 result.

## Why Q4 exists

TTF-C answers whether the relation state between paired fields changes at a transferable front. A transferable relation change is not, by itself, a breakdown. The same front can represent deterioration, recoupling, or a direction-neutral rotation.

Therefore a breakdown claim requires a separately frozen orientation rule.

## Required inputs

For system `s` and location `x`:

- paired states `A_s(x)` and `B_s(x)`;
- predeclared mismatch magnitude `M_s(x)=m(A_s(x),B_s(x))`;
- predeclared relation state `R_s(x)=r(A_s(x),B_s(x))`;
- predeclared signed orientation coordinate `g_s(x)`.

`g_s(x)` must be defined independently of the observed mismatch outcome. Examples include increasing isolation, increasing warming exposure, distance from a known boundary, or a prospectively specified intervention/exposure direction.

If no scientifically defensible signed orientation exists, Q4 is not evaluable and TTF-C remains direction-free.

## Direction score

After a TTF-C transferable relation front has been identified under the frozen TTF-C procedure, define a local oriented mismatch contrast around that front:

`Delta_M = mean(M on the positive side of g) - mean(M on the negative side of g)`

using a predeclared symmetric front neighbourhood and the same held-out systems used for transfer evaluation.

Interpretation:

- `Delta_M > 0`: deterioration / breakdown direction;
- `Delta_M < 0`: recoupling / recovery direction;
- `Delta_M ~= 0`: relation change without directional mismatch deterioration.

The sign threshold and neighbourhood width must be frozen before any empirical use.

## Mandatory synthetic families

Q4-BREAK: shared relation front plus increasing mismatch on the positive side of `g`.

Q4-RECOUPLE: the same front geometry with decreasing mismatch on the positive side of `g`.

Q4-ROTATE: shared relation rotation with constant mismatch magnitude. This must not be labelled breakdown or recoupling.

Q4-PRIVATE: system-private oriented mismatch changes with no shared front. This must not produce a transferable directional front.

Q4-FLIP: half the systems deteriorate and half recouple across the same geometric front. This may show a direction-free TTF-C front but must fail a shared directional-breakdown claim.

## Qualification logic

Q4 is a two-part gate:

1. the underlying TTF-C front must already satisfy the frozen transferable-front qualification;
2. the oriented mismatch sign must generalize to held-out systems under the predeclared direction rule.

A front is labelled transferable breakdown only when both conditions pass.

## Claim ceiling

Q4 does not establish causality. A signed exposure axis is an orientation device, not an intervention. `breakdown` here means transferable increase in predeclared mismatch across an oriented front, not demonstrated fitness loss or mechanism failure.

No empirical flower-colour, pollinator, island, or mismatch dataset is authorized by this specification.