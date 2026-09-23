# Structural / EGWE state-to-TTF handoff

Status: external predictor and transition handoff. This note specifies how a prospectively frozen Structural/EGWE state or typed-connectivity representation may enter TTF without changing TTF's core out-of-species transfer estimand or merging empirical denominators.

## 1. Roles

Structural asks whether a spatial or connectivity coordinate adds endpoint-relevant information within a system after the declared reference is supplied.

EGWE asks whether the current state is sufficient for future fate, or whether origin/history retains residual information.

TTF asks whether transition or turnover structure learned from some species predicts structure in entirely unseen species.

The core distinction is:

adequacy != portability != transferability

## 2. Upstream objects

For species s, an upstream protocol may supply:

- R_s(x): declared current state/reference;
- C_s^p(x): typed connectivity for a declared biological operator p;
- Z_s(x): prospectively defined future-transition state, when a transition-field lane is used.

These objects must be frozen before held-out evaluation-species outcomes are used to fit any shared field.

## 3. Lane A: transition-field transfer

Define a prospectively fixed within-species edge target U_se from the future-transition state at the two edge endpoints.

TTF then learns transition turnover from training species and evaluates it on species completely excluded from training.

Question: does spatial organization of future transitions transfer across unseen species?

## 4. Lane B: typed-connectivity predictor competition

Convert R and C into response-independent within-species edge features.

Use the existing predictor-space competition with two spaces:

- R: current-state/reference features;
- RC: the same R features plus typed connectivity C.

The primary transfer increment is C|R = T_RC - T_R.

Question: does typed connectivity add out-of-species transfer skill beyond the declared present-state/reference space?

This is not the same estimand as a within-species Structural C-R loss difference.

## 5. Species eligibility anti-selection

Invalid workflow:

run Structural in many species; keep only species where C-R is favorable; then run TTF.

That conditions the TTF population on the response.

Species eligibility may depend only on prospectively frozen support criteria such as provenance, geometry/admissibility, measurement support and endpoint-transition estimability.

Species with null or adverse within-species Structural results remain eligible if they pass those support rules.

## 6. Feature anti-leakage

A typed connectivity feature is valid for TTF evaluation only when it is fixed from external theory/biology, fit on TTF training species only, or cross-fitted with each evaluation species excluded.

Evaluation-species outcomes may not be used to choose radii, operator weights, graph scales, transforms, thresholds or reference variables.

## 7. Burned pilots

Burned transition pilots are feasibility-only.

They do not enter the TTF statistic, predictor-space increments or predictive denominator, and they are never pooled with confirmatory data.

## 8. Self-detectability

A qualified non-transfer result is interpretable as system- or lineage-conditioned only when within-species self structure is separately qualified and supported.

The current genetic flagship already demonstrates this distinction in a separate domain: fresh cross-species post-IBD transfer was not detected, while the separately qualified within-species self test was significant relative to its frozen structural null.

That empirical result motivates the logic here. It is not a Structural replication.

## 9. Independent TTF qualification

Passing Structural/EGWE gates does not automatically qualify TTF. Cross-species geometry, observational support, Type-I control and power remain TTF-specific gates.

A full handoff therefore requires both upstream state/transition qualification and TTF deployment qualification.

## 10. Interpretation

Within-system structure supported plus cross-species transfer supported: candidate reusable transition rule in the tested domain.

Within-system structure supported or self-detectable plus transfer not detected: system- or lineage-conditioned transition structure in the tested domain.

Upstream endpoint not evaluable: abstain; this is not a negative TTF result.

TTF transfer with upstream adequacy not earned: the TTF predictor result may stand at its own predictive level but cannot rescue the failed within-system adequacy claim.

## 11. Current empirical boundary

The genetic TTF result is its own empirical denominator.

It is not a replication of Structural A-Islands or Tanzania, and PNW/RMNP do not enter the TTF transfer denominator.

## 12. Claim ceiling

A successful future handoff can support the statement that a prospectively specified current-state/connectivity representation adds out-of-species predictive skill for within-species future-transition structure beyond a declared current-state reference.

It does not by itself prove universal causal connectivity, identical operators, identical histories, or transfer outside the frozen taxa/scale/endpoint domain.
