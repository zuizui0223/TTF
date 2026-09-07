# Mesoscopic TTF v0.3 — discrete boundary-element regime

Status: **method-development candidate; empirical genetic outcomes remain unopened**.

## Why a second spatial resolution regime is needed

Continuous TTF v0.2 learns a kernel boundary field and is qualified on dense idealized geography. The Queensland stress panel exposes a different regime: each species occupies only a subset of eight ordered drainage basins, so the common spatial support contains only seven candidate inter-basin cuts.

In such coarse geography, two independent species can choose the same transition location by chance with non-negligible probability. A null centered at zero correlation therefore does not automatically represent "no cross-species shared boundary".

Mesoscopic TTF treats the candidate boundary elements themselves as the common prediction space.

## Candidate boundary mesh

Let the common mesh contain boundary elements

\[
\mathcal C=\{c_1,\ldots,c_C\}.
\]

For Queensland, `C=7`, corresponding only to the seven anonymous adjacent positions between eight ordered basin columns. Named barriers and empirical genetic outcomes are not used to construct or select the mesh.

Each species has an observability set

\[
O_s\subseteq\mathcal C,
\]

containing only cuts for which its fixed within-species graph supplies both crossing and non-crossing edges. A cut outside `O_s` is **unobserved**, not evidence against a boundary.

## Species-specific cut evidence

For each observable cut, compare a two-level edge-turnover model (crossing vs non-crossing) with an intercept-only model using within-species rank turnover. The current implementation records the Gaussian profile log-likelihood improvement per edge,

\[
g_s(c)=\frac12\log\frac{SSE_{0,s}+\epsilon}{SSE_{1,s}(c)+\epsilon},
\]

with zero evidence assigned to unobservable cuts.

The v0.3a field averages these per-species log-evidence values opportunity-conditionally and applies a softmax over cuts.

## v0.3a predictive score

For a held-out species, restrict the learned training distribution `p_train(c)` to its observable set `O_s`. Compare its evidence under the learned distribution with a private-boundary baseline that is uniform on the same observable set:

\[
C_s=
\log\sum_{c\in O_s}p_{train,s}(c)e^{g_s(c)}
-
\log\left[\frac{1}{|O_s|}\sum_{c\in O_s}e^{g_s(c)}\right].
\]

This directly includes finite-cut chance coincidence in the baseline.

Inference still treats held-out species, not edges or records, as the sampling units.

## Prospectively registered v0.3b consensus candidate

The v0.3b candidate was implemented and tested **before v0.3a qualification results were available**.

First restore total profile evidence for species `s` by multiplying the per-edge evidence by its fixed graph edge count. Normalize within species to obtain a cut posterior-like distribution

\[
q_s(c)\propto u_s(c)\exp\{|E_s|g_s(c)\},\qquad c\in O_s,
\]

where `u_s` is uniform on `O_s`. Each species then contributes exactly one normalized distribution, so extra edges sharpen that species' location evidence but do not increase its downstream species weight.

The global consensus distribution `p(c)` maximizes the sum of observability-conditioned cross-entropies,

\[
\sum_s\sum_{c\in O_s}q_s(c)\log p_s(c)
+\lambda\sum_c\frac1C\log p(c),
\]

where `p_s` is `p` renormalized on `O_s` and the second term is a uniform pseudo-species prior.

A key invariant is that if every species is geometry-private/uninformative, `q_s=u_s`, then the fitted global field is exactly uniform even when species have different missing-cut patterns.

## v0.3b held-out score and null geometry

For held-out posterior `q_s`, define

\[
C_s=\sum_{c\in O_s}q_s(c)
\left[\log p_s(c)-\log u_s(c)\right].
\]

If a held-out species follows the geometry-private reference exactly (`q_s=u_s`), then

\[
C_s=-D_{KL}(u_s\Vert p_s)\le0.
\]

Thus the no-shared-boundary geometry is built into the score itself rather than assumed to have zero correlation. The primary test remains a one-sided held-out-species bootstrap of positive mean score.

## Frozen Queensland qualification contract

The empirical outcome remains sealed during method qualification.

- frozen topology envelope: 12 feasible ordered layouts;
- eligible species: 21, chosen only from sampling support;
- split: 10 training / 11 held-out species;
- local graph: `k=2`;
- candidate cuts: 7 anonymous inter-basin positions;
- synthetic noise SD: 0.8;
- transition width: 0.20;
- private controls: shared fraction 0 at amplitudes 0.5, 1, 2, 3, 4;
- response curve: shared fractions 0.25, 0.5, 1.0 at amplitude 2;
- qualification Gate: every layout must have maximum private rejection `<=0.10` and full-shared A=2 power `>=0.80`;
- no named Mary–Brisbane boundary and no empirical genetic distance is used in fitting, selection, or qualification.

v0.3a uses master seed `20260910`. The prospectively registered v0.3b candidate uses independent master seed `20260911` and `evidence_temperature=1.0`; that temperature is not tuned against Queensland outcomes.

## Interpretation rule

A method version that fails the frozen Gate remains a documented negative method-development result. It is not rescued by selecting favorable layouts, known barriers, species with favorable genetic results, or post-hoc evidence temperatures.

Only a version that passes on the frozen synthetic geometry may proceed to independent-seed confirmation and then to empirical genetic outcomes.
