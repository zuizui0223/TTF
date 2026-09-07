# Mesoscopic TTF candidate: geometry-conditioned consensus transfer

Status: **candidate / pre-qualification**. This note was committed while the first v0.3a Queensland mesoscopic workflow (`34124819256`) was still queued. It therefore records a prospective design issue and candidate remedy before any v0.3a calibration result was inspected.

Empirical genetic outcomes and named Queensland barriers remain sealed.

## Why a mesoscopic regime is needed

The Queensland frozen topology has eight ordered basin columns and therefore only seven candidate inter-basin cuts. The v0.2 continuous kernel field had acceptable private-world type-I behavior but insufficient full-shared power. Noise-free and oracle diagnostics then separated two problems:

1. learning a continuous field from sparse discrete topology loses transferable information;
2. a zero-centered held-out-score null does not explicitly represent finite-candidate chance commonness.

The first mesoscopic candidate (v0.3a) moved the field onto the seven shared cut elements and compared held-out evidence against a uniform private-cut mixture. Before its Gate ran, a further issue was identified: v0.3a uses a **per-edge average** profile-likelihood gain and then exponentiates it. That can over-temper within-species certainty. Conversely, simply multiplying the gain by the number of edges and pooling those values directly would let densely sampled species dominate the training field.

The candidate below separates within-species certainty from between-species weighting.

## 1. Species-specific cut posterior

For species s and candidate cut c that its fixed graph can contrast, let

\[
\ell_{sc}=\frac{|E_s|}{2}\log\frac{\mathrm{SSE}_{0s}+\epsilon}{\mathrm{SSE}_{1sc}+\epsilon}.
\]

This is the total Gaussian profile log-likelihood improvement corresponding to the already implemented per-edge gain. Let `O_s` denote the cuts observable by species s. The geometry-only private prior is

\[
u_{sc}=\frac{1}{|O_s|},\qquad c\in O_s.
\]

Define a species-specific soft boundary label

\[
q_{sc}
=\frac{u_{sc}\exp(\ell_{sc})}
{\sum_{j\in O_s}u_{sj}\exp(\ell_{sj})}.
\]

Thus more informative within-species data can make **that species' posterior** sharper, but every species still contributes one normalized distribution downstream.

## 2. Species-equal consensus field with observability conditioning

Let the global boundary field be a probability vector p over the common mesh. A species that cannot observe a cut must not contribute negative evidence against it. For each training species, compare its soft label q_s only with the field restricted to O_s:

\[
r_{sc}(p)=\frac{p_c}{\sum_{j\in O_s}p_j},\qquad c\in O_s.
\]

Estimate p by maximizing the species-equal conditional log score

\[
\sum_{s\in S_{train}}
\sum_{c\in O_s}q_{sc}\log r_{sc}(p)
+\lambda\sum_c\frac{1}{C}\log p_c,
\]

where C is the number of mesh elements. The second term is equivalent to a uniform pseudo-species prior and prevents degenerate infinite logits.

This objective has a useful null property. If every training species has exactly its geometry-private distribution, `q_s=u_s`, then the uniform global field is a stationary optimum **for every observability pattern separately**. Missing cuts therefore do not create a preferred boundary merely because they are observed by fewer species.

## 3. Held-out transfer score

For a held-out species, restrict the trained field to its own O_s and define

\[
C_s=
\sum_{c\in O_s} q_{sc}
\log\frac{r_{sc}(\hat p)}{u_{sc}}.
\]

The primary statistic remains the species macro-average

\[
T=|S_{eval}|^{-1}\sum_s C_s.
\]

If the private-boundary prior is calibrated so that

\[
E(q_s\mid O_s,H_0)=u_s,
\]

then conditional on any frozen training field,

\[
E(C_s\mid \hat p,O_s,H_0)
=-\mathrm{KL}(u_s\|r_s(\hat p))\le 0.
\]

So finite-candidate chance matching is not treated as evidence for sharedness. This is the key theoretical advantage over a generic zero-centered correlation statistic.

## 4. What still requires calibration

The inequality above depends on the species posterior being calibrated under the geometry-private null. Edge ranks are dependent and the Gaussian profile likelihood is a working likelihood, so empirical type-I control must still be demonstrated rather than assumed.

Required synthetic gates before empirical outcomes:

- exact frozen 12-layout Queensland envelope;
- same 10/11 species split and k=2 graphs;
- private worlds across multiple positive amplitudes;
- full-shared A=2 positive control;
- partial-shared response as a diagnostic;
- independent calibration seed from v0.3a;
- no tuning selected using genetic outcomes or named barrier labels.

If this candidate is pursued after v0.3a, it must be treated as a distinct estimator version and qualified on a fresh synthetic seed. The v0.3a Gate result must be retained regardless of outcome.
