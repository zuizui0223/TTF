# Identifiability limit for transferable turnover under unrestricted private strength

## Why Gate I did not reduce to a calibration bug

TTF asks whether a turnover field learned from one set of species predicts turnover
in unseen species.  Let `T` denote the upper-tail transfer statistic used for that
question.  Versions v0.2--v0.5 treated zero as the relevant null center, or tried
to remove geometry-derived positive structure from `T`.  Gate-I semi-synthetic
worlds showed that strong species-private transitions can instead produce a
positive transfer baseline on real sampling geometry.

v0.6 kept the raw transfer statistic and enlarged the null to an explicitly
species-private transition family.  This resolved type-I error when the private
signal strength was known (the non-selectable oracle test), but a worst-case null
over arbitrarily strong private transitions had zero useful power: the
noise-free private configuration was least favourable for essentially every
moderate shared world.

That failure is an identifiability result, not a reason to tune the threshold.

## Lemma 1: worst-case critical value for an upper-tail composite null

Let `P_eta` be the distribution of a scalar statistic `T` under a family of null
models indexed by nuisance parameter `eta`, and let `Q` be an alternative.  Any
upper-tail test of the form

`phi_c(T) = 1{T >= c}`

that has size at most `alpha` uniformly over the composite null must satisfy

`c >= sup_eta q_(1-alpha)(P_eta)`.

Therefore its power under `Q` is bounded by

`Q[T >= sup_eta q_(1-alpha)(P_eta)]`.

### Proof

For every `eta`, size control requires
`P_eta[T >= c] <= alpha`, hence `c >= q_(1-alpha)(P_eta)` (up to the usual
finite-sample/tie convention).  Taking the supremum over `eta` gives the first
claim.  Since an upper-tail rejection region shrinks monotonically as `c`
increases, substituting that lower bound on `c` gives the power bound.  QED.

## Corollary 1: unrestricted private strength can make moderate sharedness untestable by `T` alone

If the private-null family contains a configuration whose upper tail lies to the
right of the moderate shared alternative, then no uniformly valid monotone test
based only on `T` can simultaneously retain that private configuration in the
null and have high power for the moderate shared alternative.

This is exactly what the frozen v0.6 development result demonstrates.  The
`infinite_snr` private reference was the least-favourable configuration for all
500 v03 shared-A2 worlds and 497/500 v04 shared-A2 worlds.  The resulting
worst-case envelope rejected 0/500 shared-A2 worlds on both panels.  The same
raw statistic under the non-selectable matched-amplitude oracle controlled
private type-I error, showing that the failure is nuisance non-identifiability,
not a broken Monte Carlo implementation.

## Lemma 2: geometry with no shared-transition support cannot be repaired by inference

Fix a shared boundary and a held-out sampling design.  If the observed graph for
a held-out species contains no trait-informative contrast induced by that
boundary, then the distribution of any statistic measurable from the observed
traits and that fixed graph is unchanged by adding the unobserved portion of the
shared boundary.  No inferential procedure can recover power from support that
was never sampled.

The four-panel outcome-free observability audit is a finite-design diagnostic of
this condition.  It is not a selection rule.  In the frozen audit, mean
noise-free shared-edge contrast was about 0.1035 on v03 but 0.06995 on v04,
matching the direction of the much lower raw shared-control power on v04.

## Consequence for the method claim

The defensible target is not

> sharedness is identifiable against arbitrarily strong, completely unrestricted
> private spatial structure on every sampling geometry.

That target is false for the current statistic and, without additional nuisance
information or support assumptions, can be impossible in principle.

The next admissible target is conditional:

> given a frozen training field, sufficient held-out sampling support, and an
> observable estimate of within-species private spatial strength that does not use
> held-out cross-species alignment, does unseen-species transfer exceed the
> species-private baseline compatible with that strength?

The conditioning variable must be measurable without using the held-out
cross-species transfer outcome.  It cannot be the simulator amplitude itself.
Any profiling or conditioning rule must be fixed and validated on failed
development panels before the predeclared birds, butterflies, or RGFCA-reserve
performance is opened.

## Current frozen boundaries

- v0.2, v0.3 and v0.4 confirmatory failures remain failures.
- v0.5 selected no correction candidate.
- v0.6 selected no unrestricted private-envelope candidate.
- Birds and butterflies remain performance-unopened confirmation geometries.
- The RGFCA reserve remains performance-unopened.

The identifiability limit above narrows the scientific claim; it does not rescue
any failed gate retroactively.
