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

## v0.7 result: observable nuisance strength fixes calibration, not missing support

v0.7 measured private spatial strength from the training species only.  The
summary was the equal-species mean local coherence of edge-length-orthogonalized
turnover over nearest edge midpoints.  Held-out transfer was not used to select
the nuisance configuration.  Every private reference was split prospectively:
the first 999 worlds estimated the training-strength profile and the remaining
1,000 worlds supplied independent raw-`T` calibration values.

This profiling step did what it was designed to do.  On v03, the maximum
private-world rejection rate was 0.048 with Wilson upper 0.0704; on v04 it was
0.056 with Wilson upper 0.0797.  Matching-amplitude recovery was high for the
stronger private controls (v03: 0.882 at A2 and 0.964 at A3; v04: 0.786 at A2
and 0.960 at A3).

But power did not clear the old universal gate.  Shared-A2 rejection was 0.790
on v03 (Wilson lower 0.7522) and 0.492 on v04 (Wilson lower 0.4484).  Crucially,
the v03 value is essentially the same as the non-selectable matched-amplitude
oracle from v0.6 (0.798), while the v04 oracle is only 0.540.  Therefore the
remaining power deficit cannot reasonably be attributed to nuisance profiling.
It is a finite-design detectability limit: the common transition is insufficiently
represented by some sampling geometries.

## Consequence: validity and detectability are different gates

It is not scientifically coherent to require one method to have the same power
on every fixed sampling geometry, including geometries that barely intersect the
shared transition.  It is equally incoherent to relax type-I control whenever a
geometry is difficult.

The v0.8 architecture therefore separates two claims:

1. **Inferential validity.**  Strong species-private transitions must not create
   false sharedness.  This is a method property and is confirmed on fresh
   external geometries without dropping panels for low power.
2. **Deployment detectability.**  Before empirical trait outcomes are opened on
   a concrete deployment geometry, that exact fixed design must pass a
   prospectively specified synthetic power audit.  A design that fails is called
   non-qualifying / underpowered; its absence of rejection cannot be interpreted
   as evidence for no sharedness.

The numerical thresholds are unchanged: maximum zero-shared Wilson 95% upper
bound must be <=0.10, and a deployment geometry must additionally have
full-shared A2 Wilson 95% lower bound >=0.80.

This is not a rescue of v0.7.  v0.7 remains a failed candidate under its frozen
old rule.  v0.8 is a new qualification architecture, frozen before its fresh
external validity panels are declared.

## Current frozen boundaries

- v0.2, v0.3 and v0.4 confirmatory failures remain failures.
- v0.5 selected no correction candidate.
- v0.6 selected no unrestricted private-envelope candidate.
- v0.7 selected no candidate under the old conjunctive type-I-plus-universal-power rule.
- The old birds and butterflies holdouts remain performance-unopened and are not
  eligible for v0.8 external validity confirmation because their earlier role
  was frozen under a different conjunctive policy and their geometry-only
  observability has already been inspected.
- The RGFCA reserve remains performance-unopened and may only become a v0.8
  deployment-adequacy geometry after fresh v0.8 external validity passes.

The identifiability limits above narrow the scientific claim; they do not rescue
any failed gate retroactively.
