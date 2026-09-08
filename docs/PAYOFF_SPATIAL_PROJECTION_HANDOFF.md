# PAYOFF spatial-projection handoff

Status: external mechanistic handoff. This note specifies how a PAYOFF-derived local architecture margin can enter TTF as a prospectively fixed predictor space without changing TTF's core sharedness estimand or claim ceiling.

The motivating individual-system anchor is the Izu-island *Campanula microdonta* floral-syndrome system in `zuizui0223/payoff/docs/CAMPANULA_SYNDROME_DISASSEMBLY_HANDOFF.md`.

The full bridge is specified in `zuizui0223/payoff/docs/PAYOFF_TTF_SPATIAL_PROJECTION_BRIDGE.md`.

## 1. Separation of roles

```text
Campanula / individual-system anchor
    deep biological resolution within one species
            |
            v
PAYOFF
    local architecture payoff margin
            |
            v
spatial projection
    location -> ecological context -> payoff field
            |
            v
TTF
    cross-species transfer on held-out species
```

TTF does not estimate PAYOFF parameters by default. It evaluates whether a predictor or learned transition field transfers to unseen species.

---

## 2. Mechanism-side object

For species `s`, define an externally justified local margin

```text
M_s(e,p)=phi_s(e)+eta_s(e)(2p-1).
```

Project ecological context into geographic space:

```text
m_s(x)=M_s(e_s(x),p_s(x)).
```

Before frequency dependence is identified, use the explicitly weaker static projection

```text
m_s^0(x)=phi_s(e_s(x)).
```

The latent mechanism boundary is

```text
B_s^P={x : m_s(x)=0}.
```

This is not the TTF boundary field `B_train(x)`. TTF must keep the learned turnover field and the mechanism-derived field conceptually separate.

---

## 3. Edge feature contract

TTF operates on within-species edges. Convert the spatial payoff field into one or more edge features without using held-out trait outcomes.

A generic scalar mechanism exposure is

```text
P_se
= (1/L_se) integral_e g[m_s(x)] dl,
```

where `g` is frozen prospectively.

Possible transforms include

```text
g(m)=exp(-|m|/tau)
```

or a prespecified near-boundary indicator.

No transform is canonical. The feature definition must be calibrated for the intended sampling geometry and fixed before evaluation-trait inspection.

If several PAYOFF-derived quantities are retained, they may occupy several columns in the generic TTF feature matrix.

---

## 4. Predictor-space competition

The current `compete_predictor_spaces` implementation accepts arbitrary named column sets, so no core code change is required to add a PAYOFF-derived space.

Recommended comparison:

```python
spaces = {
    "G":   [...],
    "E":   [...],
    "GE":  [...],
    "P":   [...],
    "GEP": [...],
}

increments = {
    "P|GE":  ("GEP", "GE"),
    "GE|P":  ("GEP", "P"),
}
```

The main mechanistic increment is

```text
Delta_P|GE = T_GEP - T_GE.
```

All spaces must use the same species-disjoint train/evaluation split.

A positive increment means that the preregistered PAYOFF projection adds held-out edge-level turnover prediction beyond the declared geography/environment baseline.

---

## 5. Anti-leakage rule

The PAYOFF feature is valid for held-out evaluation only if the evaluation species' trait outcomes did not determine the feature.

Allowed constructions include:

```text
1. fixed from prior biological theory and external parameter estimates;
2. fitted only on training species;
3. cross-fitted with evaluation species excluded from every fitting step.
```

Disallowed construction:

```text
use evaluation-species turnover to choose driver weights,
thresholds, transforms or PAYOFF parameters,
then score those same species as held out.
```

If species-specific nuisance models are needed, they must follow the same cross-fitting discipline already required by TTF for intrinsic-IBD residualization.

---

## 6. Campanula is one anchor species, not eight TTF species

The Izu *C. microdonta* populations provide a deep mechanistic anchor across islands and floral modules.

They can establish:

```text
module identity;
trait direction;
allometry-adjusted disassembly;
population-history alternatives;
repeated origin versus shared ancestry;
and eventually direct reproductive payoff.
```

But the eight island populations remain populations of one species. They must not be counted as eight independent held-out species for TTF's superpopulation sharedness inference.

A valid cross-species TTF extension requires additional species with their own within-species georeferenced trait-transition graphs.

---

## 7. Two TTF lanes must remain distinct

### Lane U — unsigned transferable turnover

Canonical TTF:

```text
train-species turnover -> B_train(x)
                      -> held-out species edge scores
                      -> T.
```

Question:

> Does within-species turnover geometry transfer across species?

### Lane P — PAYOFF-informed mechanistic projection

External mechanism feature:

```text
PAYOFF receipts / ecological drivers
                -> m_s(x)
                -> P_se
                -> held-out predictor competition.
```

Question:

> Does a mechanism-derived payoff projection add transfer skill beyond baseline predictor spaces?

A result can be positive in one lane and negative in the other. That distinction is scientifically informative.

---

## 8. Direction is outside canonical TTF v0.2

TTF's primary turnover score uses within-species dissimilarity rank and is unsigned.

Therefore a shared boundary can be compatible with opposite trait directions in different species.

A pollination/selfing/island-syndrome claim needs a separate signed endpoint layer, for example a prospectively oriented trait contrast or architecture-state label.

Recommended inference order:

```text
TTF Lane U: where does turnover recur?
        |
        v
signed layer: which side changes toward which state?
        |
        v
PAYOFF / individual-system receipt: why is that direction favored?
```

Do not upgrade unsigned TTF sharedness into directional syndrome convergence.

---

## 9. Relation to TTF synthetic calibration axes

TTF separates shared fraction `pi` from within-species signal amplitude `A`.

Under the PAYOFF bridge:

```text
pi
~ fraction of species whose spatial transition geometry is sufficiently aligned
  to be transferable;

A
~ observable trait-turnover strength conditional on a transition.
```

These are analogies, not parameter identities.

In particular:

```text
A != |phi|
A != |eta|
pi != fraction of species with an identical causal mechanism.
```

The mandatory adversarial null `pi=0, A>0` remains exactly relevant: every species may have a strong private payoff/trait transition at a different location.

---

## 10. Prospective combined tests

A strong workflow freezes the following before opening evaluation traits:

```text
H1  canonical TTF transfer: T>0;

H2  PAYOFF increment:
    Delta_P|GE=T_GEP-T_GE>0;

H3  signed direction, if claimed:
    held-out transition direction agrees with the declared syndrome direction;

H4  individual-system anchor:
    direct or historical evidence supports the proposed module/payoff interpretation.
```

The four targets are distinct and should be reported separately.

---

## 11. Falsification interpretation

Useful negative combinations include:

```text
T>0 but Delta_P|GE<=0
    transferable turnover exists, but the proposed PAYOFF mechanism does not
    add predictive value beyond the baseline.

T<=0 but Delta_P|GE>0
    the generic learned field is not broadly shared, but a specified
    mechanism predictor transfers better than the trait-derived field.

T>0 and Delta_P|GE>0 but signed direction fails
    shared transition geometry exists without a common syndrome direction.

anchor fails
    cross-species spatial transfer cannot rescue the proposed individual-system
    mechanism interpretation.
```

---

## 12. Claim ceiling

A successful PAYOFF-TTF bridge can support a statement of the form:

> An independently specified architecture-payoff projection adds out-of-species predictive skill for within-species trait turnover beyond declared geographic/environmental baselines.

It does not by itself prove:

```text
a universal causal barrier;
identical mechanisms across species;
a common directional syndrome;
or direct PAYOFF parameter identification in held-out species.
```

Those require separate mechanistic and directional receipts.