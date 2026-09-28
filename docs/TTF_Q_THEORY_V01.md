# TTF-Q theory v0.1 — nested information survival and calibrated evaluable envelopes

## Scope

TTF-Q characterizes a dyadic predictor before the biological response is
opened. The quantities below are properties of the predictor and intended
design matrix. They are not partitions of outcome variance.

Let (rinmathbb{R}^n) be the centered focal source-target relation on a fixed
dyad set. Let (F) be the exact source/target fixed-effect design and
(M_F) its residual-maker. Let (C) be a matrix of declared nuisance
predictors and define (C_F=M_FC).

The focal relation remaining after endpoint identity is

[
r_F=M_Fr.
]

The relation remaining after endpoint identity and controls is

[
r_C=M_{C_F}r_F,
]

where (M_{C_F}) is the orthogonal residual-maker for the columns of (C_F).

TTF-Q reports squared-norm survival fractions because the predictor is fixed
before any response is opened.

## Proposition 1 — exact two-stage factorization

Define

[
S_F=rac{|r_F|^2}{|r|^2}
]

as endpoint-identity survival and

[
S_C=rac{|r_C|^2}{|r_F|^2}
]

as control survival conditional on endpoint identity.

Then total identifiable relational information is

[
S_{mathrm{total}}
=rac{|r_C|^2}{|r|^2}
=S_FS_C.
]

This is an algebraic identity. It is useful because it localizes information
loss: (1-S_F) is loss to endpoint identity, while (1-S_C) is the fraction
of within-endpoint relation information removed by the declared controls.

The implementation exposes these as

- `source_target_fe_retained_variance_fraction`;
- `control_unique_variance_fraction_after_fe`;
- `total_unique_variance_fraction`.

## Proposition 2 — nested-control survival

Let (C_1) be a baseline control set and let
(C_2=[C_1,G]) add one or more new controls (G) on the same dyads. Let the
corresponding residual relations be (r_1) and (r_2).

Define incremental survival

[
S_{2|1}=rac{|r_2|^2}{|r_1|^2}.
]

Then

[
S_{mathrm{total},2}
=
S_{mathrm{total},1}S_{2|1}.
]

For the climate case, (G) is continuous directed geographic co-opportunity.
Thus

[
1-S_{2|1}
]

is the exact fraction of previously unique C1 historical relation information
attributable to that additional geographic control under the declared linear
projection.

## Proposition 3 — monotonicity under nested controls

For nested control spaces,

[
operatorname{col}(M_FC_1)
subseteq
operatorname{col}(M_FC_2),
]

orthogonal projection implies

[
0le |r_2|^2le |r_1|^2.
]

Therefore

[
0le S_{2|1}le 1
]

and adding controls cannot increase TTF-Q residual information except for
floating-point tolerance.

The implementation checks this identity explicitly in the C1-to-C2
decomposition.

## Proposition 4 — scale invariance

For any nonzero scalar (a),

[
S_F(ar)=S_F(r),qquad
S_C(ar)=S_C(r),qquad
S_{mathrm{total}}(ar)=S_{mathrm{total}}(r).
]

Likewise, nonsingular rescaling of individual nuisance columns leaves their
column space unchanged and therefore leaves the residual relation unchanged.

Thus information-survival fractions describe geometry of the declared
predictor spaces rather than arbitrary measurement units.

## Signal concentration is a separate object

After the final declared control projection, let

[
w_{st}=r_{C,st}^2
]

and normalize (pi_{st}=w_{st}/sum w_{st}).

TTF-Q reports the inverse-Herfindahl signal concentration

[
N_{mathrm{dyad}}^{mathrm{signal}}
=rac{1}{sum_{st}pi_{st}^2}.
]

For source (s), aggregate
(pi_s=sum_tpi_{st}), giving

[
N_{mathrm{source}}^{mathrm{signal}}
=rac{1}{sum_spi_s^2},
]

and analogously for targets.

These values are bounded concentration summaries. They are **not effective
sample sizes**, do not determine standard errors, and do not substitute for
dependence-aware inference.

The known-truth benchmark is deliberately constructed so that unique
information can be held fixed while endpoint concentration changes. This
demonstrates empirically that information amount and information concentration
are distinct.

## Null calibration precedes detectability

Let (A) index a declared nuisance-heterogeneity regime and let
(p_0(A)) be the synthetic rejection probability under focal effect
(eta=0).

TTF-Q first evaluates a finite Monte Carlo estimate with Wilson interval
([L_0(A),U_0(A)]). For a declared Type-I ceiling (	au), define

[
K(A)=mathbf 1{U_0(A)le	au}.
]

The raw grid minimum detectable effect is

[
D_{mathrm{raw}}(A)
=
min{eta>0:L_eta(A)ge pi^*},
]

where (L_eta(A)) is the Wilson lower bound on detection probability and
(pi^*) is the target power.

The **evaluable** MDE is

[
D_{mathrm{eval}}(A)=
egin{cases}
D_{mathrm{raw}}(A), & K(A)=1,\
	ext{undefined}, & K(A)=0.
end{cases}
]

This ordering is essential. A high rejection probability is not evidence of
detectability if the same estimator is anti-conservative under the null.

## Evaluable envelope

For a fixed predictor geometry, TTF-Q defines the response-blind evaluable
envelope as the joint object

[
mathcal E=
left(
S_F,
S_C,
S_{2|1},ldots,
N_{mathrm{dyad}}^{mathrm{signal}},
N_{mathrm{source}}^{mathrm{signal}},
N_{mathrm{target}}^{mathrm{signal}},
{K(A),D_{mathrm{eval}}(A)}_{Ainmathcal A}
ight).
]

Optional measurement-repeatability quantities may be appended when
response-blind resampling of the ecological relation is available.

The envelope intentionally does not collapse to one score. The known-truth
benchmarks show why:

- endpoint loss and nuisance redundancy can reduce identifiable information
  and worsen calibrated MDEs;
- concentrated signal can leave the total unique fraction unchanged while
  invalidating null calibration;
- therefore information amount, concentration, calibration and conditional
  detectability answer different design questions.

## Known-truth validation

The frozen decomposition benchmark specifies (S_F), (S_C), and nested
geographic survival before results are generated. Across 243 scenarios the
implementation recovered all prescribed fractions with maximum absolute error
approximately (3.33	imes10^{-16}).

The calibrated detectability benchmark then holds the legacy-style reference
cell fixed at (eta=.03,A=3). Distinct constructed designs can share the same
negative binary classification while having different calibrated MDE vectors,
or can be uncalibrated for an entirely different reason. This is the direct
methodological motivation for retaining (mathcal E) rather than a single
qualification label.

## Interpretation boundary

None of these identities implies that the focal ecological relation predicts a
biological outcome. TTF-Q qualifies the information and inferential geometry
available **before** that outcome is opened.
