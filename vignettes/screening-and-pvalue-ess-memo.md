# Screening near-balanced covariates, and why the omnibus p-vs-ESS plot is jagged

Temporary memo. Date: 2026-06-17. Companion to fixed-metric-collapse-memo.md.
Two design thoughts: (1) a screening step in balanceTest() that drops
near-exactly-balanced covariates from the omnibus while still describing them;
(2) why the omnibus p-value plotted against effective sample size is not smooth,
and why "collapse" is only one of the patterns.

## Engine facts that ground thought (1)

A screen already exists. In R/Design.R:991:

    ssvar <- diag(tcov)
    zero_variance <- (ssvar <= .Machine$double.eps)
    zstat <- ifelse(zero_variance, NA_real_, ssn / sqrt(ssvar))

Here ssvar is the per-covariate within-set null variance of the aggregated,
size-scaled treated-minus-control difference (the diagonal of the test-statistic
covariance tcov, i.e. the V_d analog). Covariates with ssvar at or below machine
epsilon are dropped: their univariate z is set to NA (so they leave the
univariate family too), and only the survivors enter the SVD pseudoinverse.

The pseudoinverse applies a SECOND, relative screen on singular values
(R/utils.R:393): a direction counts only if d > tol * d[1], with
tol = sqrt(machine eps) ~ 1.5e-8, i.e. relative to the largest singular value.

The pseudoinverse ERROR we hit in the extreme-case sim
(fixed-metric-collapse-memo.md, Case 1, eps <= ~1e-6) is the all-screened case:
the surviving submatrix is 0x0 and XtX_pseudoinv_sqrt() calls stop() at
R/utils.R:381.

So the proposal is not a new mechanism. It is: change the TOLERANCE on a screen
that already exists, degrade gracefully instead of erroring, and report the
screened covariates descriptively.

## (1) The screening step

### Screen on the design variance, never on the observed imbalance

ssvar (and tcov generally) is a function of the covariate values, the strata,
and the group sizes --- NOT of the realized treatment vector Z. For pairs,
V_d = (1/S^2) sum_s c_s c_s' depends only on the within-pair differences c_s, not
on who was treated. So screening on ssvar (or on the singular values of tcov) is
assignment-invariant: it cannot be gamed by peeking at the direction the
imbalance points. Screening instead on the observed zstat or dbar ("drop
covariates that look balanced") would snoop the very thing being tested and make
the p-value manipulable. Keep the screen on the design side. This is the property
that makes it legitimate rather than p-hacking.

### "Relative to something" = the complete-randomization variance

The natural denominator is the variance the SAME difference would have with no
stratification: sigma_pool^2 * (1/n1 + 1/n0). The ratio

    ssvar / ssvar_complete  =  fraction of the covariate's variance surviving stratification

goes to ~0 under exact matching. So the screen is "drop direction j if the design
removed more than (1 - tau) of its variance." This denominator is exactly the
quantity in MODE 1 of the failure-mode checklist; the screen and the MODE-1
diagnostic compute the same design ratio.

### Do not pick one tau --- profile it

The arbitrary-cutoff objection is the same defect the coarsening rule had, and it
has the same fix: report (df, d^2, p) as a function of tau. That profile IS the
resolution-profile idea applied to the screening tolerance. At tau near machine
epsilon you get today's behavior (the lone perturbed pair contributes its
spurious 1); as tau rises, directions drop out and df steps down. The shape of
that curve is the honest report, not a single number.

### Screen directions, not just raw covariates

The diagonal screen at line 991 catches a single exactly-stratified covariate.
Exact matching on an interaction or a coarsened combination produces a degenerate
LINEAR COMBINATION --- collinearity within strata --- which only the SVD-stage
relative tol catches. Apply the substantive tolerance at both stages and report
which DIRECTIONS (not just which named covariates) were dropped.

### Graceful degradation + descriptive-only reporting

Replace the stop() with "omnibus not computed: all directions balanced within
tolerance tau; see descriptive table." Screened covariates become
descriptive-only, excluded from the multiplicity family --- and the right
descriptive is exactly the percentile-of-within-set-differences summary, because
the screened covariates are precisely the ones where magnitude is the only
remaining story and the omnibus is structurally blind to it (degree-zero
homogeneity, see fixed-metric-collapse-memo.md). The screen and the percentile
summary are two halves of one coherent report: test the directions that vary,
describe the magnitude distribution for the rest.

### Framing / prior art

Declaring "balanced within tolerance tau" rather than testing "difference = 0" is
an equivalence-test move. Hartman and Hidalgo (AJPS 2018, "An Equivalence
Approach to Balance and Placebo Tests") argue balance assessment should be
equivalence testing, not difference testing; this screen is in that family but
used as a pre-filter on the omnibus rather than a wholesale replacement. VERIFY
that citation before it goes in writing.

### Caution

The screen changes what the omnibus TESTS --- it becomes "no systematic imbalance
among the non-degenerate directions" --- so the report must name what was
removed. With the profile, that is automatic.

## (2) Why p-vs-ESS is jagged, and why collapse is not the whole story

Collapse is one pattern, not a law. The clean statement: the omnibus has power
only against standardized SYSTEMATIC (directional) imbalance, scaled by ESS, and
ZERO power against magnitude.

Under an alternative with systematic per-comparison bias mu_d, d^2 is noncentral
chi-square with noncentrality

    lambda = mu_d' V_d^{-1} mu_d  ~  S * (standardized systematic bias)^2.

Scale every covariate difference by t and both mu_d and V_d^{1/2} scale by t, so
lambda is unchanged --- no power against magnitude (the degree-zero fact). But
lambda grows with S (with ESS), so power against a DIRECTIONAL alternative does
increase with sample size. That splits the world into two regimes:

- Regime A --- detectable systematic residual. The matching left a directional
  bias it could not remove. The omnibus behaves like an ordinary test: worse
  balance gives smaller p, improving balance raises p toward uniform, more ESS
  sharpens it. NOT collapsed. This is the regime where "p does not just sit
  there" examples live.

- Regime B --- already within-set-null, shrinking magnitude with null signs. p is
  flat regardless of magnitude. This is the collapse, and it is where the
  extreme cases were constructed to sit.

So "something else is going on" in a p-vs-ESS plot is three non-smooth effects
stacked together:

1. df = rank(V_d) is an integer and JUMPS as covariates (or directions) cross the
   degeneracy threshold while you prune. P(chisq_df >= d^2) steps discontinuously
   every time df changes --- the same mechanism as the screening discussion.

2. d^2 is direction-driven, not magnitude-driven. As you prune units the residual
   DIRECTION swings around, so d^2 jitters even when the magnitude of imbalance is
   falling monotonically.

3. The reference itself is ESS-dependent. The chi-square approximation to the
   discrete within-set permutation null degrades at small ESS (few effective
   strata, skewed/lumpy support), so the p-mapping moves even at fixed d^2.

Underneath all three: pruning to change ESS also changes the population, the
weights, and the standardized-bias content, so the curve is really p plotted
against a product lambda ~ ESS * bias^2(ESS) in which BOTH factors move as you
slide along the path, often in opposite directions. That alone guarantees
non-monotone, non-smooth behavior.

### Experiment to separate the effects

Hold the sample and strata fixed and vary one thing at a time:

  (a) a magnitude scalar t on a fixed residual    -> predict p flat (pure collapse)
  (b) sign-alignment fraction at fixed magnitude  -> predict p moves
  (c) ESS at fixed standardized systematic bias   -> predict p moves through lambda

Any real p-vs-ESS plot is a path through that 3-D space; the jaggedness resolves
into "df steps + direction swings + reference degradation."

## Open items / next steps

- Reproduce Jake's actual p-vs-ESS curve (need: which dataset, how the sequence of
  matchings was generated, whether ESS changed by pruning or by changing set
  sizes) and decompose it against the 3-D space above. More useful than inventing
  a fresh curve.
- Prototype the relative-tolerance + profile screen against R/Design.R:991 and
  plot the (df, d^2, p)-vs-tau curve on the Case 1 data.
