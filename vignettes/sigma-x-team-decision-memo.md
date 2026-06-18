# Decision memo: the omnibus balance test after matching

**To:** Ben, Mark, Josh
**From:** Jake (with Claude)
**Re:** `devel-sigma-x-omnibus` --- what the sigma_x work does and does not buy us, and the decisions we have reached
**Date:** 2026-06-14

This is the short version. The long working memo is
`devel-sigma-x-omnibus-memo.qmd`; the one table below is reproducible with
`vignettes/three-test-comparison.R`.

## The question

The branch added an alternative omnibus statistic `T = d' Sigma_x^{-1} d` that
standardizes the treatment-control difference vector `d` by a covariance of the
covariates, instead of HB08's permutation covariance `V_d`. The motivation was
Rosenbaum's concern (relayed by Ben): under tight matching the released test
rejects on imbalances that are substantively trivial. We wanted to know whether
the new statistic actually answers that concern, and how to set its defaults.

## The one table

Memo's Rosenbaum fixture: 4 strata, 2 treated each, a fixed gap of 0.5 on `x1`,
within-stratum residual scale `shrink` shrinking left to right (tighter
matching). p-values for the omnibus balance test computed five ways.

```
                 shrink=0.5  shrink=0.2  shrink=0.05  shrink=0.01
HB08 (chi-sq)        0.262      0.015        0.003        0.003
HB08 (within-perm)   0.287      0.006        0.002        0.002
FixedCov (within)    0.358      0.006        0.002        0.002
Pooled  (within)     0.287      0.006        0.002        0.002
Complete-rand        0.937      0.899        0.898        0.898
```

Anchor for magnitude: at `shrink=0.05` the `x1` imbalance is 0.49 in raw units,
which is 0.11 of a full-sample SD --- small --- yet the within-block test gives
p = 0.002 while the complete-randomization benchmark gives p = 0.90.

## What the table says

1. **The over-rejection is real and lives in the REFERENCE, not the metric.**
   As matching tightens, the within-block test drives a 0.11-SD imbalance to
   p = 0.002. That is because the within-block permutation reference `V_d`
   shrinks with the matching.

2. **A fixed metric does not change the within-block test.** `FixedCov` (full-
   sample `cov(X)`) and `Pooled` (the branch default) give the *same* p-value as
   HB08 to simulation precision once matching is tight (`shrink <= 0.2`). They
   differ from HB08 only when stratification is loose (`shrink=0.5`: 0.358 vs
   0.287) --- the metric helps least exactly where Rosenbaum's concern lives.
   (In one dimension the metric cancels exactly: both the observed statistic and
   its reference scale by the same constant.) **Consequence: the sigma_x metric
   is not a different test. Its value is as a descriptive natural-scale
   magnitude (an omnibus Mahalanobis imbalance), not as a p-value.**

3. **Only changing the reference makes the test insensitive to tightness.** The
   complete-randomization benchmark stays near 0.90 regardless of `shrink`. But
   it gets there by ignoring the blocks --- it answers "is balance better than a
   coin-flip experiment on these units?" (a matching-quality benchmark), not
   "does within-block as-if-randomization hold?" (the identifying assumption the
   matched analysis relies on).

4. **The R<1e4 warning is validated.** Here R = prod choose(4,2) = 1296, and the
   chi-square p (0.015) overstates the within-permutation p (0.006) at
   `shrink=0.2` --- the small-support regime the warning is meant to catch.

**Reading of Rosenbaum 2025 (Ch. 6, pp. 149-161).** He proposes exactly the
complete-randomization benchmark and defends ignoring blocks by definition:
"covariate balance refers to the distribution of covariates in treated and
control groups, not to who is paired with whom" (p. 151). He never argues *why*
the principle of attending to pairs is wrong; our read is that he is conflating
a test-statistic/reference pathology (column-by-column above) with the principle
of conditioning on the design. He himself flags the disagreement (Problem 6.3,
p. 159: "not everyone agrees, and you should form your own opinion"). His
purpose is a stopping rule ("when to stop improving the match"), not a test of
the identifying assumption; his implementation uses per-covariate tests combined
by minimum-P / truncated product / Fisher (the `iTOS::evalBal` function), not a
d^2 omnibus.

**Our interpretation.** The within-block rejection is arguably *correct*: a
consistent same-sign gap across strata genuinely violates within-block
as-if-randomization. The real issue is significance-vs-magnitude --- a 0.11-SD
imbalance is detectable but may be substantively trivial. The fix that keeps the
design is to report the natural-scale magnitude alongside the p-value (which the
package already does per covariate), not to swap the test for one that ignores
the blocks.

## Decisions reached

1. **Primary test stays HB08 / within-block.** It tests the identifying
   assumption of the matched analysis. `balanceTest()`'s default is unchanged.

2. **RCT vs matched.** When within-block randomization is a *fact* (an RCT,
   including a matched-pair randomized experiment), `V_d` / HB08 is correct and
   the tight-matching sensitivity is a feature --- it flags broken protocols,
   buggy assignment, differential attrition. The discussion above is entirely
   about the case where within-block randomization is an *assumption* (matched
   observational study). Trigger: "is the randomization real?", not "did you
   match?"

3. **Default null backend for the fixed-Sigma path = `satterthwaite_finite`.**
   Warn, and suggest `simulate` (or `exact`), when
   R = prod_s choose(n_s, n1_s) < 1e4. Document `simulate` as the escape hatch.

4. **Third-moment skewness check deferred.** It would need a from-scratch
   Finucan sixth-order derivation (no code exists to port; the
   `i113-highermoments` branch stops at the variance). Runtime would be fine
   (analytic, ~10-50 ms), but the derivation effort is high. Revisit only if
   someone invests the algebra.

5. **Complete-randomization benchmark: ship as opt-in calibration.** Frame it as
   a matching-quality benchmark, not a reject/accept test (a large p-value does
   not certify the identifying assumption). Attribute the idea faithfully to
   Rosenbaum 2025 Ch. 6 and the works it cites (Brumberg et al. 2024,
   Hansen-Bowers 2008, Pimentel et al. 2015, Yu 2021), and state in the docs
   that applying it to our d^2 omnibus is in the spirit of, not identical to,
   his per-covariate implementation.

## Open questions for the team

1. Given that the fixed metric does not change the within-block *test*, do we
   keep the sigma_x statistic at all --- repositioned as an omnibus *descriptive*
   magnitude (natural-scale Mahalanobis imbalance) rather than a rival test? If
   yes, the default-Sigma choice (full-sample `cov(X)` vs within-stratum-pooled)
   matters only for that descriptive number, and full-sample `cov(X)` is the
   natural pick.

2. Do we want an equivalence-testing path (test "imbalance within tolerance
   delta") as the principled way to keep the design *and* avoid rejecting
   trivial imbalances? Set aside for now, but it is the clean answer to
   significance-vs-magnitude.

3. Confirm the complete-randomization benchmark ships as calibration only, with
   no reject/accept language anywhere in the output or docs.
