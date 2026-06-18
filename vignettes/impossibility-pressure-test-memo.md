# Pressure-testing the balance-calibration impossibility

To: Jake Bowers, Ben Hansen
Date: 2026-06-17

## 1. The question

When we match well, our omnibus balance test can turn against us. Tighten a match so that treated and control units inside each set look nearly identical on the covariates, and the d^2 test of Hansen and Bowers (2008) can return a small p-value --- the design looks alarmingly imbalanced --- even though the actual covariate gap between treated and control never changed. The reason is that d^2 measures the realized mean difference against the within-set randomization variance V_d, and a tight match shrinks V_d toward zero. A fixed gap divided by a vanishing yardstick looks enormous. So we ask one question: is there a single number, computed only from the realized matched design (the within-set randomization), that reports how imbalanced a design is, rises when the real imbalance rises, and does NOT climb toward its maximum merely because the match got tighter? Call such a number a design-internal, magnitude-sensitive, non-collapsing calibration. We ran an adversarial search to find one or to show that none exists.

## 2. The dimensional argument (the conjecture under test)

Here is the reasoning that says no such number exists. A design-internal calibration sees only two things: the aligned mean difference dbar and the within-set randomization covariance V_d. Anything it reports is some function f(dbar, V_d). Now ask what happens under a change of units --- rescale the covariates, or equivalently shrink the match. A calibration that does not care about the arbitrary scale of the covariates must be scale-free (degree zero), and a degree-zero function of (dbar, V_d) can depend on them only through the combination dbar / sqrt(V_d), the standardized imbalance. As the match tightens, V_d -> 0 with dbar fixed, so dbar / sqrt(V_d) -> infinity. Any scale-free design-internal number therefore diverges to its maximum: it collapses. The only way to stay bounded as V_d -> 0 is to not depend on V_d at all --- but a function of dbar alone is a bare magnitude (a ruler with no reference null), not a calibration.

A note on the limit. The divergence dbar / sqrt(V_d) -> infinity is asymptotic; in the finite harness sqrt(V_d) has a finite-sample floor (the group-center jitter floors it), so the standardized imbalance climbs (1.52 -> 7.67) but does not literally diverge. This is consistent with treating the dimensional argument as a heuristic rather than a theorem.

State this precisely. CONJECTURE (the impossibility): no single scalar can be design-internal, magnitude-sensitive, and non-collapsing at the same time. The dimensional argument above is a heuristic proof, not a theorem --- it assumes the calibration is exactly degree-zero in scale and a smooth function of (dbar, V_d), and "magnitude-sensitive" and "non-collapsing" are stated as limiting behaviors rather than as a closed inequality. We treated it as a conjecture and tried to break it.

## 3. The adversarial search

Nine inventors each took a distinct strategy and built the most honest design-internal calibration that strategy allows, then ran it through a shared harness: a collapse curve (gap held fixed at 0.2, match tightened from shrink 0.6 down to 0.012) and a magnitude curve (match held fixed, gap grown from 0 to 1). The strategies were:

1. d^2 divided by its own within-set re-randomization null (a ratio built so V_d would cancel).
2. A balance number from the generalized eigenvalues of within-stratum vs total covariance.
3. A studentized omnibus using across-stratum spread of per-set differences instead of V_d.
4. A Bayesian / empirical-Bayes shrinkage posterior that the standardized within-set imbalance exceeds a threshold.
5. An e-value / betting martingale accumulated across matched sets.
6. An aligned-rank Wilcoxon statistic at the matched sets, calibrated by within-set re-randomization.
7. A within-strata adjusted R^2 of treatment on covariates, corrected for the K/(n-S) chance floor.
8. An explicit tightness penalty: subtract lambda * log(V_d) from the collapsing within-set z.
9. A between-set-standardized within-set re-randomization percentile (the wildcard).

Result: nine claimed escapes, zero survived. Every candidate that was honestly design-internal and magnitude-sensitive collapsed; every candidate that did not collapse had either dropped V_d (becoming a bare magnitude) or normalized the magnitude away (becoming magnitude-blind). No breaker found a flaw that rescued any candidate, and no survivor remained.

The instructive failures:

- The ratio (idx 1) was the most direct attempt to cancel V_d. It failed because the within-set null distribution of d^2 is scale-free --- a chi-square with df = rank(V_d), whose mean stayed pinned near 3 and whose 95th percentile stayed near 7.4 across every shrink level. Dividing a degree -2 observed quantity by a degree-0 reference leaves a degree -2 quantity. The V_d cancels inside the reference, never in the numerator. The raw ratio climbed 1.60 -> 5.63 -> 17.30 -> 19.22 as the match tightened at a fixed gap. (These exact figures belong to an idx-1 construction not in the shared harness; the direction and mechanism --- the ratio climbs as the match tightens because V_d cancels only in the reference --- are confirmed in the harness, the specific numbers are not.)

- The tightness penalty (idx 8) is the cleanest illustration of the bind. The alarm is log(max d^2) - log(V_d), and subtracting lambda * log(V_d) cancels the divergence only at lambda = 1 exactly. Any lambda < 1 under-cancels and still collapses; at lambda = 1 the penalty deletes V_d entirely and leaves log(max d^2), a bare magnitude. To pick a non-trivial cancellation point you need a reference tightness that does not shrink with the match --- and the only such anchor is a fixed pre-match scale, which is external by definition.

The closest near-miss was the W-vs-Tot spectrum (idx 2), the only inventor who built and ran a non-collapsing design-internal sibling. Its "Candidate B" --- a degree-zero weighted average of surviving variance fractions along the imbalance direction --- is design-internal and does not collapse; the alarm even falls as the match tightens (0.062 -> 0.009 -> 0.0005). But its magnitude curve is flat noise (0.0008, 0.0011, 0.0003, 0.0006, 0.0015): it achieved non-collapse by normalizing the magnitude away, so it depends only on the direction of dbar, not its size. That is exactly the bounded degree-zero function the dimensional argument says must ignore V_d. The near-miss broke on the magnitude axis, precisely where the conjecture predicts.

Zero survivors out of nine adversarial attempts is strong support for the conjecture, not a proof of it. The search is bounded by the strategies the nine inventors imagined; a tenth strategy could in principle escape. But every failure broke for the reason the dimensional argument names, which raises our confidence that the argument identifies the real obstacle rather than an accident of these nine designs.

## 4. The escape is an external reference

If the obstacle is V_d -> 0, the escape is to stop dividing by V_d --- to compare the realized design against a fixed reference that does not shrink with the match. Jake proposed comparing the matched design to a whole-pool complete randomization (CRE): a coin flip over all units, ignoring the strata. We tested four versions on six criteria.

| Reference | Collapses? | Magnitude-sensitive? | Floor-blind? | Valid percentile? | Unbuilt-design cost? | Oracle-free works? |
|---|---|---|---|---|---|---|
| Pool-CRE, max\|SMD\| (raw) | no | yes | YES | no (conservative) | yes | yes |
| Pool-CRE, fixed-Sigma Mahalanobis | no | no | YES | no | yes | yes |
| Pool-CRE + covariance adjustment (oracle group) | no | yes | no | no | yes | yes |
| Pool-CRE + covariance adjustment (oracle-free k-means) | no | yes | no | YES | yes | yes |

Read the table top to bottom. The raw whole-pool CRE with max|SMD| does escape collapse --- at fixed gap, the alarm stays flat near 0.005-0.018 as the match tightens, and on the harder, less-separated DGP the alarm even falls as the match tightens (0.266 -> 0.033), the opposite of collapse, so the no-collapse result is not an artifact of well-separated groups. But it is NOT a valid percentile as coded: the function `cal_maxsmd_poolCRE` computes the observed statistic with within-set centering while computing the reference with pool centering, so they are different statistics and the alarm is pinned near 0.000 under its own reference (mean 0.000, KS D 0.992 vs uniform), not uniform. This is the same observed-vs-reference mismatch flagged for rows 2-3 below; the alarm is conservative and downward-biased rather than an honest tail probability. (An apples-to-apples pool/pool version --- observed AND reference both pool-centered --- IS genuinely uniform, mean 0.495, KS D 0.040, but that is a different contender from the set-vs-pool function in this row.)

And it is so floor-blind that it clears almost every design. The whole-pool reference scrambles treatment across groups that nobody randomized over (group centers 0, 8, 16, 24 in the test data), so the reference max|SMD| distribution is enormous (median 0.169) relative to any realistic matched-set imbalance (observed 0.034). The alarm crosses 0.05 only at a residual gap near 0.5 SD-direction units and crosses 0.50 only near 2.7. A design has to be badly broken before this reference notices.

The fixed-Sigma Mahalanobis version is worse on two counts: it is even more floor-blind (one between-group eigenvalue of order 250 dominates the pool covariance --- 262.8 in the memo's seed, near 245 in another, the same order of magnitude --- so the reference distances are huge and the observed within-set contrast sits permanently in the left tail), and it is an invalid percentile as coded --- the observed statistic uses within-set centering while the reference uses pool centering, so they are different statistics and the alarm sits at 0.000, not 0.5, under its own null.

What fixes the floor-blindness is a covariance adjustment: residualize the covariates on the coarse group structure before computing max|SMD| against the pool-CRE null. With the oracle group label, this clears collapse and floor-blindness (magnitude rises 0.002 at gap 0 to 0.926 at gap 1) but is not a valid percentile as implemented, because the observed statistic is set-adjusted while the reference is pool-only --- the non-collapse comes from comparing two different statistics, which yields a conservative, downward-biased number rather than an honest tail probability.

The oracle-free caveat is the one we can report as resolved. The covariance adjustment does NOT need the true group label. Recovering the coarse structure unsupervised by k-means on the observed covariates (k-selection recovered the coarse group count and did not chase the 20 matched sets) approximately reproduces the oracle-adjusted alarm (agreement depends on group separation and the k-selection rule), AND --- when the observed and reference statistics are matched apples-to-apples --- gives a valid percentile (mean alarm 0.510 at gap 0, KS-vs-uniform p 0.057). But two oracle-free choices fail, and the difference is not obvious from the outside. Residualizing each covariate on the others (the naive move) is non-collapsing but FLOOR-BLIND, because the common-direction gap u = (1,1,1)/sqrt(3) is collinear and gets absorbed. Residualizing on matched-set fixed effects COLLAPSES (alarm climbs toward its maximum, 0.832 -> 1.000, read off the curve) and is floor-blind. Note that the collapse verdict for set-FE rests on the climb-to-maximum read directly off the curve, not on the harness `collapse_curve(...)$collapsed` boolean: because the loose-match value is already near 0.9, the rise falls under the harness's `a_tight - a_loose > 0.1` threshold and the boolean reports FALSE. That boolean has a false-negative for any contender that starts high at the loose match. So oracle-free adjustment works, but only when the recovered structure captures the coarse confounding the match exploited without absorbing the treatment-direction signal. That is a real condition on the analyst's choice, not an automatic gain.

## 5. Recommendation and honest scope

The search supports a three-way bind. Of these three properties --- a single calibrated omnibus number, a design-internal reference, and non-collapse as the match tightens --- you can have at most two. The nine inventors confirm you cannot keep all three with a within-set null. The external references show the only escape from collapse requires a fixed reference that does not shrink, and that escape always incurs the unbuilt-design cost: a perfectly balanced matched design (gap 0) is scored against a whole-pool randomization it would never have used and reads as beating more than 99 percent of those draws (alarm 0.004). You are no longer asking "is this design balanced given how it was built"; you are asking "is this design more balanced than a coarser randomization nobody ran."

What to report:

- Lead with a bare descriptive magnitude --- observed max|SMD| on the residualized covariates --- with a substantive gate the reader can check. This is Ruler A: collapse-immune by construction, magnitude-sensitive, but it requires you to choose the cutoff c, and the same design flips flagged-or-cleared on c alone. It is not a calibration and should not be dressed as one.
- Add the structure-adjusted pool-CRE percentile as a SEPARATE number, with the covariance adjustment recovered from the observed covariates by k-means and pre-specified. Label it honestly: the percentile is against an unbuilt whole-pool randomization, its power depends on the recovered structure capturing the same variation the matching used, and (in the set-adjusted-vs-pool form) it is conservative rather than exactly uniform.
- Do not present any single number as a design-internal calibration that escapes the impossibility. None did.

Two scope limits. First, everything here concerns OBSERVED covariates. A design that beats every randomization reference on the covariates we measured says nothing about the covariates we did not --- this is balance reporting, not a defense of ignorability. Second, the impossibility is a conjecture supported by a bounded adversarial search and a heuristic dimensional argument, not a proved theorem. Turning the dimensional argument into a theorem --- stating the exact regularity and degree-zero conditions under which no f(dbar, V_d) can be magnitude-sensitive and bounded as V_d -> 0, and proving it --- is the open problem this memo does not solve.