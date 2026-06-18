# Does a fixed-metric, center-then-Mahalanobis balance statistic collapse?

Temporary memo. Date: 2026-06-17. Context: the balance-calibration
collapse / scale-invariance impossibility (see the
collapse-calibration-strengthened-impossibility and
balance-calibration-failure-mode-checklist notes).

## Plain-language summary (read this first)

The thing that decides whether the test rewards a better matching is what you
compare your statistic against. It is not the Mahalanobis rescaling, and it is
not whether you standardize before or after.

You compute one number: the imbalance after matching, measured with a fixed
covariance metric. Call it the observed statistic. To say whether that number is
"good," you need something to compare it to. You have two choices, and they give
opposite results.

Choice A: compare to a chi-square. A chi-square curve does not change when you
change the matching --- it is the same curve no matter how good or bad your
matching is. So if a better matching makes your observed number smaller, it
lands further down the same fixed curve and you get a bigger p-value. The number
does track matching quality. The re-inflation intuition is exactly right here:
the Mahalanobis metric blows the small differences back up, and because the
comparison curve is sitting still, that inflation shows up in the answer. The
catch: a chi-square is the right comparison curve only if you ignore the blocks
and pretend the whole sample was randomized at once. It is not the right curve
for "randomize within the matched sets." So this version answers the wrong
question (complete randomization, Rosenbaum's marginal benchmark) and, used as a
within-block test, has the wrong error rate. It rewards better matching only
because it threw the blocks away.

Choice B: compare to within-set re-randomization. Now you shuffle treatment
labels within each matched set, recompute the same statistic many times, and ask
where your observed number falls. The problem: when the matching is finer, the
units inside each set are closer together, so EVERY reshuffle produces a small
imbalance too. The comparison distribution shrinks right alongside your observed
number. A better matching gives you a smaller observed value AND a smaller pile
of reshuffled values, and your observed value sits in the same relative spot in
the pile. In fact, if the matching has done its job (treatment looks
as-if-random within sets), your observed value is just one ordinary draw from
that pile, so its p-value is uniform --- no matter how fine or coarse the
matching is. So this version cannot tell a better matching from a worse one. And
the Mahalanobis metric does not save you, because you apply the same metric to
your observed value and to every reshuffled value: whatever inflation it does to
your number, it does identically to the whole pile, so it does not move where
your number sits. The metric changes WHICH directions of imbalance count, but it
cannot make the pile stop shrinking, because the pile's size is set by how tight
the matched sets are, not by your choice of metric.

So the two fixed-metric ideas are not two flavors of one thing. They differ in
one place only --- the comparison distribution --- and that one difference is the
whole story:

- The chi-square version (2a) rewards better matching, but only by ignoring the
  blocks and answering the complete-randomization question with the wrong error
  rate.
- The within-set permutation version (2b) honestly respects the blocks and has
  the correct error rate, but cannot reward better matching, because its
  comparison pile shrinks with the matching.

You cannot have both at once. That is the impossibility, in plain terms: a
comparison that respects your blocked design shrinks as the design improves and
so cannot grade the improvement; a comparison that grades the improvement has to
be one that ignores the blocks.

## The question

Suppose we have two matchings, one with lower delta (better balance) than the
other, and we build a test statistic in this order:

1. Subtract off the covariate means of the strata (center within strata --- NOT
   divide by within-stratum variance).
2. Compute a Mahalanobis-type quadratic form using a FIXED covariance metric.
3. Summarize.
4. Use within-strata randomization to generate a null.

The intuition to test: the Mahalanobis distance rescales all variables, so even
if the centered differences become small, the metric should re-inflate them.
Does the collapse come from the ORDER in which we standardize --- specifically,
does $d^2$ collapse because it standardizes LAST (by the null variance), and
could standardizing FIRST by a fixed metric escape it?

## Two clarifications that sharpen the construction

Clarification 1. Step 1 is centering only (subtract stratum means), not
within-stratum studentizing. And we are not taking a nonnegative distance
WITHIN each stratum and then summarizing. The sign-loss objection (a per-stratum
norm is swap-invariant, so within-pair randomization has nothing to act on, and
the null becomes a point mass for matched pairs) applies only to the
norm-within-stratum-first construction, which is not the intended one.

The power-preserving version of the idea is: center within strata, aggregate the
SIGNED within-stratum mean differences across strata, and apply a single
fixed-metric quadratic form at the end. That is exactly the two fixed-metric
variants below.

Residual caution on "summarize via median/mean": if you literally median the
per-stratum Mahalanobis DISTANCES, you are back to nonnegative-per-stratum and
you reintroduce the sign loss (robust, but sign-blind and low-power). The
version that keeps the permutation signal is aggregate-signed-then-one-norm,
which has no natural "median" form. The robustness a median would buy trades
directly against the sign information the permutation null needs.

Clarification 2. The two fixed-metric ideas already considered and "shot down":

- (2a) Compute Var-Cov of X using all units and covariates BEFORE matching. Do
  the matching (may drop units, certainly stratifies). Compute the d-statistics.
  Compute $d^2$ with the fixed Var-Cov and compare to a chi-square (large
  sample).

- (2b) Same fixed pre-matching Var-Cov. Do the matching. Build the DISTRIBUTION
  of $d^2$ by permuting within sets to get the d-statistics, each time using the
  SAME fixed Var-Cov to get $d^2$, and read a p-value off the permutations
  rather than the large-sample chi-square.

## The resolution

(2a) and (2b) differ in exactly ONE thing --- the reference distribution --- and
that single difference is the entire trilemma boundary. Neither the fixed metric
nor the order of standardization is the lever.

### (2b): fixed Sigma + within-set permutation reference --- valid test, collapses as a quality score

This is an exact permutation test. Under the within-set as-if-random null on the
matched covariates, $D_\text{fixed,obs}$ is LITERALLY a draw from its own
permutation distribution, so its permutation p-value is Uniform(0,1) by
construction --- at EVERY value of delta, for ANY fixed Sigma. A finer matching
that achieves within-set randomness and a coarser matching that achieves
within-set randomness on its looser sets both return $U(0,1)$. The p-value
carries zero information about delta. That is the collapse, stated as sharply as
it can be: when the null holds, the p-value is uniform by construction, so it
cannot rank designs by balance quality.

Why the fixed metric does not rescue this. In (2b), Sigma sits inside the
SCORING function, and the same Sigma scores the observed value and every
permutation draw. So the overall SCALE of the imbalance cancels in the rank
comparison --- it hits $D_\text{fixed,obs}$ and the reference identically. What
does NOT cancel is direction: Sigma reweights which directions of imbalance
count, so it changes WHAT is tested. But the reference's scale is carried
entirely by $V_d$, the within-set permutation covariance of the aggregated mean
difference, and $V_d$ shrinks with the matching ($V_d \sim \delta$). No choice
of fixed scoring metric can stop the reference from shrinking, because the metric
is a property of the score, not of the reference.

This is where the re-inflation intuition nets out: it is correct about the
numerator AND equally correct about the denominator. Inflating both sides of a
rank comparison by the same fixed linear map leaves the comparison's
delta-dependence exactly where it was.

### (2a): fixed Sigma + chi-square reference --- does NOT collapse, and that is the problem

The chi-square$_k$ reference is DESIGN-BLIND: it does not depend on the matching
at all. So as the matching improves, $\bar d_\text{obs}$ shrinks,
$D_\text{fixed,obs} = \bar d_\text{obs}' \Sigma^{-1} \bar d_\text{obs}$ shrinks
against a fixed yardstick, and the p-value rises. It DOES track delta. The
Mahalanobis inflation shows up undiluted here because there is no permutation
reference getting inflated alongside it. The intuition is exactly right for
(2a).

But that delta-tracking is bought with a design-blind reference, which has two
faces:

- As a within-set TEST it has the wrong size. $D_\text{fixed}$ is not
  chi-square$_k$ under the within-set null; it is a weighted sum of
  chi-square-ones, $\sum_j \lambda_j Z_j^2$, with $\lambda_j$ the eigenvalues of
  $\Sigma^{-1/2} V_d \Sigma^{-1/2}$. Comparing it to chi-square$_k$ is the wrong
  reference, so the size is off.

- Read charitably as a SCORE, the chi-square label is cosmetic and what you are
  really reporting is the fixed-Sigma length of the imbalance --- a bare
  magnitude. With the proper variance scaling, chi-square$_k$ is exactly the
  COMPLETE-RANDOMIZATION (block-ignoring) benchmark: Rosenbaum's marginal
  reference, the very thing the design was supposed to respect. This is the
  "forbidden = complete-rand via pooled scale" corner.

## Epitaphs (they were shot down for DIFFERENT reasons)

- (2a) does not collapse only because it silently ignores the blocks --- it
  answers the complete-randomization question and has the wrong size for the
  within-block design. Forbidden / bare-magnitude corner.

- (2b) honestly respects the blocks and has correct size, but collapses as a
  balance-quality score.

## Answer to "is it the order of standardization?"

No. What decides the answer is the comparison distribution, not where the
standardization sits or which fixed metric you pick. (2a) and (2b) hold the metric, the centering, and the
order all fixed and identical; they differ ONLY in
design-blind-reference vs. design-respecting-reference, and that one switch flips
you between "tracks delta but ignores blocks" and "respects blocks but
collapses." That is the trilemma made concrete in two lines of the same
construction.

The Mahalanobis-inflation instinct is right about the mechanism, and that is
exactly why it is not free: the inflation produces delta-tracking ONLY when
paired with a design-blind reference (2a). Pair it with the honest within-set
reference (2b) and the inflation cancels in scale. The impossibility is just the
statement that you cannot get the inflation to track delta AND keep the
design-respecting reference at the same time.

## Extreme-case checks (run through balanceTest, 2026-06-17)

Script: vignettes/extreme-cases-sim.R (devtools::load_all, matched pairs,
strata(pair)). Two extremes confirm that the within-set omnibus d^2 is blind to
imbalance MAGNITUDE by construction: d^2 is homogeneous of degree zero in the
within-pair differences (scale every c_s by t and dbar scales by t, V_d by t^2,
so the ratio dbar' V_d^{-1} dbar is unchanged).

### Case 1: 99 identical pairs + 1 pair off by eps on one covariate

Predicted: V_d is rank 1, so df = 1, and d^2 = (eps/100)^2 / (eps^2/100^2) = 1
EXACTLY for any eps > 0; large-sample p = P(chisq_1 >= 1) = 0.317; permutation p
= 1 (the single flippable pair gives a quadratic form invariant to the flip --- a
point mass at 1). Observed (pair-stratified omnibus):

    eps=1e+00 : chisq=1  df=1  p=0.317311
    eps=1e-01 : chisq=1  df=1  p=0.317311
    eps=1e-02 : chisq=1  df=1  p=0.317311
    eps=1e-04 : chisq=1  df=1  p=0.317311
    eps=1e-06 : ERROR: Cannot calculate pseudoinverse: ... all covariates
                constant (within strata)?
    eps<=1e-06: ERROR (same)

The nonlinearity is a flat line plus a cliff, NOT a ramp. d^2 does not slide to 0
as delta shrinks; it is pinned at the rank (= 1 here) for every detectable eps,
then falls off a cliff once the lone singular value eps^2/1e4 drops below the SVD
tolerance (between eps=1e-4 and 1e-6) and V_d can no longer be inverted. The
pseudoinverse error IS the rank-degeneracy discontinuity; at exactly eps=0, V_d
is the zero matrix and d^2 is 0/0.

Same-printout trilemma. For the eps=1 data, $overall had two rows:

    pair  chisq=1.000   df=1  p=0.317     <- block-respecting (within-pair)
    --    chisq=0.0046  df=3  p=0.9999    <- block-ignoring (unstratified)

The block-ignoring "--" row scales with eps^2 (one perturbed unit among 200
barely moves the marginal mean -> p ~ 1, and grows with eps); the block-
respecting "pair" row is frozen at d^2 = 1. That is the (2a)-tracks-magnitude vs
(2b)-collapses dichotomy printed on a single object.

### Case 2: anti-match (100 pairs), who-gets-treated is the crucial variable

Univariate, within-pair gap = 5:

    adversarial (always treat the high unit): chisq=100  df=1  p=1.52e-23
    random within-pair (200 reps):  mean chisq=1.04 (theory ~1), mean p=0.51,
                                    p approximately uniform
    multivariate random anti-match (k=3): mean chisq=3.02 (theory ~ k=3)

- Random assignment within the awful pairs: signs cancel, dbar ~ 0, d^2 ~
  chisq_k, p ~ uniform --- IDENTICAL to a great match. Pairing dissimilar units
  does not create bias under within-pair randomness; it inflates V_d, which d^2
  divides out. Match dissimilarity != imbalance.
- Adversarial assignment (signs aligned): d^2 = S = 100, p ~ 1e-23. Detected
  hard --- but d^2 = S^2 abar^2 / sum(g^2) = S regardless of the gap size. What is
  detected is that all 100 within-pair coin flips landed the same way (prob
  2^-100), i.e. the SIGN structure, not the magnitude.

### Lesson

Across both extremes the within-set d^2 and its p-value read the direction/sign
of imbalance relative to within-pair reshuffling and are exactly blind (degree-
zero homogeneous) to its magnitude. That is why neither extreme behaves like a
balance-quality score, and it is the cleanest possible statement of the collapse:
the perfect match is pinned at d^2 = 1 (then errors at machine precision), the
random anti-match looks like a perfect match, and only sign-aligned assignment
moves the statistic.

## Possible next exhibit

A small simulation filling in the 2x2 directly: build a fine and a coarse
matching on the same data and show (2a) tracks delta but has badly wrong size
under the within-block null, while (2b) has correct size but a p-value flat in
delta. Clean exhibit for the trilemma write-up.
