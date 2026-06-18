# E-values for balance assessment: a short memo

Author: Jake Bowers (with Claude). Companion to `sigma-x-rosenbaum-memo`.
Draft for Ben, Mark, and Josh. Status: conceptual proposal plus a small
numerical check. No code in the package yet.

## Why this memo

The sigma_x memo ends at an impossibility: no balance measure computed
from the matched sample alone can reward absolute balance, because the
statistic and its within-sample reference distribution rescale together.
That argument was about p-values, and the natural next question is
whether e-values escape it. They do not escape the scale problem. But
they answer a different question we care about and have not solved: how
to give a user one trustworthy omnibus number after they have looked at
many covariates, and after they have searched over many candidate
matches. That is a multiplicity problem, and e-values are built for it.

This memo separates the two issues, states what e-values buy, and is
honest about what they cost.

## What an e-value is, in one paragraph

An e-value is a nonnegative statistic `E` whose expected value under the
null is at most 1: `E_{H0}[E] <= 1`. Read it as the payoff of a bet
against the null that breaks even under the null; a large realized payoff
is evidence against the null. By Markov's inequality,
`P_{H0}(E >= 1/alpha) <= alpha`, so rejecting when `E >= 1/alpha` is a
level-alpha test for any sample size --- no asymptotics. The price for
that finite-sample validity is conservativeness: the Markov bound is
loose, so an e-value test has less power than a well-calibrated exact
test when the model is right. Everything below trades power for
robustness and for protection against multiple looking.

## What e-values do not fix: the scale problem

For completeness, because it is the first thing a reader will hope for.
Any e-value built from the matched sample alone and made "equivariant"
--- meaning a within-stratum rescaling of the covariates multiplies `E`
at every assignment by the same factor --- is forced to be exactly
scale-invariant. The e-value validity constraint must hold at every
scale at once, and a single rescaling factor cannot satisfy it unless
that factor is 1. So such an e-value sees relative imbalance only, exactly
like HB08. The only escape is a bet that carries an absolute scale (a
fixed betting parameter, or a prior, in units of the covariate), which is
the sigma_x memo's "external yardstick" wearing different clothes. The
pool-calibrated bet in the last section does see absolute imbalance,
precisely because its scale comes from the pre-matching pool. E-values
relocate the scale question; they do not remove it.

## What e-values do fix: multiplicity, two kinds

This is the reason to care, and it maps onto the "omnibus measure to
avoid multiple testing" goal.

### Across covariates: combine by averaging, under any dependence

A user testing balance on `p` covariates faces `p` tests. The HB08
chi-square already combines them into one omnibus, but it does so by
inverting `Cov(d)` (or `Sigma_x`), which is the part of the machine that
struggles with collinear covariates and small `n`, and it leans on a
chi-square approximation that is imperfect at small `n` (the 0.044 gap in
the handoff notes).

E-values combine more cheaply and with a stronger guarantee. If
`E_1, ..., E_p` are e-values for the same null, then their arithmetic
mean `(1/p) sum_j E_j` is again an e-value:

    E_{H0}[(1/p) sum_j E_j] = (1/p) sum_j E_{H0}[E_j] <= 1,

by linearity of expectation alone. No independence assumption. The
covariates can be arbitrarily correlated and the averaged omnibus is
still valid (Vovk and Wang 2021). This matters here because balance
covariates are correlated by construction, and the chi-square omnibus
needs their covariance to combine them, while the e-value omnibus does
not.

The check on the sigma_x toy (two correlated covariates, within-stratum
randomization null) confirms the averaged omnibus has null mean exactly 1:

```
                    delta=1 delta=0.5 delta=0.1 delta=0.01
omni_linear        2.000000  2.000000  2.000000   2.000000   # within-sample, scale-invariant
omni_linear_nullmu 1.000000  1.000000  1.000000   1.000000   # valid e-value
omni_pool_exp      1.086751  1.021244  1.000844   1.000008   # pool-calibrated, sees delta
omni_pool_nullmu   1.000000  1.000000  1.000000   1.000000   # valid e-value
```

The within-sample omnibus (`omni_linear`) is flat in delta, as the
impossibility result requires. The pool-calibrated omnibus
(`omni_pool_exp`) is a single number that both (a) combines the two
covariates without a multiplicity correction and without inverting a
covariance, and (b) moves toward 1 (no evidence of imbalance) as the
match tightens, because its scale is the pool SD. Both are valid e-values
under the dependence between `x1` and `x2`.

For per-covariate flagging rather than a single omnibus, the e-BH
procedure (Wang and Ramdas 2022) controls the false discovery rate under
arbitrary dependence, which is a cleaner guarantee than the Holm
adjustment `balanceTest` now applies to the per-covariate z-scores.

### Across specifications: multiply, and stop whenever you like

The multiplicity that the sigma_x memo raises but leaves open is the
caliper search. A matcher tries several calipers, several distances,
several rules for dropping units, checks balance each time, and reports
the match they keep. Each choice changes the surviving sample and the
reference distribution, so the reported p-value is a minimum over a
search and is not valid.

E-values are the standard tool for this. If each analysis in a sequence
produces an e-value that is valid conditional on everything before it,
the running product `M_T = prod_{t<=T} E_t` is a nonnegative
supermartingale starting at 1, so by Ville's inequality
`P_{H0}(sup_t M_t >= 1/alpha) <= alpha`. The user may peek after every
candidate match and stop whenever they like; rejecting when the product
crosses `1/alpha` still controls the type I error over the whole search
(Shafer 2021; Ramdas, Grunwald, Vovk, and Shafer 2023). This is the
anytime-valid guarantee, and it is exactly the protection a matcher who
iterates on the design needs and currently does not have.

One technical obligation makes this honest rather than automatic. As the
caliper changes, the matched set changes, so the within-set
randomization null changes, and "valid conditional on the past" has to be
arranged against a moving target. The clean way to fix this is to hold
the null fixed --- the complete-randomization null on the pre-matching
pool --- and let each specification enter through the bet (which units
and covariates you wager on from the current match), not through the
null. The matched set then shows up in the numerator, the pool in the
denominator, which is the same split the sigma_x memo arrived at for the
scale problem. I believe this construction works; I have not proved the
conditional-validity step, and I am flagging it as the open piece.

## What I am proposing, concretely

1. Offer a per-covariate e-value and an averaged omnibus e-value as an
   optional output of `balanceTest`, alongside the existing chi-square.
   The averaged omnibus is one number, valid under covariate dependence,
   finite-sample valid, and needs no matrix inverse. Default scale:
   within-sample (scale-invariant, tests the randomization null like
   HB08). When a `candidate_pool` is supplied, switch the bet to the
   pool scale, so the same omnibus also rewards absolute balance.

2. Offer e-BH per-covariate flagging under arbitrary dependence as an
   alternative to Holm.

3. Treat the sequential caliper-search e-value as a research direction,
   not a v1 feature. It needs the fixed-pool-null construction worked out
   and the conditional-validity step proved. If it holds, it is the most
   useful thing here: a single number a matcher can trust no matter how
   long they searched.

## Costs, stated plainly

E-value tests are conservative; expect less power than HB08 when the
chi-square model is accurate and the user did exactly one analysis. The
case for them is not power. It is (a) validity that does not depend on
the chi-square approximation at small `n`, (b) an omnibus across
correlated covariates that needs no covariance inverse, and (c) honest
inference for a user who looked many times. Where none of those apply,
HB08 is the better choice, and the documentation should say so.

## Proof-status summary

- Averaged e-value is a valid omnibus under arbitrary dependence:
  proved (one line above), standard (Vovk and Wang 2021), and confirmed
  numerically on the toy.
- Equivariant within-sample e-value is forced scale-invariant: proved
  in the prior analysis; consistent with `omni_linear` flat in delta.
- Pool-calibrated bet sees absolute imbalance: follows from the
  external-scale characterization; confirmed by `omni_pool_exp`.
- Sequential caliper-search product is anytime-valid: conditional on a
  conditional-validity step I have not yet proved for the moving-null
  setting. Direction, not result.

## References (verify metadata before circulating)

- Vovk, V. and Wang, R. (2021). E-values: calibration, combination and
  applications. Annals of Statistics 49(3), 1736--1754.
- Wang, R. and Ramdas, A. (2022). False discovery rate control with
  e-values. Journal of the Royal Statistical Society Series B, 84(3).
- Shafer, G. (2021). Testing by betting: a strategy for statistical and
  scientific communication. Journal of the Royal Statistical Society
  Series A, 184(2), with discussion.
- Ramdas, A., Grunwald, P., Vovk, V., and Shafer, G. (2023).
  Game-theoretic statistics and safe anytime-valid inference. Statistical
  Science 38(4).
- Grunwald, P., de Heide, R., and Koolen, W. (2024). Safe testing.
  Journal of the Royal Statistical Society Series B, 86, with discussion.
