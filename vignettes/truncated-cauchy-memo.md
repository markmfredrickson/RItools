# Truncating the Cauchy combination in balanceTest()

On 2026-10-25 Jake chose to have `acat_pvalue()` use riposte's truncated
Cauchy combination with truncation share t = 0.9. `acat_pvalue()` computes
`cauchy_comb_p` when `balanceTest()` is called with
`cauchy.combination = TRUE`, and commit 3f80d9d makes the change. This memo
records the evidence he had when he chose.

We had first agreed to make the switch and to write the tests first. Before
writing them I simulated how often each version gives a p-value below 0.05
when the null hypothesis is true. The truncated version does so more often.
We had not weighed that, so Jake chose again from the three options at the end
of this memo, with the results below in hand.

## The problem the truncation fixes

The Cauchy combination of Liu and Xie (2020) converts each per-covariate
p-value p to tan((0.5 - p) pi) and averages the converted values. Its omnibus
p-value is the probability that a standard Cauchy variable exceeds that
average. A p-value of 1 converts to minus infinity. The branch capped p at
1 - 2.2e-16, and the cap converts to about -1/(2.2e-16 pi) = -1.4e15. So with
one covariate at p = 1, the omnibus p-value was still about 1.

In `balanceTest()` a covariate's p-value is exactly 1 when its adjusted
difference is exactly 0. We built a 40-unit, 4-block example with three
covariates. The first is binary, with the same treated and control share in
every block, so its p-value is 1. The second is imbalanced, with p-value
1.5e-05, and the third is noise, with p-value 0.63. Before the change, the
branch's `cauchy_comb_p` was 0.99999999999999933. Without the binary
covariate it was 3.1e-05.

riposte 0.0.0.9008 uses the truncated conversion of Gui, Jiang, and Wang
(2025, Biometrika 112(4), asaf038) instead. Each p becomes
tan((0.5 - t p) pi), with t = 0.9 by default, so p = 1 becomes
tan(-0.4 pi) = -3.08. If the n converted values sum to S, the omnibus p-value
is min(1, n P(C > S) / t), where C is a standard Cauchy variable. In the
example this gives 4.6e-05.

## Rejection rates under the null hypothesis

When the null hypothesis is true, a test at the 0.05 level should give a
p-value below 0.05 in 5 percent of samples. We drew 40,000 sets of K z-scores.
Each set came from a normal distribution with mean 0 and correlation rho
between every pair of z-scores. We turned each z-score into a two-sided
p-value and combined each set. The standard error of a rejection rate near
0.05 is sqrt(0.05 x 0.95 / 40000) = 0.0011, and near 0.01 it is 0.0005. In
the table, t = 1 is the branch's untruncated combination before the change,
and t = 0.9 is riposte's default.

| K  | rho | t = 1 at 0.05 | t = 0.99 at 0.05 | t = 0.9 at 0.05 | t = 1 at 0.01 | t = 0.9 at 0.01 |
|----|-----|------|------|------|-------|-------|
| 6  | 0   | 0.0513 | 0.0529 | 0.0573 | 0.0103 | 0.0106 |
| 6  | 0.5 | 0.0583 | 0.0585 | 0.0604 | 0.0118 | 0.0119 |
| 40 | 0   | 0.0497 | 0.0550 | 0.0626 | 0.0100 | 0.0108 |
| 40 | 0.5 | 0.0647 | 0.0664 | 0.0691 | 0.0128 | 0.0129 |

1. At the 0.05 level the rejection rate with t = 0.9 is higher than with
   t = 1 in every row. With 40 independent covariates it goes from 0.0497 to
   0.0626. That increase, 0.0129, is about 12 standard errors and the largest
   in the table.
2. When the z-scores are correlated, the rejection rate at 0.05 is already
   above 0.05 with t = 1: 0.058 with 6 covariates and 0.065 with 40. With
   t = 0.9 it is 0.2 to 0.4 percentage points higher than that.
3. At the 0.01 level every version's rejection rate is within about 0.003 of
   0.01.

The value of t also matters when one covariate has a moderately small p-value
and another has p = 1. For two covariates with p-values 0.03 and 1, the
omnibus p-value is 1 with t = 1, 1 with t = 0.99, and 0.081 with t = 0.9. For
five covariates with p-values 0.01, 0.01, 0.01, 0.01, and 1, it is 1, 0.017,
and 0.013.

## Which error matters for balance checking

In a balance test, a researcher who rejects the null hypothesis flags the
design as imbalanced. When the rejection rate is too high, the researcher
flags some well-balanced designs and re-matches or adjusts when it was not
needed. The p = 1 failure goes the other way. When one covariate's adjusted
difference is exactly 0, the omnibus p-value is about 1 even if another
covariate is badly imbalanced, and the researcher goes ahead with an
imbalanced design. In the table, the first error adds at most 1.3 false flags
per 100 balanced designs tested at 0.05. In the example above, the second
error reports an imbalance with p = 1.5e-05 as p = 1. I judged the second
error the worse one for a balance check, and that judgment is the argument for
truncating.

## Tests that changed

Two existing tests in `tests/testthat/test.cauchycomb.R` checked properties
that hold for Liu and Xie's combination and not for the truncated one:

- The idempotence test expected that combining five copies of 0.3 returns
  0.3. The truncated combination returns 0.394. Each copy converts to
  tan((0.5 - 0.27) pi) = 0.882, the five sum to S = 4.41,
  P(C > 4.41) = 0.0710, and 5 x 0.0710 / 0.9 = 0.394. For a single p-value p
  it does return p, because P(C > tan((0.5 - t p) pi)) = t p and
  t p / t = p.
- The calibration test expected the combined p-value to be uniform under the
  null hypothesis, with mean 0.5. With 6 independent uniform p-values, the
  truncated combination has mean 0.61, and 41 percent of its values are
  exactly 1. The formula min(1, n P(C > S) / t) gives 1 whenever
  6 P(C > S) / 0.9 exceeds 1, and that happens when S is below 1.96.

Both tests now run with `truncation = 1`, where they still hold. New tests
check that a single p-value combines to itself, that the rejection rate at
0.01 is 0.01 +/- 0.003, and that the 40-unit example gives a
`cauchy_comb_p` below 0.001. The end-to-end null test in
`tests/testthat/test.notforCRAN.R` expected a mean of 0.5. It now requires a
mean between 0.5 and 0.7, and the mean in that test was 0.615.

## A note for riposte

riposte's test of the truncated combination's level uses 1,000
re-randomizations. It accepts any rejection rate below
0.05 + 3 x 0.0069 = 0.071, where 0.0069 = sqrt(0.05 x 0.95 / 1000) is the
standard error of a rejection rate near 0.05 from 1,000 draws. A rate of 0.06,
about what the table shows with t = 0.9, would pass that test.

So we ran riposte's own design with many more draws. With Wallsten and
Nteta's 36 outcomes held fixed and treatment reassigned 100,000 times, the
truncated Cauchy combination rejected at 0.0487 at the 0.05 level and the
hybrid at 0.0471. Over 50,000 fresh data sets with no effect (36 units, a
5-point outcome), the rates were 0.0507 and 0.0495. riposte's six
representations of one outcome are strongly dependent, and in these designs
the truncated combination kept its level. The script is
`dev/level-truncated-cauchy.R` in riposte. riposte's level test now uses
10,000 reassignments, and its test of the permutation hybrid uses 5,000 data
sets instead of 300.

## Code

```r
source("~/repos/riposte/R/truncated-cauchy.R")
source("~/repos/riposte/R/asymptotic-engine.R")  # riposte_cauchy_log_tails(), used when t = 1
set.seed(12)
B <- 40000
out <- NULL
for (K in c(6, 40)) for (rho in c(0, 0.5)) {
  Z <- sqrt(rho) * rnorm(B) + sqrt(1 - rho) * matrix(rnorm(B * K), B, K)
  P <- 2 * pnorm(-abs(Z))
  for (t in c(1, 0.99, 0.9)) {
    pc <- apply(P, 1, riposte_truncated_cauchy, truncation = t)
    out <- rbind(out, data.frame(K = K, rho = rho, t = t,
                                 r05 = mean(pc < 0.05), r01 = mean(pc < 0.01)))
  }
}
out
```

## Decision

Jake chose the first of these three options on 2026-10-25.

- t = 0.9, which I recommended. It matches riposte and the value Gui, Jiang,
  and Wang recommend. With it, a moderately small p-value next to a p-value of
  1 still gives a small omnibus p-value. It raises the rejection rate at 0.05
  by up to 1.3 percentage points (0.0497 to 0.0626 with 40 independent
  covariates).
- t = 0.99. It raises the rate by at most about 0.5 percentage points. But
  with one p-value of 1, the omnibus p-value is 1 even when another
  covariate's p-value is 0.03.
- Keep Liu and Xie's conversion. Under the null hypothesis, with independent
  covariates, its combined p-value is exactly uniform. It leaves the p = 1
  failure in place.
