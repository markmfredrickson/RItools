# Truncating the Cauchy combination in balanceTest()

The decision this memo asks for: should `acat_pvalue()`, which computes
`cauchy_comb_p` when `balanceTest()` is called with `cauchy.combination = TRUE`,
switch to riposte's truncated Cauchy combination, and with which truncation
share t? I recommend t = 0.9, riposte's default.

On 2026-10-25 we decided to make the switch, writing the tests first. Before
writing them I simulated how often each version rejects a true null. The
truncated version rejects more often at the 0.05 level, which we had not
weighed when we decided, so the decision is open again. The options are at the
end.

## The problem the truncation fixes

The Cauchy combination of Liu and Xie (2020) converts each per-covariate
p-value p to tan((0.5 - p) pi), averages the converted values, and reads the
omnibus p-value off the tail of a standard Cauchy distribution. A p-value of 1
converts to minus infinity. The branch caps p at 1 - 2.2e-16, which converts
to about -1/(2.2e-16 pi) = -1.4e15, so one covariate with p = 1 still sets the
omnibus p-value to about 1.

In `balanceTest()` a covariate's p-value is exactly 1 when its adjusted
difference is exactly 0. In a 40-unit, 4-block example with a binary covariate
that has the same treated and control share in every block (p = 1), an
imbalanced covariate (p = 1.5e-05), and a noise covariate (p = 0.63), the
branch's `cauchy_comb_p` is 0.99999999999999933. Without the binary covariate
it is 3.1e-05.

riposte 0.0.0.9008 uses the truncated conversion of Gui, Jiang, and Wang
(2025, Biometrika 112(4), asaf038) instead: p becomes tan((0.5 - t p) pi),
with t = 0.9 by default, so p = 1 becomes tan(-0.4 pi) = -3.08. If the n
converted values sum to S, the omnibus p-value is min(1, n P(C > S) / t), where
C is a standard Cauchy variable. In the example this gives 4.6e-05.

## How often each version rejects a true null

The truncation removes the p = 1 failure, but it changes how often the
combination rejects when every covariate is balanced in expectation. Under a
true null, a test at the 0.05 level should reject 5 percent of the time. I drew 40,000 sets of K z-scores, each set from a normal distribution
with every pair of z-scores correlated rho, turned each z-score into a
two-sided p-value, and combined each set. The standard error of a rejection
rate near 0.05 is sqrt(0.05 x 0.95 / 40000) = 0.0011, and near 0.01 it is
0.0005. In the table, t = 1 is the branch's current, untruncated
combination, and t = 0.9 is riposte's default.

| K  | rho | t = 1 (branch now) at 0.05 | t = 0.99 at 0.05 | t = 0.9 (riposte) at 0.05 | t = 1 at 0.01 | t = 0.9 at 0.01 |
|----|-----|------|------|------|-------|-------|
| 6  | 0   | 0.0513 | 0.0529 | 0.0573 | 0.0103 | 0.0106 |
| 6  | 0.5 | 0.0583 | 0.0585 | 0.0604 | 0.0118 | 0.0119 |
| 40 | 0   | 0.0497 | 0.0550 | 0.0626 | 0.0100 | 0.0108 |
| 40 | 0.5 | 0.0647 | 0.0664 | 0.0691 | 0.0128 | 0.0129 |

Three things show in the table.

1. Truncation with t = 0.9 raises the rejection rate at 0.05. With 40
   independent covariates it rises from 0.0497 to 0.0626, an increase of
   0.0129, about 12 standard errors. That is the largest increase in the
   table.
2. When the z-scores are correlated, the branch's current combination already
   rejects too often at 0.05: 0.058 with 6 covariates and 0.065 with 40.
   Truncation adds 0.2 to 0.4 percentage points to that.
3. At the 0.01 level every version stays within about 0.003 of 0.01.

The value of t also matters when a moderate signal meets a p-value of 1. For
two covariates with p-values 0.03 and 1, the omnibus p-value is 1 with t = 1,
1 with t = 0.99, and 0.081 with t = 0.9. For five covariates with p-values
0.01, 0.01, 0.01, 0.01, and 1, it is 1, 0.017, and 0.013.

## Which error matters for balance checking

In a balance test, rejecting the null means flagging the design as
imbalanced. A combination that rejects too often flags some well-balanced
designs, and the researcher re-matches or adjusts when it was not needed. The
p = 1 failure goes the other way: one exactly balanced covariate hides real
imbalance on another, and the researcher goes ahead with an imbalanced design.
The table puts the first error at no more than 1.3 extra false flags per 100
balanced designs tested at 0.05. The example above puts the second at an
imbalance with p = 1.5e-05 reported as p = 1. I think the second error is the
worse one for a balance check, which is the argument for truncating.

## Tests that would change

Two existing tests in `tests/testthat/test.cauchycomb.R` encode properties
that hold for Liu and Xie's combination and not for the truncated one:

- The idempotence test expects that combining five copies of 0.3 returns 0.3.
  The truncated combination returns 0.394: each copy converts to
  tan((0.5 - 0.27) pi) = 0.882, the five sum to S = 4.41, P(C > 4.41) = 0.0710,
  and 5 x 0.0710 / 0.9 = 0.394. For a single p-value p it does return p,
  because P(C > tan((0.5 - t p) pi)) = t p and t p / t = p.
- The calibration test expects the combined p-value to be uniform under the
  null, with mean 0.5. With 6 independent uniform p-values, the truncated
  combination has mean 0.61, and 41 percent of its values are exactly 1.
  The formula min(1, n P(C > S) / t) gives 1 whenever 6 P(C > S) / 0.9
  exceeds 1, which happens when S is below 1.96.

Porting would replace them with riposte's single-p-value test and with
rejection-rate checks at 0.05 and 0.01. Neither test would be dropped without
a replacement.

## A note for riposte

riposte's test of the truncated combination's level draws 1,000
re-randomizations and accepts any rejection rate below 0.05 + 3 x 0.0069 =
0.071, where 0.0069 = sqrt(0.05 x 0.95 / 1000) is the standard error of a
rejection rate near 0.05 from 1,000 draws. A rate of 0.06, the size of the inflation above, would pass it.

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

Should `acat_pvalue()` use the truncated conversion, and with which t?

- t = 0.9 (my recommendation). It matches riposte and Gui, Jiang, and Wang's
  recommended value, and it keeps a moderate signal visible next to a p-value
  of 1. It raises the rejection rate at 0.05 by up to 1.3 percentage points
  (0.0497 to 0.0626 with 40 independent covariates).
- t = 0.99. It raises the rate by at most about 0.5 percentage points, but a single
  p-value of 1 still hides a p-value of 0.03.
- Keep Liu and Xie's conversion. It stays exact for independent covariates and
  leaves the p = 1 failure in place.
