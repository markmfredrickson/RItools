# Comparing a matched design to a randomized standard: what the omnibus tells you, and what it hides

Date: 2026-06-18. This is the answer we reached after testing six ways to make
the comparison we want. The six candidates and the adversarial tests are in the
workflow record; the numbers below come from one simulated data set, checked
directly against balanceTest() -- read them for the pattern, not the digits.

## The question

You stratify an observational study so that treated and control units sit
together on observed covariates, and you drop units that have no good match.
Then you want to know one thing: is the design I built comparable to a randomized
experiment? You match carefully, the standardized differences look small, and
balanceTest() returns a tiny omnibus p-value that calls the design imbalanced.
Which number do you believe, and what does either one license you to claim?

## The central result: the same-strata randomized comparison cannot track the size of imbalance

The omnibus d^2 compares the imbalance you see to the imbalance you would see if
treatment were randomized within the strata you built. To do that it divides the
observed imbalance by the spread of imbalance those within-strata reshuffles
produce. As the strata get tighter -- as they approach exact matching, where the
covariates are constant within each set -- that spread shrinks toward zero. A
fixed imbalance, even a substantively tiny one, divided by an ever-smaller
number, looks ever more extreme. Multiply every within-stratum imbalance by any
constant and d^2 does not change: the omnibus is homogeneous of degree zero in
the imbalance.

Hold a covariate's treated-control gap fixed at half a standard deviation and
tighten the match (shrink the within-pair noise):

    within-pair noise   omnibus chi-sq   omnibus p    SMD    M (fixed scale)   P = chisq/M
    1.00                36.8             1.3e-09      0.43   0.19              195
    0.50                108.6            2.0e-25      0.49   0.24              455
    0.20                170.9            4.8e-39      0.51   0.22              768
    0.10                192.6            8.6e-44      0.50   0.22              893
    0.05                198.2            5.3e-45      0.52   0.22              921

The gap, read as a standardized difference, stays near 0.5 the whole way. The
omnibus p falls from one-in-a-billion to 5e-45. Nothing about the imbalance
changed. The match got tighter, and the p-value chased the shrinking reference.

We tried six ways to keep the comparison to a randomized standard while making it
track the size of imbalance: a coarser grouping for the randomized standard; an
inflation factor placing the design between exact matching and randomization; a
0.25 SD caliper built into the comparison; a percentile against the analyst's own
design search; a fixed-scale randomization check; and plain interpretation. An
independent agent tried to break each one against five failure modes, with its
own simulations. Five of the six collapse, reduce to a bare magnitude, or quietly
ignore the strata you built. Only the sixth survives, and only as interpretation,
not as a new statistic. This is not a gap in our cleverness. It is a property of
the question: you cannot compare a design to a block-randomized experiment on
these same sets and have the comparison track the size of imbalance.

## What you can do: read the omnibus as size times precision

The omnibus chi-square factors -- exactly with one covariate, approximately with
several -- into two numbers that answer two different questions:

    omnibus chi-square  =  SIZE  x  PRECISION

SIZE is the imbalance measured on a FIXED scale that does not shrink as you
match: the standardized differences balanceTest() already prints, and their
multivariate companion M = dbar' Sigma^{-1} dbar, with Sigma the pooled covariate
covariance fixed BEFORE matching. Exact matching makes M = 0, but not the reverse
(within-stratum mean differences can cancel across strata): M is a necessary, not
sufficient, sign of exact matching. SIZE tells you how far
the design is from exact matching on observed covariates, and it does not
collapse (in the table it sits near 0.2 throughout).

PRECISION is how much the matching tightened the randomization reference:
P = Sigma_pool / V_d, the fixed pooled variance over the within-strata
randomization variance of the mean difference. P grows without bound as the match
tightens -- in the table, from 195 to 921 -- because matching moves covariate
variance from within strata to between strata. P grows with the effective
number of strata times the share of variance moved between strata. It has NO
closed form as simple as 1/(1 - eta^2); that expression omits the sample-size
factor and is wrong as a formula for P.)

So a small omnibus p-value is the PRODUCT of a stable size and an exploding
precision. When the size is small and the p is small, the p is small because of
the precision, not the imbalance.

## How to read the output: three numbers, in one order

1. Read SIZE first -- the standardized differences and M, on the fixed pre-match
   scale. This is the number that says how far you are from exact matching on
   observed X. The p-value does not say this.
2. Read the omnibus p as DIRECTION, knowing the precision P. It answers a real and
   separate question: is the residual a persistent, same-sign gap that
   re-randomizing within your strata would rarely produce?
   - small p, small SIZE, large P: a tiny gap seen through a tight match. The
     match worked; report the size. Do not read the small p as failure.
   - small p, large SIZE: a real, persistent gap. Here the small p is substantive.
   - large p, small SIZE: balanced, and not even a persistent direction.
3. Close with scope, every time. All of this concerns OBSERVED covariates.
   Comparability to a randomized standard on observed X is not ignorability.
   Unobserved confounding is a separate question, for sensitivity analysis
   (Rosenbaum's Gamma).

The rule of thumb for students: the p-value tells you the DIRECTION of the
residual gap relative to chance at your strata; the standardized difference tells
you its SIZE. A tight match makes the p small on purpose. When the p is small and
the standardized differences are small, read the standardized differences.

## Where the 0.25 SD caliper fits, and where it does not

Use 0.25 SD, if at all, as a tolerance on the SIZE -- "is each standardized
difference within 0.25?" -- never as a gate on the p-value, and never as a claim
that being within 0.25 SD justifies ignoring confounding. It is a yardstick for
reading a magnitude, not a balance verdict, and it has no substantive provenance
that would make it one. (It is a matching-literature convention from Rubin 2001
via Stuart 2010, not the Rosenbaum-Rubin 1985 propensity-score caliper it is
usually credited to, and it mixes conventions with RItools' pooled-SD
denominator, whose medical lineage uses 0.1; see balance-threshold-provenance.md.)
Every attempt to build the caliper INTO the randomized comparison failed: the
collapse reappears in the covariate direction the caliper does not control.

## A correction about Tukey: do not teach the platinum/gold framing as his

Tukey's "platinum standard" ranks ANALYSES by how few unverifiable assumptions
the inference needs. Platinum means the p-value rests only on the randomization
the experimenter actually performed; gold, silver, and baser metals lean
progressively more on assumed distribution shapes. His ladder does not rank
DESIGNS by how well they remove confounding. Two consequences:

- "Platinum = exact matching, gold = randomized experiment" is our internal
  shorthand, not Tukey's. Do not attribute it to him.
- By Tukey's own definition, an observational matched design cannot reach the
  platinum standard, because it has no verifiable assignment mechanism -- its
  inference rests on the unverifiable assumption that treatment is as-if-random
  within strata. His distinction puts randomization ABOVE non-randomized
  matching, not below it.

What does survive translation, and is worth citing him for: exact balancing is a
stricter balance TARGET than randomization, because matching forces zero
imbalance while randomization only balances on average -- Tukey calls relying on
randomization's average balance "an inadequate scientific method" (p. 270). And
his motto, "balance what you can, randomize the rest, analyze by rerandomization"
(p. 271), is the spirit of the omnibus. Cite him for those two points, and
disclaim the metal-tier reading for designs.

## For the package

Report the three numbers together: the standardized differences and M (size,
distance from exact matching, on a fixed pre-match scale); the omnibus p
(direction, the genuine within-strata randomization comparison); and P
(precision, the reason a tight match makes the p small). Label M and P as
descriptions, not as a new calibration -- read as a calibration, M is just the
magnitude and P just carries the collapse. The machine-tolerance fixes stand:
drop covariates with no within-stratum variation, with a message rather than an
error; and offer the Cauchy combination when there are more covariates than the
within-strata degrees of freedom can support. For teaching, one plot earns its
place: M against both the within-strata reference (whose percentile is flat across
imbalance size -- the collapse) and a complete-randomization reference (which
tracks size but ignores your strata). The gap between the two curves is the
impossibility, drawn.

## What we did not solve

There is no non-collapsing comparison to a block-randomized experiment on these
same sets. That is settled, not open. What remains open: whether to print M and P
or keep them in teaching material; how to report the precision across several
covariate directions rather than as one number; and how to define "pre-match
pool" when units were dropped. Unobserved confounding stays out of scope here. It
belongs to the sensitivity analysis, which is where the real threat to a causal
claim lives once the observed covariates are handled.
