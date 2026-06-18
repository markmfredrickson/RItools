# Judging a matched design against block-randomized experiments

## A plain writeup of the resolution-profile idea and what the simulation showed

This memo proposes one way to say whether a matched design is good enough, and
reports a simulation that tests the proposal. I have tried to write it plainly.
Where a claim is proved, shown only by simulation, or argued heuristically, I say
so.

---

## 1. What we are trying to do

We have a matched design. Units are grouped into sets. Each set holds treated and
control units (or, in a dose study, higher- and lower-dose units) that are
similar on the measured covariates. We want to argue that this design supports
the assumption behind the analysis: that within sets, treatment is as good as
randomly assigned. The evidence we have is the set of covariate differences
between treated and control units, summed across sets. We want a single summary
of those differences across all covariates at once, and we want to read that
summary against a standard --- what a randomized experiment would produce.

A randomized experiment is the right standard because it is the design we trust.
The question is which randomized experiment. Complete randomization of the whole
sample is one candidate, but it is the wrong one here: matching reorganizes the
sample into similar sets rather than discarding units, so the experiment it
imitates is one that randomizes treatment within blocks, not one that randomizes
the whole sample at once.

---

## 2. The problem we have to solve

The d^2 test compares the covariate imbalance we observe to the imbalance we
would see if treatment were re-randomized within the realized matched sets. As
the match gets tighter, this comparison turns against the analyst. A small, real,
persistent difference between treated and control units gets driven toward a
p-value of zero, so a design that matched better looks worse.

Your team proved this cannot be fixed by changing the statistic. Any statistic
computed from the matched sample, compared to a reference distribution built from
that same sample, returns the same verdict whether the within-set differences are
large or tiny. The proof is in the team's memos; I take it as given.

Here is the mechanism in one line of algebra. Write d for the vector of
treated-minus-control covariate differences summed across sets, and V for the
covariance of d under re-randomization within the matched sets. The d^2 statistic
is

    T = d' V^+ d .

Now imagine shrinking every within-set distance by a factor delta (a tighter
match). The differences d shrink in proportion to delta. The covariance V is
built from within-set sums of squares, so it shrinks in proportion to delta^2.
The ratio T = d' V^+ d then has delta in the numerator and delta^2 in the
denominator, twice over, so the two cancel and T does not change at all. The test
literally cannot tell a match with delta = 1 from a match with delta = 0.001.
That is the collapse.

---

## 3. The idea
<!-- TODO: I'm not sure about the "Strict" versus "loose" distinction. What does this mean and why does itmatter here?-->
Re-randomizing within the realized matched sets is one choice of reference, and it
is the strictest one. Complete randomization of the whole sample is the loosest.
Between these two extremes is a range of references: re-randomize treatment within
blocks that are coarser than the matched sets but finer than the whole sample.

The proposal is to compute the design's percentile against each reference in that
range, from the matched sets up to the whole sample, and report the result as a
curve. I call it the resolution profile. At one end (the matched sets) the curve
reproduces the d^2 test, which collapses. At the other end (the whole sample) it
reproduces the complete-randomization comparison. The interesting part is in
between.

The single number to read off the curve is the coarsest-to-finest crossing point.
Define

    l* = the finest reference blocking at which the design's imbalance is no worse
         than the median block-randomized experiment at that blocking.

In words: l* is the finest block-randomized experiment the design balances
covariates as well as. A design that matches the balance of a finely blocked
experiment has a small l* and is well supported on the measured covariates. A
design that only matches the balance of a coarsely blocked experiment has a large
l*. The number is on a scale a researcher understands: the size of the blocks in
the experiment the design imitates.

---

## 4. Why coarsening the reference rescues the comparison

The collapse in section 2 came from the within-set covariance V shrinking as fast
as the differences d. A coarser reference breaks that link.

When the reference blocks are coarser than the matched sets, each block contains
units from several different matched sets, so the covariate spread within a block
includes the differences between set centers. Those between-set differences are
fixed by where the sets sit; they do not shrink when you tighten the match. So the
coarse reference's covariance stays bounded away from zero as delta goes to zero,
while the imbalance the coarse reference measures shrinks toward whatever residual
the matching could not remove. The numerator shrinks; the denominator does not;
the ratio falls. The design looks better as the match improves, which is what we
want.

Stated as the lesson: a matched design's value is the between-set covariate
balance it achieved, and the within-set reference removes exactly that balance
from view before measuring. To see what the matching accomplished, compare against
a randomization that keeps the between-set variation in view, which is any
blocking coarser than the matched sets.

Status of this claim: the collapse at the finest reference is proved (by the
team). The immunity of coarser references is an algebraic scaling argument plus
the simulation in section 6. I have not proved it in general.

---

## 5. The construction

Notation follows balanceTestEngine and sigma_x_test.

- Units i = 1..N, covariates X (N by K), treatment contrast c.
  For a bipartite design c_i = z_i - pi_b, the treated indicator centered within
  its block b. For a dose design c_i = D_i - Dbar_b, the dose centered within b.
- For a blocking B, d_B = sum_i c_i x_i is the within-block-adjusted difference.
- V_B = sum over blocks b of [ (sum_{i in b} c_i^2) / (n_b - 1) ] S_{xb}, with
  S_{xb} the within-block covariate scatter, is the covariance of d_B under
  re-randomization within the blocks of B. This is randomization_cov_d.
- T_B = d_B' V_B^+ d_B is the d^2 statistic at blocking B.
- P_B = the fraction of within-block-randomized assignments whose T_B is at least
  the observed T_B. This is the percentile: small means the design's imbalance is
  large for that reference, large means the design beats that reference.

The ladder of references runs from the matched sets (finest) to one block
(coarsest). The continuous version indexes the ladder by a caliper h on a 1-D
matching score (a propensity score, say): block together matched sets whose score
locations fall within h of each other. At h = 0 the blocks are the matched sets;
at h larger than the score range the block is the whole sample. The profile is
P(h), and l* is the block size at which P(h) first reaches 1/2.

One rule must be fixed before the profile is read, and it is forced by the idea of
the ladder itself. A block in a real block-randomized experiment holds units that
resemble each other; it blocks on a covariate. So the only honest way to coarsen is
to merge the matched sets that are nearest on the matching score the design was
built from, and that rule must be declared in advance. An analyst who instead merges
sets that are far apart on the score builds a reference whose blocks are internally
heterogeneous, which lets the design look as-good-as-random at a finer resolution
than it earns. Section 6 shows this is not a small effect: the choice of coarsening
rule moves l* by a full rung, and always in the direction that flatters the design.
The script vignettes/resolution-profile-threebenchmark.R implements the
geometry-respecting ladder (merge nearest sets on the score) and the interleaved
ladder that violates it, so the gap is reproducible.

---

## 6. What the simulation showed

Two scripts: resolution-profile-demo.R (the collapse and the discrete ladder) and
resolution-profile-stability.R (the continuous caliper and the stability study).
The data: matched sets spread along a latent score, a fixed per-covariate
treated-control gap g (the residual the match leaves), and within-set noise scaled
by a parameter that shrinks as the match tightens. Percentiles use the analytic
chi-square form, checked against permutation.

The collapse, and that it sits only at the finest reference. As the match
tightens, the percentile against the matched-set reference falls from 0.47 to
0.0002, while the percentile against complete randomization stays near 0.99. Same
design, opposite verdicts. Stepping one rung coarser than the matched sets removes
the collapse: that percentile holds near 0.74 no matter how tight the match. The
collapse is confined to the single finest reference (simulation, plus the scaling
argument of section 4).

l* tracks design quality (validity). Median l*, measured as the mean reference-
block size at the crossing, rises monotonically with the persistent gap:

    gap   0.0   0.1   0.2   0.4   0.8
    l*    3.1   4.9   8.1   14.5  21.9   (units per reference block)

A gap of zero means the design already is a within-set random assignment; its
profile sits near 1/2 at every resolution, and l* is correctly ill-defined (no
crossing to find). A large gap pushes l* toward the whole sample, meaning only a
coarsely blocked experiment would balance as poorly.

l* holds up against the choices I get to make (robustness). The bin offset of the
caliper grid moves the raw percentile at a fixed h by about 0.25; averaging over
offsets removes that. Building the ladder from a covariate instead of the true
score gives essentially the same l* (within-design median difference 0.21 block-
size units). The one real weakness: a heavily noisy score shrinks l* (from 11.4
to 8.9), because a bad score scrambles the blocks and makes the design look
as-good-as-random at a finer resolution than it earns. Build the ladder from the
score the matching used.

The coarsening rule is the larger lever, and it deserves to be stated plainly
rather than buried with the bin offset. On one fixed design, merging the nearest
matched sets on the score (the geometry-respecting rule) and merging far-apart sets
(interleaved blocks) give l* of about 11.5 and 6.4 units per block --- a full rung
apart, larger than the bin-offset wiggle, and always in the direction that flatters
the design. The interleaved ladder is not a strawman: any caliper or clustering that
ignores the score geometry drifts toward it. The discipline that closes this degree
of freedom is the one in section 5: pre-specify the ladder as nearest-set merges on
the matching score. (Reproduced in resolution-profile-threebenchmark.R; this is a
real defect the earlier draft understated, and the fix is cheap.)

A second result cuts the other way and answers the standard objection to a
design-respecting benchmark --- that balance is about the covariate distributions in
the treated and control groups, not about who is paired with whom. l* is invariant
to the pairing. Holding the same units and the same treated and control covariate
distributions, three different valid pairings of the same 48 units give resolution
profiles that agree to about 0.001. The reason is structural: at every reference
coarser than the matched sets, the d^2 statistic depends on the matched-set means
and on which units sit in each resolution band, not on who is paired with whom inside
a set. So l*, which lives above the finest rung, reads off the treated and control
distributions sorted into resolution bands --- exactly the object the objection says
balance should be about. Only the single finest rung uses the pairing, and that is
the rung we already know collapses.

The analytic backend is accurate. The chi-square percentile matches a 4000-draw
within-block permutation (0.018 vs 0.012 at the finest checked rung, 0.534 vs
0.552 at the coarsest). The rung where chi-square is least accurate is the finest,
which sits below l* and does not drive the crossing. Use simulation or enumeration
at the finest rung, as the existing R < 1e4 gate already does.

---

## 7. How to think about uncertainty: randomization, not resampling

The inference model is the assignment mechanism. Conditional on the units, their
covariates, and the blocks, the only randomness is which units are treated. The
percentile P(h) is therefore an exact, finite-sample quantity: it is the position
of the observed design in the distribution of d^2 over within-block re-randomized
assignments at resolution h. We compute it exactly by enumeration when the number
of assignments is small, or estimate it by drawing assignments, in which case the
only error is Monte Carlo error in the percentile, which we control by the number
of draws. There is no superpopulation and no resampling of units, so there is no
bootstrap.

For a single design, then, the curve P(h) is the complete answer and l* is one
reading of it. Report the curve.

The spread I reported in the stability study --- l* varying by 30 to 40 percent
across replications --- is a different thing. It is the variation of l* across many
designs generated at the same nominal quality. It tells me that l* is a somewhat
noisy summary of design quality, so two designs with the same true residual can
give l* values that differ by a third. That is worth knowing when comparing
designs, but it is not an uncertainty band for one analyst's one design.

---

## 8. What is established, and what is not

Proved (by the team): the d^2 collapse at the matched-set reference, and that no
statistic from the matched sample compared to a same-sample reference escapes it.

Shown by simulation and supported by a scaling argument, not proved: that coarser
references do not collapse; that l* rises monotonically with the residual gap;
that l* is stable to the bin offset and to a reasonable score; that l* is invariant
to the within-set pairing (three pairings of the same units with the same marginals
agree to about 0.001); and that the geometry-respecting ladder (merge nearest sets
on the score) removes the coarsening-rule sensitivity that an interleaved ladder
otherwise opens.

Observed, not proved: that P(h) is near-monotone in h once there is a residual to
detect (Spearman correlation about 0.99). I used isotonic regression to read off a
single crossing, justified by this observation, not by a theorem. The honest
object is the curve.

Open: a formal statement and proof that coarser references escape the collapse; a
formal version of the pairing-invariance result and of why the geometry-respecting
ladder is the right pre-specification (both shown only by simulation here); the unit
for l* (units per block is not comparable across set structures, so
matched-sets-per-block or h/sd(score) may be better for comparing designs); and a
prototype of the continuous caliper, with the nearest-set-merge ladder, in
balanceTest().

---

## 9. What the method does not do

It speaks to the measured covariates only. A small l* says the design balances the
measured covariates as well as a finely blocked experiment would; it says nothing
about unmeasured confounders, and "as well as experiment l*" is a statement of
consistency, not a certificate of balance. Report it next to a sensitivity analysis,
not in place of one.
