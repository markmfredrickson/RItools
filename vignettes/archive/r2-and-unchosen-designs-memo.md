# An R^2 omnibus, and how to write about comparing to designs you did not choose

Jake Bowers, with Claude. 2026-06-16. Branch: devel-sigma-x-omnibus.
Companion to sigma-x-team-decision-memo.md, resolution-profile-memo.md, and
e-values-memo.md. Every quantitative claim below was checked in R by the agents;
where a claim is proved, shown only by simulation, or argued heuristically, I say
so.

This memo explores two ideas you raised. The first: could an R^2-style "how well
do the covariates explain treatment" metric be the omnibus balance summary? The
second: how do we write about l* and the resolution profile, where we compare a
design we chose on purpose against block-randomized designs we did not choose?

The short version. (1) An R^2 of treatment on covariates is the d^2 test on an
interpretable [0,1] scale -- the same number you already have, relabeled -- so it
does the same thing as the p-value and inherits the same collapse. The part of
your instinct that pays off is the model-selection part: a chance floor of
K/(n-S), a high-dimensional regularization story d^2 currently lacks, and flexible
predictors that catch the non-additive imbalance d^2 is blind to. (2) The writing
worry dissolves once we notice you already published the move (Rabb et al. 2022).
l* generalizes one accepted comparison into a graded family. Two new simulation
results sharpen the writeup: l* is invariant to who-is-paired-with-whom (a direct
answer to the standard objection), but it depends on the coarsening rule (a real
defect we must pre-empt by pre-specifying the ladder).

---

## Part 1. An R^2 omnibus

### 1.1 The intuition, and what it gets right

The LR and F tests embody one idea: in a real experiment the covariates should
not predict treatment assignment. An R^2 makes that idea an effect size. If the
covariates explain treatment, the stratification did a poor job; if they do not,
it did well. R^2 lives on [0,1], so it reads as a quantity, not a p-value. That is
the appeal, and it is real: a reader who cannot interpret a chi-square on 7 degrees
of freedom can interpret "the covariates explain 4% of the within-block variation
in treatment."

### 1.2 R^2, T^2, F, and d^2 are one statistic on four scales

Regress the within-stratum-centered treatment indicator on the within-stratum-
centered covariates (strata fixed effects partialled out). Call the resulting
coefficient of determination R^2. Then, writing T2 for a Hotelling/Mahalanobis
quadratic form in the same within-block-adjusted difference vector d that d^2 uses,

    R^2 = T2 / (T2 + df),        df = n - S   (T^2 form)
    R^2 = K F / (K F + (n-S-K)), F the test of "all covariate coefficients zero"

with n the sample size, S the number of strata, K the covariate rank. So R^2, T^2,
F, and the HB08 d^2 = d' V_d^+ d are monotone transforms of one another. (Proved,
verified to 1e-17 in the no-strata case and to 8 digits with strata; the no-strata
df is n-2, the S=1 special case.)

The consequence for inference. Reading observed R^2 against its within-block
re-randomization distribution gives a percentile. Because R^2 is a monotone
transform of the same quadratic form, that percentile equals the d^2 percentile.
R^2-against-its-randomization-distribution is not a second, independent check; it
is the d^2 test wearing a [0,1] costume.

One caveat the adversarial pass forced, and it matters. The monotone equivalence
is *exact* (identical percentile) only when

  - K = 1 (a single covariate), always; or
  - the per-block weight w_b = n1_b n0_b / (n_b (n_b - 1)) is constant across
    blocks -- which holds for any pure pair match (w_b = 1/2 for every pair) and
    for any design where all blocks share a size and a treated:control ratio.

The reason: R^2 uses the unweighted within-stratum scatter Sx = sum_b S_xb as its
metric, while d^2 uses V_d = sum_b w_b S_xb. These agree up to a permutation-fixed
scalar only when w_b is constant. For K >= 2 with heterogeneous block sizes (full
matching mixing pairs, triples, 1:3 sets), the two are no longer exact monotone
transforms: in a heterogeneous K=2 fixture the agent found genuine rank inversions
on 2132 of 5000 permutations, Spearman 0.99984, and percentiles that differ (0.137
vs 0.133; gaps up to 0.03 across fixtures). So: exact for pair matching and other
homogeneous designs; equal-to-simulation-precision otherwise. If you want exact
equality to HB08 in a heterogeneous design, compute the R^2 in the V_d metric
(d' V_d^+ d / SST) rather than the regression's default scatter.

### 1.3 The bad news: raw R^2 does not escape the collapse

Your scale-invariance result applies to R^2 unchanged, and transparently. Write the
explained sum of squares as ESS = d' (Xc'Xc)^+ d and the total as the variance of
the treatment contrast, which does not involve X at all. Under a within-stratum
rescaling Xc -> delta Xc, the difference d scales by delta and (Xc'Xc)^+ by
delta^-2, so ESS is unchanged -- for *every* assignment. The whole randomization
distribution of R^2 is therefore unchanged, and so is the percentile. (Derived;
verified: R^2 = 0.7554 identical at delta = 1, 0.5, 0.1, 0.001; the within-block
percentile flat at 0.0003 across shrink in {0.5, 0.2, 0.05, 0.01}.)

So the within-stratum-metric R^2 buys a scale, not information. It is exactly as
blind to absolute imbalance as d^2 is. The only escape is the one the team already
found: compute the R^2 (or its Mahalanobis numerator) in a *fixed external* metric
-- full-sample cov(X) from a pre-match pool. Then the value falls as delta^2 as the
match tightens (verified: R2_fixed 0.013 -> 0.0016 across the shrink ladder), so it
sees absolute balance. But that fixed metric cancels in the within-block
permutation percentile, so, exactly as in the sigma_x analysis, it changes the
reported *number*, not the *test*. R^2 relocates the scale question onto a
friendlier axis; it does not remove it, the same conclusion the e-values memo
reached for e-values.

A warning for the printout. The within-metric R^2 *value* (as opposed to its
percentile) rises toward 1 as the match tightens (0.18 -> 0.99 across the ladder),
because it is the delta^0 ratio of a shrinking d to a faster-shrinking metric.
Do not display the within-metric R^2 value as a balance-quality score: larger means
*better relative detection of imbalance*, the opposite of the experiment intuition.
Only the percentile is the inferential quantity, and it is flat.

### 1.4 Your two reference standards are three benchmarks

You named two standards. The first, randomization, is the percentile above: it does
exactly what the p-value does. The second -- "the R^2 you would see if the
covariates were constant within strata, so treatment can explain nothing" -- is the
perfect-stratification ideal. If X is constant within every stratum, within-stratum
centering sends X to zero, d = 0, and R^2 = 0 for any assignment. (Proved; verified
R^2 = 0 to ten digits over 500 reassignments.)

There is a third benchmark between them, and it is the one your model-selection
instinct is reaching for. Under the within-block null, the *chance* level of R^2 is
not 0. K covariates explain, on average, K of the n-S within-stratum residual
degrees of freedom of the treatment contrast:

    E[ R^2 | within-block null ] approx  K / (n - S).

(Simulation, matched across four (n,S,K) configurations: 0.068/0.050/0.098/0.033
observed vs 0.067/0.050/0.100/0.033 predicted. This is K/(n-S), the residual-df
form -- not the often-quoted K/n, which ignores the degrees of freedom the block
means consume.) So the three benchmarks order as

    ideal (0)  <=  chance (K/(n-S))  <=  observed R^2,

and observed exceeds chance only when a real residual imbalance survives the
matching. The quantity an analyst wants -- how far observed balance sits above
chance -- is exactly the adjusted R^2,

    adjusted R^2 = 1 - (1 - R^2) (n - S) / (n - S - K),

which centers at 0 under the null (verified: null mean -0.0006 with the block-
honest df, vs +0.008 if you wrongly use n-1) and rises under genuine imbalance.
This three-line picture -- ideal floor, chance floor, observed -- is a clean,
design-independent way to plot the resolution profile that a reader can interpret
without simulating a reference. It is worth adding to the resolution-profile memo.

### 1.5 Where the model-selection instinct genuinely pays off

The adversarial novelty check was blunt: the *linear* R^2 framing reinvents
received wisdom (the Hotelling-T^2-to-R^2 identity is textbook; "regress treatment
on covariates as a balance check" is the propensity-score diagnostic and the
classifier two-sample test). If we write any of it up, those sources go in the first
paragraph, or we repeat the "reinventing a settled idea" problem the paper-decision
memo already flagged. Three pieces survive as genuine, if narrow, contributions:

N1. The chance floor K/(n-S) for the *blocked* d^2 omnibus. Correct, useful as the
benchmark line on an R^2 plot, and apparently unstated in this literature -- but a
corollary of the standard ratio-of-quadratic-forms argument, so claim it modestly.

N2. High dimensions. When K approaches n-S, the d^2 pseudo-inverse truncates at
rank min(K, n-S) and the chi-square reference on K df is simply wrong (verified:
d^2 saturates at 24 with rank 24 when n-S = 24, regardless of K = 22 or 28). The
current d^2 has no good high-dimensional story. A ridge- or cross-validation-
regularized R^2 stays finite and can be calibrated against the within-block
permutation distribution (verified: ridge R^2 well-behaved at 0.26 -> 0.94 as K
goes 5 -> 28). This is a real practical win; credit the shrinkage-discriminant
literature (Ledoit-Wolf; ridge LDA) for the idea and claim only the wiring into the
within-block balance test.

N3. Non-additive imbalance, the most defensible contribution. The linear d^2 is
blind to imbalance that lives in second moments or interactions: with equal
within-block treated/control means but treated SD 2 vs control SD 0.5, the linear
d^2 gives a permutation p around 0.7 (misses it), while d^2 on X augmented with
squares and cross-products gives p around 0.01 (catches it). The squared feature
turns a variance difference into a mean difference the linear machinery then sees;
a random forest or gradient booster discovers such features automatically, and its
held-out AUC or Tjur's coefficient of discrimination summarizes the separability.
The new sentence for the HB08 lineage is precise: the linear omnibus cannot see
non-additive imbalance by construction, and flexible features recover it inside the
within-block randomization null. Credit the classifier two-sample test (Friedman
2004; Lopez-Paz and Oquab 2017; Gagnon-Bartsch and Shem-Tov 2019) and Tjur (2009).
(All citations above need metadata verification before they leave the desk.)

### 1.6 Recommendation for Idea 1

Do not add a linear penalized or cross-validated R^2 as a rival omnibus *test*: it
is provably the same test as d^2, and offering it as a second check would mislead.
Two narrower additions are worth it.

  - Report within-block R^2 or adjusted R^2 as an interpretable effect size beside
    the d^2 p-value: "the covariates explain X% of the within-block treatment
    variation, against K/(n-S) by chance." This is descriptive, like the sigma_x
    natural-scale magnitude the decision memo already endorses, and it is the
    cleanest way to honor your "interpretable scale" goal.
  - Offer a nonlinear or high-dimensional variant as an explicitly *different*
    diagnostic, documented as targeting a different alternative (non-additive
    imbalance; the K-near-(n-S) regime where the chi-square df is unreliable), with
    a within-block permutation reference. Warn that under the design metric it too
    is scale-invariant: it answers "does within-block as-if-randomization hold,
    including in higher moments?", not "is absolute balance good enough?". For the
    latter the user must supply a pre-match pool -- the sigma_x external-scale path.

Net: the R^2 idea is a presentation layer over d^2 plus two real extensions
(high-dim, non-additive). It does not rescue the design-quality question; that
remains the resolution profile's job, which is Part 2.

---

## Part 2. Writing about comparing to designs you did not choose

### 2.1 The worry, and why it dissolves

You chose a paired design on purpose, and l* compares it against a 2:2 design and
complete randomization -- designs you did not choose. The worry is that this looks
arbitrary or even like reference-shopping. It dissolves once we notice you already
published the move, at PNAS, and it passed review without comment.

In the vaccination social-norms study (Rabb, Bowers, Glick, Wilson, Yokum 2022,
PNAS 119(29) e2118770119) you matched survey respondents into pairs, reported the
Hansen-Bowers omnibus as chi^2(7) = 5.1, P = 0.65, and wrote that "the
covariate-to-perceptions relationship in our nonrandomized paired design is
consistent with what we would see in a pair-randomized experiment." That
pair-randomized experiment was never run. The sentence compares a design you chose
to an experiment you did not run, and it is a routine thing to write. l* changes one
thing: instead of checking the single closest comparator, it reads a whole family
of block-randomized experiments ordered by block size and reports the finest one
the design matches. The strange-looking move is the accepted move, graded.

### 2.2 The legitimacy argument

Three established facts carry it, none of which you have to defend afresh.

First, the comparison is licensed by the logic of matching itself. A match is, by
its own stated purpose, an attempt to reconstruct an experiment you could not run.
Rabb et al. say exactly this before the omnibus appears: "We build on the intuitions
from this idealized pair-randomized scenario in creating our own nonrandomized
design." If the design's whole justification is "I built it to imitate experiment
E," then comparing its balance to E's expected balance is not importing a foreign
standard -- it is checking whether the design did the one thing it was built to do.
Hansen and Bowers (2008) state the same premise: the d^2 statistic asks whether the
imbalance we observe is larger than the imbalance randomizing treatment within the
strata would produce.

Second, the comparison is a reference distribution, not a claim about what happened.
A z-score reads a number against the standard normal curve to learn how unusual it
is; nobody objects that the data were not drawn from a normal distribution, because
the curve is a measuring instrument, not a description of the data. The family of
block-randomized experiments plays the same role here, and its graduations -- block
sizes -- are known qualities of balance: fine blocks of nearly identical units
balance covariates well, one coarse block of everyone balances them poorly, and the
sizes in between trace the range. We never assert we ran any of them; we use their
known balance properties to grade the balance we did achieve.

Third, testing against a *specified* design, not just "an experiment," is already
established. Branson (2021, Observational Studies) tests a matched dataset against a
fully specified assignment mechanism, including block randomization, and reads off
whether the data are consistent with it. Once a design can be tested against one
specified randomized standard, testing it against a row of standards ordered by
block size is the same claim repeated; l* is a summary of where, along that row,
the "consistent with" verdict turns over.

The honest boundary, which the precedent observes scrupulously: each of these reads
a non-extreme result as *consistency*, never proof. Rabb et al. write "consistent
with," "compares favorably with," not "is balanced." l* inherits both the
permission and the boundary: it licenses "the design balances the measured
covariates as well as experiment l* would," and forbids "the design is balanced" or
"the design is valid."

### 2.3 Two new facts from simulation

The adversarial pass, role-playing a Rosenbaum-style skeptic, produced two results
that change how the writeup should read. Both are simulation, from-scratch, matching
the randomization_cov_d arithmetic.

A strength to state out loud: l* is invariant to who is paired with whom. Rosenbaum
(2025, Ch. 6, p. 151) objects that "covariate balance refers to the distribution of
covariates in treated and control groups, not to who is paired with whom." The
agent built three radically different valid pairings of the same 48 units --
nearest-neighbor, random-within-group, and adversarial farthest-within-group --
holding the treated and control marginal covariate distributions byte-identical, and
computed the resolution profile for each. The profiles agreed to a maximum absolute
difference of 0.001 (Monte Carlo SE per cell about 0.008). At every rung coarser
than the finest, l* depends on the matched-set *means*, not on the pairing inside a
set. So l* reads off the treated and control covariate *distributions* sorted into
resolution bands, which is exactly the object Rosenbaum says balance is about. This
is a direct rebuttal, and it belongs in the paper with the simulation behind it, not
left implicit.

A defect to concede, and fix: l* depends on the coarsening rule. The agent took one
fixed design of 16 tight pairs and built two equally legitimate coarsenings to the
same block counts -- merging sets that are *near* each other on the matching score,
versus merging sets that are *far* (interleaved). Both are valid block-randomized
experiments you did not run. They cross the median at different rungs: nearest-merge
l* about 4 blocks, far-merge l* about 8 blocks, consistently across five seeds. The
far-merge rule systematically flatters the design by a full rung (a percentile swing
of 0.39 at 8 blocks, larger than the bin-offset swing the resolution-profile memo
currently reports). This is the kind of investigator degree of freedom a careful
referee will find, and the current draft understates it.

The fix is the discipline the memo already gestures at, and it is *forced* by the
measuring-instrument framing rather than bolted on. A coarse block in a real
block-randomized experiment holds *similar* units -- it blocks on a covariate. So
the only honest coarsening merges the nearest matched sets on the matching score
used to build the design, and that ladder must be named before balance is examined.
An analyst who coarsens against the score geometry is building a dishonest
instrument. State plainly that l* is interpretable only relative to a pre-specified,
geometry-respecting ladder. With that rule the profile is reproducible and
pairing-invariant (the strength above is what you get *because* you coarsen by
similarity).

### 2.4 The four-step structure, and draft prose

The writeup should reproduce the four-step arc the reader already accepts from Rabb
et al., then extend only the last step. The draft paragraphs below are in your
conventions and flagged as draft; verify the Rabb et al. citation and the Rosenbaum
and Branson page numbers before circulating.

Step 1 -- open with the experiment you could not run.

> [DRAFT] When we build a matched design we are trying to recover, after the fact,
> an experiment we could not run: one that assigns treatment at random among units
> that resemble each other. So the honest question about a finished match is not
> whether the covariates are balanced but whether they are balanced as well as a
> randomized experiment would balance them -- and the whole question turns on which
> experiment.

Step 2 -- report the like-for-like comparison and show it breaks.

> [DRAFT] Re-randomizing treatment within the realized matched sets is the
> experiment most like the design -- the same sets, the same units, only the labels
> reshuffled -- and a single test against it is the natural thing to report. When
> that comparison is favorable it is enough: in a study of vaccination social norms
> we matched respondents into pairs, reported the omnibus balance test as
> chi^2(7) = 5.1, P = 0.65, and concluded the design was consistent with a
> pair-randomized experiment (Rabb et al. 2022). But the like-for-like comparison
> can turn against a good design. As a match tightens it shrinks the within-set
> spread the finest reference measures imbalance against, so a small, persistent,
> substantively trivial difference is driven toward a percentile of zero -- not
> because the design got worse, but because the standard got stricter. The reference
> most like the design is the one that collapses.

Step 3 -- introduce the ladder as more of the same comparison, not a new object.

> [DRAFT] That favorable sentence already compares a design we chose to an
> experiment nobody ran; the pair-randomized experiment is a measuring instrument,
> not a claim about what happened, in the same way a z-score reads a number against
> the normal curve without asserting the data are normal. So compare to more than
> one. Order the block-randomized experiments by block size, from tiny blocks of
> nearly identical units, which balance covariates well, out to one block holding
> everyone -- complete randomization -- which balances them poorly. The block sizes
> are the graduations, and each is a known quality of balance. For each experiment
> compute the design's percentile in that experiment's within-block re-randomization
> distribution of the same d^2 statistic, coarsening by merging the nearest matched
> sets on the score the match was built from. Read off l*: the finest
> block-randomized experiment whose ordinary luck the design's balance matches.

Step 4 -- report l* in the "consistent with" register, with the honesty constraints.

> [DRAFT, model results paragraph] Our paired design balances the seven measured
> covariates about as well as a block-randomized experiment with blocks of roughly
> [N] units would (l* = [N]; Figure X plots the percentile across the full range of
> block sizes, coarsening by similarity on the matching score). At the finest
> comparison the design sits at the [k]th percentile; the curve rises smoothly to
> the [m]th against complete randomization and first reaches the median experiment
> at blocks of [N] units. This is a calibration of measured-covariate balance, not a
> test of unconfoundedness, and we report it beside the sensitivity analysis in
> Section Y, not in place of it. We also checked that l* does not depend on which
> units were paired within sets: three alternative pairings with the same treated
> and control covariate distributions give the same profile.

### 2.5 What a referee will say, and the reply

  - "A percentile against an experiment you did not run is meaningless." Reply in
    the text, do not dodge it: the inference model is the within-block assignment
    mechanism. Conditional on units, covariates, and blocks, the only randomness is
    which units are treated, so the percentile is an exact finite-sample position in
    a well-defined distribution -- the reference distribution, exactly as the normal
    curve is for a z-score. If the paper does not say this, l* reads as a
    superpopulation or bootstrap claim, which it is not.
  - "This is a p-value stopping rule, the Imai-King-Stuart error." Reply: l* is
    computed on a fixed sample at each fixed resolution, no units are dropped, and
    the threshold is a quality label (the median experiment), not a reject/accept
    boundary. The object is the whole curve; l* is one reading of it.
  - "Balance is about distributions, not who is paired with whom" (Rosenbaum p.151).
    Reply with the pairing-invariance simulation of 2.3, and cite Branson (2021) and
    Rosenbaum's own Problem 6.3 (p.159, "not everyone agrees, and you should form
    your own opinion") to show the field treats the choice of reference as open.
  - "l* is an artifact of an arbitrary ladder." This one has teeth (2.3). Concede the
    coarsening-rule dependence, state the pre-specified geometry-respecting rule, and
    report the whole curve so l* is a summary of a fixed object, not a cherry-picked
    rung. Do not headline l* with the curve hidden.
  - Avoid "certifies." The precedent supports "is consistent with" and "balances as
    well as," not "certifies balance." And keep absolute-magnitude language ("the
    imbalance shrinks as the match tightens") out of the l* section: l* is built from
    the scale-invariant within-block percentile and reads *relative* balance; absolute
    magnitude is the sigma_x descriptive number's job, in a different section.

---

## Part 3. The calibration trilemma (added after the 2026-06-16 discussion)

Two threads from the discussion --- "why compare my paired design to a 2:2 or a
completely randomized design I did not build?" and "can the distance from the R^2
baselines be a balance measure?" --- turn out to be one constraint. It reframes the
whole recommendation, so it goes here rather than getting folded silently into
Parts 1 and 2. The reasoning, with runnable code, is in r2-balance-reasoning.qmd.

### 3.1 The distance from the chance floor is adjusted R^2, and a randomization
distribution of it adds nothing

The distance from the chance floor K/(n-S), normalized by the room above it, is
exactly adjusted R^2: adjusted R^2 = (R^2 - K/(n-S)) / (1 - K/(n-S)) (identity gap 0
in code). So the "distance from a well-balanced baseline" Jake reached for already
has a name, and it is the cleanest one-number effect size the R^2 reframing
produces. But a randomization distribution *of that distance* gives no new test: the
baselines K/(n-S) and 0 are constants (they do not depend on the assignment), and
subtracting a constant preserves ranks, so the within-block percentile of the
distance equals the percentile of R^2 equals the d^2 percentile --- verified
identical to the digit (all 0.0002 at the matched-set blocking). The baselines are
analytic reference points (chance = the randomization mean, ideal 0 = perfect
stratification), not the seed of a second test.

### 3.2 The trilemma

Against the within-pair standard --- re-randomize treatment within the realized
matched sets, the standard the matched analysis actually assumes and the one Jake
trusts (the Rabb et al. 2022 comparison) --- *any* calibration collapses as the
match tightens, regardless of the statistic or the metric, because the reference
distribution's spread is built from the within-set covariate spread and tightening
the match shrinks exactly that. (Team-verified: within-block p runs 0.287 -> 0.002
across the tightening ladder; the fixed metric does not change it.) So you cannot
have all three of:

  1. the within-pair standard (re-randomize within your realized sets),
  2. a calibration ("better than X% of randomized versions") that does not collapse,
  3. no comparison to designs you did not build.

Pick two. The l* resolution profile keeps (2) by relaxing (1): a coarser reference
is a weaker form of the same within-similar-units exchangeability assumption, which
is the only honest reading of "as good as a 2:2 experiment" --- not a foreign design
but a looser version of your own assumption. The distance/magnitude work keeps (1)
and (3) and gives up (2): it reports a magnitude (adjusted R^2, or the fixed-metric
sigma_x distance), not a non-collapsing calibration. Equivalence-within-delta is the
same corner with the "how big is too big" judgment made an explicit tolerance.

### 3.3 Where this leaves the recommendation (revised)

Parts 1 and 2 leaned on the resolution profile as the answer to design quality.
Given Jake's stated preference --- the within-pair re-randomization is the standard
that makes sense to him, and comparing to unbuilt designs is the confusing part ---
the recommendation now leads with the corner that keeps his standard:

  - Lead with the **within-pair comparison plus a magnitude**: report the d^2
    p-value (or its R^2 relabeling) and, beside it, adjusted R^2 and the
    per-covariate standardized differences / fixed-metric sigma_x distance. When the
    p-value collapses on a trivial imbalance, the magnitude is the honest reading.
  - Offer an **equivalence test** against a tolerance delta as the principled "how
    big is too big" answer, on the covariate scale, never invoking an unbuilt design.
  - Demote the **l\* resolution profile** to an optional sensitivity analysis on the
    *resolution* of the exchangeability assumption, presented in that language (never
    "a ruler of experiments you did not run"), and only for analysts who find that
    relaxation a useful question. If it is not useful to them, report the within-pair
    comparison plus the magnitude and stop.
  - Do **not** present a randomization distribution of the R^2 distance as a new test
    (3.1).

The coarsening question Jake flagged --- combine each matched set with its *nearby*
sets on the matching distance scale (a set's location summarized by, say, the
average or median within-set distance) --- is exactly the geometry-respecting ladder
in resolution-profile-threebenchmark.R. Whether that relaxation is a question worth
asking is the open decision, not a settled recommendation.

---

## What is established, and what is not

Proved (and numerically verified):

  - adjusted R^2 = (R^2 - K/(n-S)) / (1 - K/(n-S)) exactly; and a within-block
    randomization distribution of (R^2 - constant baseline) has the same percentile
    as R^2, hence as d^2 (all 0.0002 at the matched-set blocking) --- the distance
    is not a new test (r2-balance-reasoning.qmd).

  - R^2 = T^2/(T^2 + df) and R^2 = KF/(KF + (n-S-K)); R^2, T^2, F, d^2 are monotone
    transforms.
  - Exact monotone equivalence (identical percentile) holds for K=1 or for constant
    per-block weight w_b (pair matching, equal-size/equal-ratio designs); otherwise
    monotone to simulation precision (Spearman > 0.9998, percentile gaps up to ~0.03
    in heterogeneous designs).
  - Within-stratum-metric R^2 is exactly scale-invariant; it does not escape the
    impossibility. A fixed external metric makes the value see absolute balance but
    not the within-block test.
  - X constant within strata gives R^2 = 0 (the ideal floor).
  - l* is invariant to within-set pairing (three pairings, identical marginals,
    profiles within 0.001).

Shown by simulation, not proved:

  - E[R^2 | within-block null] approx K/(n-S) (the chance floor); approximate,
    error O(1/n), under a Gaussian-ish fixture.
  - Flexible/squared features detect non-additive imbalance the linear d^2 misses
    (p ~ 0.7 vs ~ 0.01 on an equal-means/unequal-variance fixture).
  - High-dimensional ridge R^2 stays finite where the chi-square df fails; not yet
    shown to have higher power at fixed level.
  - l* depends on the coarsening rule (nearest-merge vs far-merge differ by a rung,
    in the design-flattering direction); the geometry-respecting ladder fixes it.

Open:

  - The exact within-block permutation expectation of R^2 (vs the Gaussian-df
    heuristic); whether a WLS within-R^2 with block weights restores exact d^2
    equivalence in heterogeneous designs.
  - A CV-calibrated permutation reference for the ridge/forest R^2, sized for
    type-I error.
  - The unit of l* for cross-design comparison (units-per-block is not comparable
    across set structures); a band for l*.
  - Citation metadata for everything not re-read this session: Hansen-Bowers 2008
    wording, Branson 2021, Imai-King-Stuart 2008, Lopez-Paz-Oquab 2017,
    Gagnon-Bartsch-Shem-Tov 2019, Tjur 2009, Rosenbaum 2025 page numbers. The only
    quotation verified against its source this session is the Rabb et al. PNAS
    passage.
