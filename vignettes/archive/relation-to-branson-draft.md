# Relation to Branson (2021) --- draft related-work prose

Jake Bowers, with Claude. 2026-06-16. Branch: devel-sigma-x-omnibus.
Status: DRAFT for Jake to edit. Companion to r2-and-unchosen-designs-memo.md,
sigma-x-team-decision-memo.md, and resolution-profile-memo.md. Two agents read
the Branson PDF and the randChecks CRAN source; the prose was then drafted in two
framings and adversarially vetted (a referee-who-knows-Branson pass plus a
writing-rules pass). This file holds the synthesized draft and the credit/claim
split it rests on.

Branson (2021) is the closest prior art for this work. The section must foreground
him, or a referee from the same community will read the paper as reinventing a
settled idea. The boundary the prose holds:

**Credit to Branson (do not claim as new):**

- testing whether a matched dataset is consistent with a SPECIFIED design's
  randomization, with the design made an explicit, swappable null H0: W ~ P*(W|X),
  re-run under complete / block-paired / constrained (rerandomization) designs;
- putting several designs on ONE univariate scale (a fixed-covariance Mahalanobis
  distance) and reading which design a dataset approximates by where the single
  observed value falls in each reference density, reading BOTH tails;
- the design-then-analyze workflow and the caution (p. 19) that assuming a precise
  design can harm inference when post-matching bias remains;
- comparing a matched dataset to designs it was not built as --- a published,
  refereed precedent for the move (alongside Rabb et al. 2022).

**Ours (claim, framed as extending Branson):**

- the within-block collapse / exact scale-invariance, which Branson's fixed-metric,
  single-resolution comparison does not exhibit (he uses an unblocked difference and
  never refines the blocking);
- the point that a fixed metric does NOT cure the collapse (the metric cancels in
  the within-block percentile), so the cure is to read a magnitude, not swap metrics;
- reading the pool-covariance Mahalanobis as an absolute MAGNITUDE (= sigma_x) plus
  per-covariate pool-standardized differences read with substantive judgment, rather
  than as a design-selection p-value, in the design-search-then-fixed-analysis
  setting;
- declining the complete-randomization end of the design menu as the
  matching-as-pruning question;
- the resolution profile and l* as an OPTIONAL sensitivity analysis --- a continuous,
  nested generalization of Branson's discrete menu --- not the centerpiece.

---

## Draft prose

Consider what an analyst actually does. She has a pool of treated and control
units, she matches them into sets that are similar on the measured covariates, and
she builds those sets to approximate a block-randomized experiment --- in the
simplest case, a paired experiment in which one unit of each pair would have been
assigned to treatment by the flip of a coin. She then holds that matched design
fixed and analyzes the outcomes as the design she built, conditioning on the sets
she formed. The covariate balance she reports is evidence for the assumption her
analysis rests on: that within sets, treatment is as good as randomly assigned. She
has already specified the experiment her matching imitates, so her question is not
which experiment her data resemble, but whether the match she built supports that
specification well enough on the covariates she can see.

The work we build on here is Branson (2021), who formalized the test of whether a
matched dataset is consistent with a specified design's randomization. Testing
balance against a design-induced reference is already what Hansen and Bowers (2008)
do, conditioning on the within-block design; Branson's sharpening is to make the
design an explicit, swappable null. His H0 is that treatment was drawn from
P*(W | X) for a particular P*, and he re-runs the same test under each candidate P*
in turn: complete randomization, block or paired randomization, and constrained
(rerandomized) assignment. For an omnibus summary he recommends a Mahalanobis
distance, M = (xbarT - xbarC)' [cov(xbarT - xbarC)]^{-1} (xbarT - xbarC), in which
xbarT - xbarC is the overall, unblocked treated-minus-control mean difference and
cov(xbarT - xbarC) = (N / (NT NC)) cov(X) is the full pre-match pool covariance of
the covariates, computed once and held fixed across every design. Because the
metric is fixed, the reference distributions of the several designs can be drawn on
one scale and overlaid, with a single line at the observed M; Branson reads which
design a dataset approximates by where that one value falls in each density ---
preferring the design under which it is most central, flagging a dataset as too well
balanced for one design (evidence for a more constrained one) or too imbalanced to
be plausible under another. He calls his test a generalization of Hansen and Bowers
(2008), with their permutation test a special case in which draws from P* are
permutations of the observed assignment. Making the design an explicit, swappable
null and putting several designs on one fixed scale are Branson's contributions, and
we use both.

One construction choice separates his test from the Hansen-Bowers d^2, and it is
two-fold. Branson's input is the unblocked overall mean difference and his
denominator is the fixed pool covariance; the Hansen-Bowers omnibus is
T = d' V_d^+ d, where d is the within-block-adjusted treated-minus-control
difference vector and V_d is the covariance of d under within-block permutation ---
a quantity that depends on the design and on the realized match, and that shrinks as
the matched sets tighten. So Branson holds one metric fixed and slides the design
into the reference, while the Hansen-Bowers test studentizes a within-block-adjusted
difference by a denominator that moves with the match. Branson treats the
Hansen-Bowers test as a Mahalanobis-distance test and does not remark on this
difference in input and denominator.

That moving denominator produces a behavior Branson's fixed-metric, single-
resolution comparison does not encounter. Studentizing by V_d makes the within-block
percentile collapse as the match tightens: a fixed, substantively trivial imbalance
is driven toward p ~ 0, not because the imbalance grew but because V_d shrank
beneath it. The collapse is an exact scale invariance, and the argument is one line:
any statistic positive-homogeneous in the within-stratum deviations, compared
against a permutation reference drawn from the same matched sample, has a percentile
unchanged when those deviations are rescaled by a positive constant, because the
observed statistic and every value in its reference rescale together and the rank is
preserved. A fixed metric does not cure it. The metric cancels in the within-block
percentile --- rescaling the within-stratum deviations multiplies the observed
statistic and its whole reference by the same factor, whatever covariance sits in
the middle --- so the cure is not a better denominator but a different object: read
the imbalance as a magnitude rather than as a within-block percentile. Branson's M
is outside this result for two reasons worth keeping separate: his statistic is the
unblocked overall difference, not a within-stratum-deviation quantity, and he
compares designs at a single blocking rather than refining it. His fixed metric is
not what protects him; an unblocked input and a single resolution are.

This shapes what we report. In a matched design search the analyst is not choosing
among complete, block, and constrained randomizations; she has fixed the design from
the outset and chooses the match itself by craft --- the per-covariate
treated-minus-control differences within sets, their percentiles across sets, and
substantive judgment on the covariates that matter ("half the sets show no age
difference, but two sets differ by more than ten years --- is that acceptable
here?"). So we report a fixed-pool Mahalanobis magnitude of the same family as
Branson's M --- the pool covariance is his default cov(X) (his Footnote 3, p. 9),
applied to the within-block-adjusted difference d our blocked design produces ---
read as a magnitude rather than a design-selection p-value. Read this way it does
what the within-block percentile cannot: it shrinks as the match improves, because d
falls while the pool denominator stays fixed. We decline the design-menu
calibration, and one end of it in particular. Reading a paired design against
complete randomization of the whole sample asks whether the treated and control
covariate distributions balance as a coin-flip experiment on these units would,
ignoring who is matched with whom; that is the matching-as-pruning question ---
discard units until the groups balance, then analyze as if completely randomized ---
and the analyst we have described does not adopt it. Branson's own caution is
consistent with holding the built design fixed: he warns (p. 19) that "it may be
harmful to assume a precise experimental design or even any design if there are
substantial biases that remain after matching." One can sweep the within-block
reference across a nested ladder of coarsenings and read off the finest blocking at
which the design's percentile crosses one-half; we treat that summary, l*, only as
an optional sensitivity analysis on the resolution of the within-set exchangeability
assumption, and we keep its coarser-block-randomization middle distinct from its
complete-randomization end, which returns to the pruning question.

---

## Verify before circulating

1. Pagination: cite Observational Studies 7(2), 2021, pp. 44-80 (page headers); the
   Project MUSE cover says pp. 1-36. Resolve against the journal TOC. DOI
   10.1353/obs.2021.0031; arXiv:1804.08760.
2. The verbatim quotes (p. 9 "special cases ..."; Footnote 3, p. 9; p. 19 "may be
   harmful ...") come from an agent's read of the PDF; spot-check against the file.
3. The sigma_x / M relationship: the prose uses Branson's pool covariance with OUR
   within-block-adjusted d; Branson carries an N/(NT NC) factor on the unblocked
   difference, and the two magnitudes coincide only in the unstratified case.
   Confirm the definition you want for the reported magnitude.
