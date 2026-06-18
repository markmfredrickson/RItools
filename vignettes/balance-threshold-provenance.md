# Where the 0.25 standardized-difference balance rule comes from

Date: 2026-06-17. A short reference note. We had been calling the "flag a
covariate when its standardized difference exceeds 0.25" rule the
"Rosenbaum-Rubin rule." A check of the primary sources shows that the
attribution is wrong, that two research traditions use two different numbers,
and that those same two traditions use two different standard deviations in the
denominator. Before any of this goes into a paper, confirm the two sources
marked "secondary" below against their originals.

## The 0.25 attached to Rosenbaum and Rubin (1985) is a caliper, not a balance cutoff

Rosenbaum and Rubin (1985, *The American Statistician* 39:33-38) propose 0.25 as
a caliper width for matching on the linear propensity score --- how far apart two
units' propensity scores may be and still be matched. They do not propose 0.25 as
a threshold for judging balance after matching. Stuart (2010) quotes them
directly: "Rosenbaum and Rubin (1985b) generally suggest a caliper of 0.25
standard deviations of the linear propensity score." Citing Rosenbaum and Rubin
(1985) for a balance cutoff confuses two different uses of the same number.

## The 0.25 balance cutoff is Rubin (2001), as compressed by Stuart (2010)

Stuart (2010, *Statistical Science* 25:1-21) writes: "the absolute standardized
differences of means should be less than 0.25 ... (Rubin, 2001)." Rubin's own
text does not say 0.25. Rubin (2001, *Health Services and Outcomes Research
Methodology* 2:169-188) states his first condition as a mean difference of less
than half a standard deviation, and for the propensity score rather than for
individual covariates. The flat "0.25" is Stuart's restatement of Rubin, not
Rubin's own number.

## Two traditions use two different numbers

- Matching and political methodology (Rubin 2001; Stuart 2010; Ho, Imai, King and
  Stuart 2007 and the MatchIt package) use 0.25.
- Pharmacoepidemiology and medicine (Normand et al. 2001, *Journal of Clinical
  Epidemiology* 54:387-398; Austin 2009, *Statistics in Medicine* 28:3083-3107;
  Austin 2011, *Multivariate Behavioral Research* 46:399-424) use 0.1. Austin
  (2011): "a standardized difference that is less than 0.1 has been taken to
  indicate a negligible difference in the mean or prevalence of a covariate
  between treatment groups (Normand et al., 2001)."

Cohen's "small" effect size of 0.2 (Cohen 1988) is a separate lineage and is not
the source of either rule, though it sometimes gets cited as if it were.

## Two traditions use two different standard deviations in the denominator

- Stuart and MatchIt divide by the standard deviation of the treated group, fixed
  at its pre-matching value.
- Austin and the medical literature divide by the pooled standard deviation,
  sqrt((s_treated^2 + s_control^2)/2).

RItools computes its `std.diff` with the pooled standard deviation (Design.R
line 643), which is the Austin convention. The Austin convention is the one that
goes with 0.1, not 0.25. So pairing a 0.25 cutoff with RItools' `std.diff`
draws the cutoff from one tradition and the denominator from the other. We
should choose the default cutoff and the denominator on purpose and write down
the reason for each; `covariate.scales` lets a user supply a fixed denominator
(for example the pre-matching treated-group standard deviation) to get the
Stuart convention instead.

## Read the magnitudes, not the p-values

Ho, Imai, King and Stuart (2007, *Political Analysis* 15:199-236) argue that
balance is a property of the matched sample in hand, and that hypothesis tests
and their p-values should not be used to decide whether a sample is balanced,
because whether a test rejects depends on the sample size, which matching itself
changes. That argument supports reading the standardized differences directly
rather than reading the omnibus p-value.

## Confidence

The quotations from Rubin (2001), Stuart (2010), Austin (2009), and Austin (2011)
come from the primary text. The Rosenbaum and Rubin (1985) caliper wording and
the Normand et al. (2001) wording come from reliable secondary quotation; the
originals were behind paywalls at the time of writing. Confirm those two against
their originals before citing them in a paper.
