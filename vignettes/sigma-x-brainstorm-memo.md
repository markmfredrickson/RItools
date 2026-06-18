# Brainstorm memo: omnibus balance assessment as design-respecting calibration

**To:** Ben, Mark, Josh
**From:** Jake (with Claude)
**Re:** Alternatives to the d^2 "over-rejection" -- what we actually want, and the menu of statistics that can deliver it
**Date:** 2026-06-14

Companion to `sigma-x-team-decision-memo.md` (which settled the d^2 backend decisions). This one is the substantive/methodological direction. Citations were gathered by three parallel literature agents; page-level details flagged "verify" should be confirmed before any of this enters a paper.

## The goal

We use RItools to SEARCH and SCREEN matched designs. We want to say:

- "This design is more balanced than 80% of randomized versions of this design."
- "It would be very strange to see imbalance this large under randomization of this design."

We keep the DESIGN question -- respect the blocks. We are NOT adopting Rosenbaum's (2025, Ch. 6) ignore-the-blocks complete-randomization benchmark, and NOT the King-style pruning view (discard units to hit a threshold), which wastes data. Matching builds a block-randomized-experiment analog; we assess how randomization-like the realized design is.

## The reframe that delivers it

Both target sentences are the two tails of ONE object: the percentile of the design's omnibus imbalance in the WITHIN-BLOCK re-randomization distribution. The released d^2 with the within-block reference already computes it -- the one-sided p = P(T_rand >= T_obs) IS that percentile:

- large p (e.g. 0.80) = "more balanced than 80% of randomized versions of this design" (the good case)
- small p (e.g. 0.002) = "stranger than 99.8% of them" (the over-rejection, read honestly)

So the change is INTERFACE, not a new statistic or reference:

1. report it as a continuous, two-tailed CALIBRATION ("better than X%"), not a reject/accept threshold;
2. compute it by PERMUTATION/EXACT (the backend we are keeping) -- the chi-square percentile is wrong exactly in the small-support (R<1e4) regime where screening happens;
3. always show per-covariate EFFECT SIZES beside the percentile.

Anchor: Branson (2021, Observational Studies 7(2):1-36) formalizes testing a matched dataset against a SPECIFIED design, INCLUDING block randomization -- the design-respecting counterpart to Rosenbaum's complete randomization.

## The honest tension

Under the within-block reference, a persistent same-sign gap is flagged even when it is small on the covariate's natural scale, because it genuinely is not what block randomization produces. That is correct detection, not an artifact. We handle it by (a) reporting effect sizes alongside the percentile, and (b) not optimizing blocks tighter than we would defend -- benchmark against the stratification we actually believe is as-good-as-random. Only abandoning the blocks (Rosenbaum) makes a persistent gap look fine, and that gives up the design question.

## We checked the alternatives on our own examples

We ran the two leading alternatives -- Chen-Small's graph-based tests and Branson's design-targeted randomization test -- on the same fixture that exposes the d^2 over-rejection (K=4 strata, 2 treated each, a gap of 0.5 on x1 = 0.11 full-sample SD; reproducible in `vignettes/chen-small-comparison.R`).

Tight matching (shrink = 0.05), each statistic applied to RESPECT the design (within-stratum centering -- the graph analog of what d^2 does -- plus the within-block permutation reference):

| Statistic (design-respecting) | p-value |
|---|---|
| d^2 (HB08) | 0.0015 |
| Chen-Small CrossNN | 0.0020 |
| Chen-Small CrossMST | 0.0065 |

All three reject the same 0.11-SD persistent gap at essentially the same strength; at loose matching (shrink = 0.5) all three are fine (0.28, 0.49, 0.14). Branson is a framing whose verdict tracks the target design you name: block-randomization target + d^2 IS the 0.0015; complete-randomization target IS 0.90 (our d2_completeRand). The only configurations that do NOT reject either abandon the blocks (Rosenbaum's marginal + complete randomization = 0.90) or fail to condition on them (CrossMST on raw, un-centered covariates = 0.35 -- swamped by between-stratum geometry, low power, not a principled judgment that the gap is small).

So no statistic -- mean-based, graph-based, or distributional -- makes a persistent gap look like balance while respecting the design.

## The deeper point: this is about p-values, not the statistic

The reason no statistic escapes is that the over-rejection is a property of the p-value, not of d^2. A p-value answers "how surprising is this under the null?", which conflates the SIZE of the imbalance with the PRECISION it is measured at. Tighter matching shrinks the within-block reference (more precision), so a fixed 0.11-SD gap yields an ever-smaller p-value -- the number goes to zero because the ruler got finer, not because the imbalance grew (d^2 within-block: 0.28 at shrink=0.5, 0.0015 at shrink=0.05, with the gap held fixed). Every statistic that feeds a p-value inherits this; swapping the statistic cannot fix it.

Two consequences for what we report:

1. Read the percentile as a CONTINUOUS calibration, never a reject/accept threshold. "99.8th percentile" is information, not a verdict.
2. Always pair it with a PRECISION-INVARIANT MAGNITUDE -- the imbalance on a fixed, design-independent scale that does NOT shrink as matching tightens. This is exactly the discarded sigma_x statistic's VALUE: T = d' cov(X)^{-1} d with FULL-SAMPLE cov(X) is a natural-scale Mahalanobis magnitude that stays put (~0.11 SD here) no matter how tight the match. We dropped sigma_x as a TEST because its p-value cancels to d^2's; the p-value critique resurrects it as the right MAGNITUDE. The per-covariate standardized differences the package already prints, and the energy/MMD distance VALUE (as opposed to its p-value), play the same role.

The honest report for a screened design is a pair: "at the Xth percentile of block-randomizations of itself (calibration), and Y on the natural covariate scale (magnitude)." Either number alone misleads; together they separate detectability from importance.

## The statistic is a separate, lower-stakes choice

The calibration framing is statistic-agnostic: any omnibus imbalance statistic can be reported as a within-block percentile. Choose the statistic for its power profile. Ranked for our use:

1. **d^2 (have it).** Mean imbalance; optimal against dense shifts; decomposes into per-covariate z-scores (interpretable). PRIMARY.
2. **ACAT / Cauchy combination** (Liu & Xie 2020, JASA 115(529):393-402; ACAT variant Liu et al. 2019, AJHG 104(3):410-421). Combines the per-covariate p-values; approximately valid under arbitrary dependence; closed-form Cauchy null; powerful against SPARSE imbalance that d^2 misses. Cheap add -- uses the p-values `balanceTest()` already produces. Caveats to document: tail approximation (not exact); degenerates under strong negative dependence. Report ALONGSIDE d^2 (dense vs sparse are complementary, not redundant).
3. **Graph-based: Chen & Small (2022, Biometrics 78(1):202-213)** -- CrossNN / CrossMST, R package `BalanceCheck`. Catches the FULL multivariate distribution (variance, interactions), uses within-matched-set permutation (respects the design), clean asymptotic nulls. Two independent agents converged on this as the most deployable distributional-balance option. Lineage: Rosenbaum cross-match (2005, JRSS-B 67(4):515-530, exact null); Chen-Friedman (2017, JASA 112(517):397-409, chi^2_2 null, location+scale power).
4. **Energy distance / kernel MMD** (Szekely-Rizzo; Gretton et al. 2012, JMLR 13:723-773). Elegant unifier -- d^2 is the Mahalanobis-kernel special case, and the null is an INFINITE weighted sum of chi-squared(1) (the infinite-dim cousin of our finite mixture), so our permutation + moment-matching backends transfer. BUT: extra power over d^2 is modest for the mean-shift alternatives that matter for confounding (Langmore 2025); NO off-the-shelf stratified version exists (a research project; templates: Ho 2024 restricted block permutation; Ozier-Lafontaine et al. 2024 RKHS designs); loses per-covariate interpretability; curse of dimensionality. RESEARCH DIRECTION, not near-term.
5. **E-values / testing-by-betting** (Shekhar & Ramdas 2024, IEEE Trans. Inf. Theory 70(2):1178-1203). A CONTINUOUS, anytime-valid evidence measure with no pre-set tolerance -- the "metric beyond p-values" without an equivalence test. NOT yet applied to covariate balance: a genuine OPEN PROBLEM. The ambitious build; connects to Vovk-Wang p-merging under arbitrary dependence (Biometrika 2020; Annals 2021).

Unifying thread: most of these share the weighted-sum-of-chi-squares null, so the finite-sample machinery is shared; and "combine evidence under arbitrary dependence" (ACAT <-> Vovk-Wang p-merging <-> e-values) is the coherent spine on the inference side.

## Recommended next steps

1. **Near-term, concrete:** ship the design-respecting CALIBRATION interface on d^2 -- within-block permutation/exact percentile, reported two-tailed and ALWAYS paired with a precision-invariant natural-scale MAGNITUDE (the per-covariate standardized differences already printed, and/or the sigma_x statistic VALUE with full-sample cov(X)). The magnitude is the part the p-value cannot give -- the reason sigma_x survives, as a magnitude not a test. This is mostly wiring the backend we already agreed to keep into a calibration-plus-magnitude report.
2. **Cheap complement:** add ACAT as a secondary omnibus alongside d^2 (sparse-imbalance sensitivity), with documented caveats.
3. **Interop, don't reimplement:** point users to `BalanceCheck` (Chen-Small) for full-distribution balance; consider a thin bridge rather than porting.
4. **Research:** an e-value / betting calibration for "design vs its own block-randomization" is the novel contribution if we want to build one.
5. **Caveat for SEARCH:** selecting the best of many matches on the percentile biases it (winner's curse); condition downstream inference on the acceptance/search criterion (Branson 2021).

## Citation hygiene

Well-verified across agents: Branson 2021; Chen-Small 2022; Liu-Xie 2020 (ACAT); Gretton et al. 2012; Rosenbaum 2005; Rosenbaum 2025 Ch.6; Hansen-Bowers (Stat. Sci. 23(2):219-236 -- one agent miswrote 219-248, use 219-236). Verify before a paper: Shekhar-Ramdas volume/pages; Vovk-Wang 2021 Annals pages; Arias-Castro & Pelletier cross-match consistency venue; Zinger et al. 1992 pages.
