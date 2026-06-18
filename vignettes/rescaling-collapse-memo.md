# Re-inflating the collapsed balance test: rescaling, known variances, and the overfitting analogy

*Memo to Jake Bowers. ASCII only. Each claim is tagged with its proof status: PROOF / SKETCH / HEURISTIC / SIMULATION / CONJECTURE. Simulation scripts: `vignettes/rescaling-check.R` (the 2x2 of metric x null) and `vignettes/rescaling-prototype.R` (the proposed rescues --- magnitude, Stein pooling, shrinkage dial, pool null), run from the repo root with `R_LIBS=.local Rscript vignettes/<script>.R`. The Stein-pooling conjecture in sections 5 and 6 has since been simulated; the result corrected the recommendation --- see the updated text and the ledger.*

---

## 0. Direct answer to your question

Short answer: you can re-inflate, and the re-inflated *magnitude* genuinely does not collapse --- but re-inflation alone does not give you a *test*, because the collapse lives in the null, not in the statistic.

Here is the one-sentence reason, which the simulation confirms exactly. A p-value is a percentile: where the observed statistic falls in its own reference distribution. If you keep the within-matched-set re-randomization reference and merely divide both the observed statistic and every value in its reference distribution by the same fixed number Sigma_0, the percentile does not move at all --- a common positive rescaling of a statistic and its own null is invisible to the percentile. So `M = d' Sigma_0^+ d` referred to the within-set null gives you *exactly the same p-value* as the collapsing `T = d' V^+ d`. (Verified: under the within-set null with X fixed, the Spearman correlation between the two statistics across permutations is exactly 1; see section 4 and the script.)

The re-inflation does one real thing, and it is the thing the HANDOFF already calls the "magnitude": `M = d' Sigma_0^+ d` computed *only on the observed assignment* is a precision-invariant number that stays put as the match tightens (in the toy it sits at 1.65--1.72 while the realized `T` climbs from 5 to 18). That is a fix for the *reporting* problem --- "is the imbalance large on a fixed natural scale?" --- but it is not, by itself, a fix for the *testing* problem.

To turn re-inflation into a non-degenerate *test* you must also change the null away from within-matched-set re-randomization to a fixed/pool reference (a large-pool block-randomized experiment with the same stratum structure and harmonic weights). When you do, the p-value stops collapsing (in the toy it holds near 0.20 across all match tightnesses). But that is a change of null hypothesis: you stop testing "as-if-randomized *within the realized matched sets*" and start testing "balance no worse than a large-pool block-randomized experiment." Your stress-test of your own understanding is correct on every point.

So: re-inflation is **both** a genuine fix (for the magnitude) **and** a reframing (for the test) --- and the reframing is exactly the resolution-profile / pool-reference move you already have. There is no third option that keeps the within-set null and escapes the collapse; section 4 connects that to the team's impossibility claim. The overfitting analogy (section 5) is real and productive, and it points to the same two survivors: an out-of-sample (pool) variance, and selective/conditional inference (Branson). Recommendation in section 6.

---

## 1. Setup and notation

Units `i = 1..N`, one binary treatment `z`, covariate matrix `X` (`N` by `K`), matched sets (strata) `b`. Within set `b`: `n_b` units, `n_1b` treated, `n_0b` control, treated fraction `pi_b = n_1b / n_b`. The within-set-centered treatment contrast is `c_i = z_i - pi_{b(i)}`.

The adjusted difference vector and its within-set re-randomization covariance are

```
d   = sum_i c_i x_i                         (K-vector; harmonic-type weights are absorbed into c and the set sums)
V_d = sum_b [ (sum_{i in b} c_i^2) / (n_b - 1) ] * S_{xb}
```

where `S_{xb}` is the within-set covariate scatter matrix (`sum_{i in b} (x_i - xbar_b)(x_i - xbar_b)'`, suitably normalized). This is exactly `randomization_cov_d` in `balanceTestEngine.R`: in the engine, `dv` is the within-stratum treatment variance, `tmat` is `X` centered within strata and scaled by `wtratio` (the harmonic-mean stratum weight from `harmonic_times_mean_weight()`), `ssvar` is the per-covariate diagonal of `V_d`, and the chi-square is built from the SVD of `tmat * sqrt(dv)` --- the matrix square root whose crossproduct is `V_d`. The HB08 omnibus is

```
T = d' V_d^+ d,    df = rank(V_d),    p = P(chi^2_df >= T)   [or a permutation/finite-sample percentile].
```

**The collapse (PROOF, scalar case; SKETCH, general).** Scale every within-set deviation `x_i - xbar_{b(i)}` by a factor `delta` (a tighter match), holding set centers and the persistent treated-minus-control gap fixed *as a separate residual*. Two regimes:

- *Pure scaling* (gap and noise scale together): `d ~ delta`, `V_d ~ delta^2`, so `T = d' V_d^+ d ~ delta^2 / delta^2 = O(1)` --- invariant. In the one-covariate case `T = d^2 / V` and the cancellation is exact. (SIMULATION: the `oned-invariance` chunk in `resolution-profile-memo.qmd` shows `T` constant to six digits across delta = 1, 0.1, 0.01, 0.001.)
- *Fixed residual* (gap fixed, only noise shrinks): `d -> d_0` (a nonzero constant set by the gap), `V_d -> 0` (the within-set scatter that anchors it vanishes), so `T -> infinity` and `p -> 0`. The better-matched design looks worse. (SIMULATION: the toy in `rescaling-check.R` shows `T_VV` climbing 5.4 -> 13.0 -> 17.6 -> 18.0 as `shrink` goes 0.5 -> 0.01.)

Both are the same fact: the within-set null measures imbalance *in units of the within-set noise*, and matching shrinks that noise.

---

## 2. Formalizing the re-inflation idea

"Re-inflate by a known var(d)" means: stop standardizing `d` by the *realized* within-set covariance `V_d` (which shrinks with `delta`), and instead standardize by a *fixed, design-external reference covariance* `Sigma_0` that does not depend on the realized within-set scatter. Define

```
M = d' Sigma_0^+ d.
```

What should `Sigma_0` be? The natural design-respecting choice you proposed: the covariance of `d` under a *block-randomized experiment with the same stratum structure and harmonic-mean weights, but with units drawn from a large pool*, so that the within-stratum covariate scatter `S_{xb}` is the *pool/population* scatter `S_pool`, not the shrunken realized scatter. Concretely,

```
Sigma_0 = sum_b [ (sum_{i in b} c_i^2) / (n_b - 1) ] * S_pool
```

with `S_pool` the within-stratum covariate scatter you would see if the strata were filled by random draws from the pool rather than by tight matching. (HANDOFF section 6 / the `rct-vs-matched` memo note that the design-*independent* version of this is the full-sample `cov(X)`, which equals the complete-randomization `V_d` up to a scalar; that scalar is irrelevant to a percentile, so "full-sample `cov(X)`" and "pool block-randomized `V_d`" define the same standardization up to scale. The *block* version above keeps the stratum structure and harmonic weights, which is the design-respecting refinement.)

`M` is precisely the precision-invariant Mahalanobis magnitude in HANDOFF #6. As the match tightens (`delta -> 0`, fixed residual), `d -> d_0` and `Sigma_0` is held fixed, so `M -> d_0' Sigma_0^+ d_0`, a finite nonzero constant. (SIMULATION: `M_S0` in `rescaling-check.R` is 1.72, 1.68, 1.66, 1.65 across `shrink` 0.5 -> 0.01 --- flat.) **As a magnitude, re-inflation works. PROOF for the one-covariate case (the limit is immediate from `d -> d_0`, `Sigma_0` fixed); SIMULATION for the general case.**

The question is whether `M`, referred to a reference distribution, gives a *test* that does not collapse. That depends entirely on which reference distribution, which is section 3.

---

## 3. The crucial subtlety: separate the statistic from the null

A p-value needs two ingredients: a statistic and a reference distribution. We have two choices of standardizing metric crossed with two choices of null. Here is the full 2x2, with the verified verdict in each cell.

Standardizing metric (rows): **realized `V_d`** vs **fixed `Sigma_0`**.
Reference null (columns): **within-matched-set re-randomization** vs **pool / large-sample block-randomized reference**.

```
                          NULL = within-set re-rand        NULL = pool block-rand reference
                          (X fixed, z permuted in sets)    (units redrawn from pool)
  ---------------------------------------------------------------------------------------------
  METRIC = realized V_d   (A) HB08 as released.            (C) Observed T standardized by the
  (T = d' V_d^+ d)         COLLAPSES. p -> 0 as match           realized (shrunken) V_d, but
                          tightens. This is the bug.            compared to a pool reference
                                                                whose own V_d does not shrink.
                                                                Incoherent: numerator and
                                                                reference use different scatter.
                                                                Does not collapse but is hard to
                                                                interpret; not recommended.
  ---------------------------------------------------------------------------------------------
  METRIC = fixed Sigma_0  (B) M = d' Sigma_0^+ d compared   (D) M compared to its distribution
  (M = d' Sigma_0^+ d)        to the WITHIN-set null.            under the POOL block-rand null.
                          COLLAPSES --- and gives the          DOES NOT COLLAPSE. p stable as
                          IDENTICAL p-value to (A).            match tightens (~0.20 in the toy).
                          A fixed rescale is invisible         This is a real, non-degenerate
                          to the percentile.                   test --- but of a DIFFERENT null.
```

**Verdict, cell by cell (all SIMULATION-confirmed in `rescaling-check.R`; the (A)=(B) identity is also PROOF):**

- **(A) realized `V_d`, within-set null --- COLLAPSES.** The released test. `p` = 0.012, 0, 0, 0 across `shrink` = 0.5, 0.2, 0.05, 0.01.

- **(B) fixed `Sigma_0`, within-set null --- COLLAPSES, identically to (A).** This is the heart of your question and the answer is a clean negative. Under the within-set null the covariate matrix `X` is held fixed and only `z` is permuted, so `Sigma_0` is a *constant* across the entire reference distribution. Dividing `d' (.) d` by the same constant for the observed value and for every permuted value is a strictly monotone (indeed linear) rescaling of the statistic, and a monotone rescaling cannot change a rank, hence cannot change a percentile.
  - PROOF (one covariate): `M = d^2 / Sigma_0 = (V/Sigma_0) * (d^2/V) = (V/Sigma_0) * T`. But that factoring is misleading because under the within-set null `V` is recomputed per permuted `z` for the `T` statistic, whereas `Sigma_0` is fixed for the `M` statistic --- so `M` and `T` are *not* the same function of the permuted assignment. What IS true and decisive: `M` itself, as a function of the permuted `z`, is just `d(z)^2 / Sigma_0`, a fixed positive multiple of `d(z)^2`. Its percentile equals the percentile of `d(z)^2`, the bare squared imbalance. And `d(z)^2` is exactly what the collapse drives to its extreme tail (as `delta -> 0`, every permuted `d(z)` converges to the gap-determined discrete support, and the observed sits at its max). So `M` under the within-set null has the same percentile as the bare imbalance, which collapses.
  - SIMULATION confirmation, the cleanest single number in this memo: under the within-set null with `X` fixed, the Spearman rank correlation between `T_VV = d^2/V_realized` and `M_S0 = d^2/Sigma_0` across 5000 permutations is **exactly 1.0**. They order every permutation identically, so they return the identical p-value. The columns `p_realizedV_withinNull` and `p_fixedS0_withinNull` in the toy output are equal to four decimals (0.0118, 0, 0, 0). **The metric is irrelevant to the within-set p-value.** This is the same fact the `rct-vs-matched` memo recorded ("a fixed metric does not fix the over-rejection; the metric cancels under the within-set reference"); here it is stated as a percentile-invariance theorem.

- **(C) realized `V_d`, pool null --- does not collapse, but incoherent.** You would be comparing an observed statistic standardized by the shrunken realized `V_d` against a reference distribution whose statistics are standardized by the un-shrunken pool `V_d`. The observed and the reference are on different scales by construction, so the percentile is dominated by the scale mismatch rather than by imbalance. Not recommended; listed only for completeness.

- **(D) fixed `Sigma_0`, pool null --- DOES NOT COLLAPSE. A real test.** Now both the observed `M` and every value in the reference distribution use the same fixed `Sigma_0`, and the reference is generated by re-randomizing in the *pool* (so its imbalances reflect pool-scatter variability, which does not shrink with the realized match). SIMULATION: `p_fixedS0_poolNull` = 0.200, 0.207, 0.210, 0.210 across `shrink` = 0.5, 0.2, 0.05, 0.01 --- flat and non-degenerate. This is the survivor.

**Your understanding, restated and confirmed:** re-inflating by a fixed `Sigma_0` gives a precision-invariant *magnitude* (cell B's *statistic*, read on the observed assignment, is flat). But to get a non-degenerate *test* you must also move the null to the pool reference (cell D), and that is a change of null hypothesis. There is no way to keep the within-set null and escape the collapse: cells (A) and (B) are provably the same p-value. **PROOF** of impossibility within this family: any statistic of the form `g(d)` referred to the within-set null has a percentile that depends on `d` only through the within-set null distribution of `g(d)`; standardizing `d` by any `z`-independent (hence permutation-constant) matrix is a relabeling of `g` that preserves ranks and so preserves the percentile. The collapse is a property of the within-set null distribution of the *direction* `d`, and no `z`-constant rescaling touches it. This is the team's impossibility claim, sharpened: it is not merely that "no statistic from the matched sample escapes," it is that *no standardization of `d` by anything held fixed across the within-set permutations* can change the verdict, because such standardizations are percentile-invariant by construction.

(The only standardizations that are *not* percentile-invariant are ones that depend on the permuted `z` --- like the realized `V_d` recomputed per draw. But those are the studentizing denominators, and studentizing is exactly what *creates* the collapse here, because `V_d(z)` shrinks. So the within-set null is squeezed from both sides: fixed metrics are invisible, and the natural data-dependent metric is the poison.)

---

## 4. Where the collapse actually lives, in one sentence

The collapse is not in the statistic and not in the metric. It is in the **null**: the within-matched-set re-randomization reference measures imbalance in self-induced units (the within-set noise that matching minimized), so it certifies "as-if-randomized within these sets" against an ever-shrinking ruler. Re-inflating the statistic moves the ruler's *label* but not the ruler. Only changing the null changes the ruler.

This is why the resolution profile (your `resolution-profile-memo.qmd`) works: coarsening the reference blocks reintroduces *between-set* covariate variation --- the thing matching actually accomplished --- into the null's scatter, and that variation does not shrink with `delta`. Cell (D) is the extreme-coarse end of that ladder (the pool / one-block reference); the resolution profile is the whole ladder from cell (A) at the finest rung to cell (D)-like behavior at the coarsest.

---

## 5. The overfitting analogy, developed and mined

Your instinct that this is an overfitting problem is exactly right, and it is more than an analogy --- it is the same mathematical structure. Matching *fits* the within-set covariate variation (it chooses sets to minimize within-set distance), and then the within-set scatter `S_{xb}` is *reused* as the error estimate `V_d`. That is in-sample variance estimation after a fit chosen to minimize that very variance: the textbook recipe for an anticonservatively small variance. It is the design-based twin of using the fitted residual variance from a model selected to minimize residuals, or of using training error as a test-error estimate.

Formally: the realized `V_d` is `E_pi[ d d' ]` where the expectation is over within-set permutations *of the realized, fitted sample*. The matching minimized the realized within-set scatter, so `V_d` is a biased-down estimate of the variance `d` "should" have under the design you are claiming (a block-randomized experiment). The test divides a fixed signal by a variance you shrank on purpose.

Here is each transferable fix from the overfitting literatures, mapped onto `var(d)` in this stratified randomization-inference setting, with whether it breaks the collapse and what it costs.

**(a) Out-of-sample / design-based variance (use a pool's `var(d)`).** The cleanest analogue of "estimate error out of sample." Use `Sigma_0` = pool block-randomized `var(d)` instead of the realized `V_d`, AND refer to the pool null. This is cell (D). BREAKS THE COLLAPSE (SIMULATION). COST: a change of null hypothesis --- you now test "no worse than a large-pool block-randomized experiment" rather than "as-if-randomized within the realized sets," and you need a pool (the un-matched reservoir of controls, or the full sample) to define `S_pool`. In your applications you usually *have* that pool (it is the set of potential controls before matching). This is the design-based, superpopulation-free version because the "pool" is a finite, real set of units, and the reference randomization is a real block-randomized experiment over them --- no infinite population, no parametric model. CONJECTURE worth checking: with harmonic weights and a large pool, `Sigma_0` has a closed form you can compute without simulation (it is `randomization_cov_d` evaluated on the pool scatter), which would make this nearly free on top of the existing engine.

**(b) Degrees-of-freedom / effective-df correction.** In linear-model overfitting you correct `RSS/n` to `RSS/(n - df)` where `df` counts fitted parameters. The analogue here: matching "used up" degrees of freedom to align set centers, and `rank(V_d)` overstates how many independent covariate directions the within-set null can actually move once the match is tight. HEURISTIC: this does not break the collapse. The collapse is in the *scale* of `V_d` (it shrinks toward a floor set by the residual), not primarily in its *rank*; a df correction rescales the chi-square reference but cannot stop `d_0' V_d^+ d_0 -> infinity` as `V_d -> 0`. Effective-df helps the *high-dimensional* degeneracy (rank deficiency, the separate problem noted in section 2 of `resolution-profile-memo.qmd`), not the collapse. COST: modest, but it solves a different problem. Not the fix you want.

**(c) Regularized / shrinkage covariance with a known target (Ledoit-Wolf-style).** Ledoit-Wolf shrinkage (Ledoit & Wolf 2004, J. Multivariate Anal., a real reference) replaces a noisy sample covariance with `(1-lambda) S_hat + lambda Target`. Their target is usually `(trace/p) I`; *your* target is principled and known --- `Sigma_0`, the pool/block-randomized `var(d)`. So use

```
V_lambda = (1 - lambda) V_d + lambda Sigma_0,    lambda in [0,1].
```

This is genuinely attractive because it makes the choice explicit and continuous: `lambda = 0` is the within-set test (collapses), `lambda = 1` is the pool magnitude/test (cell D, does not collapse), and intermediate `lambda` interpolates. CRUCIAL SUBTLETY (and this is where it gets interesting): whether shrinkage breaks the collapse depends, again, on the *null*, not just the statistic. If you shrink the *statistic's denominator* but keep the within-set null, then `Sigma_0` and the shrinkage are `z`-constant, so by the section-3 theorem the percentile is *still* invariant to `lambda` --- it collapses for every `lambda < 1` and discontinuously jumps only at... no, even `lambda` near 1 keeps `V_lambda` `z`-dependent through the `V_d` term, so the percentile moves continuously, but the within-set null still drives `d(z)` to its collapsing support. HEURISTIC/CONJECTURE: shrinking the statistic alone, under the within-set null, does NOT break the collapse; you must shrink the *null's* covariance too, i.e. generate the reference under `V_lambda` as well, which for `lambda = 1` is just cell (D). So shrinkage is best understood as a *continuous dial between cell (A) and cell (D)*, equivalent in spirit to the resolution profile's caliper `h`. COST: same as the profile --- you must say which `lambda` (or report the curve). VALUE: it gives a covariance-level, rather than blocking-level, parameterization of the same ladder, which may be cleaner to implement on the existing engine (you already build `V_d` and could build `Sigma_0`; convex-combine and feed the result to the same SVD path). PROTOTYPED (`rescaling-prototype.R`) and confirmed: under the within-set null the dial moves the statistic (`T`: 18 -> 1.65 across `lambda`) but not the p-value (0 at every `lambda`), so it is the profile's caliper expressed in covariance form, not an independent fix.

**(d) Empirical-Bayes / Stein pooling of per-stratum variances toward the pool value.** Efron-Morris / James-Stein (both real: James & Stein 1961; Efron & Morris 1975, JASA) shrink noisy per-group estimates toward a common target. Here the per-stratum within-set scatters `S_{xb}` are each estimated from very few units (often `n_b = 4`), so they are noisy *and* biased-down by matching. Shrink each `S_{xb}` toward the pool `S_pool`:

```
S_b(alpha) = (1 - alpha_b) S_{xb} + alpha_b S_pool,   alpha_b larger when n_b small.
```

It is the most design-based of the rescues, because it keeps the within-set null in spirit and corrects a real pathology (each set's scatter is a low-`n`, selection-biased estimate of the scatter that set "should" have). I first conjectured it would break the collapse. SIMULATION (`rescaling-prototype.R`) shows it does not, and the reason is instructive. As `delta -> 0` the realized `S_{xb} -> 0` but `S_b(alpha) -> alpha_b S_pool > 0`, so `V` stays bounded away from zero and the *statistic* `T = d' V^+ d` no longer diverges --- in the toy `T_stein` flattens near 3 while `T_realized` climbs to 18. But the *within-set p-value does not move*: it collapses identically to the released test (0.012, 0, 0, 0 across `shrink` = 0.5 ... 0.01 --- the same column as realized-`V` to the last digit, and the multivariate omnibus behaves the same). This is exactly the section-3 theorem in action. `S_b(alpha)` is built from covariate scatter and the fixed treated/control counts, so it is *constant across within-set permutations of `z`*, and any `z`-constant standardization is percentile-invariant under the within-set null. Stein pooling is a `z`-constant metric; it cannot rescue the within-set test, only bound the statistic. Nor is the result an artifact of the shrinkage weight: sweeping `alpha` from 0 to 1 (`alpha = 1` replaces each set's scatter entirely with the pool's) leaves the within-set p-value collapsing at *every* value while only the statistic's size moves (`rescaling-prototype.R`, omnibus, `shrink = 0.01`: p = 0 at `alpha` = 0, 0.1, 0.5, 0.9, 0.99, 1 while `T_stein` falls from 18.0 to 1.16). No empirical-Bayes choice of `alpha` can change a percentile the metric does not enter, so the tuning question --- which a referee will ask --- does not bear on the test at all; it bears only on the magnitude `alpha` would report.

CORRECTED STATUS: fix (d) yields a precision-invariant *magnitude* --- a bounded number to report, in the same role as cell (B)'s `M` --- NOT a non-collapsing within-set test and NOT a version of cell (D). To get a non-collapsing test you must still change the null to the pool reference (cell D). COST of the magnitude reading: you assert each set's internal scatter should resemble the pool's, a mild within-set superpopulation flavor, with no assumption on outcomes. So fix (d) is a second, more design-based route to the *magnitude*, alongside the fixed-`Sigma_0` magnitude of fix (a) --- useful for reporting, but it was not the test rescue I had hoped, and the impossibility result of section 3 is why.

**(e) Sample splitting / cross-fitting.** Estimate the imbalance `d` on one part of the data and the reference variance on a disjoint part. In the overfitting literature this is the gold-standard cure (the held-out variance is honest). HEURISTIC: this does not map cleanly onto a single stratified design, because matched sets are small and splitting them destroys the within-set structure that defines `d`. You *could* split *strata* (estimate `Sigma_0` on half the sets, `d` on the other half), but that just approximates the pool `Sigma_0` with extra variance and a power loss, and it still requires the pool-null change to avoid the collapse. COST: high (lost power, awkward with small sets), benefit over (a) is nil. Not recommended here, though it is the principle that *justifies* (a): the pool is the honest held-out sample.

**(f) Post-selection / selective inference; the winner's-curse / Branson framing.** This is the deepest connection and the one most faithful to your design. Matching does not just fit within-set variation; it *selects* a configuration (which units, which sets, possibly the best of many candidate matches). Selective inference (Fithian, Sun & Taylor 2014, arXiv:1410.2597 --- real preprint; Lee, Sun, Sun & Taylor 2016, Ann. Statist., real) says: when you test a hypothesis chosen by looking at the data, condition the reference distribution on the selection event. Branson (2021, Observational Studies 7(2):1-36 --- the design-respecting anchor in your HANDOFF) does exactly this for matching: test the matched dataset against a *specified* design, and when the match was *searched*, condition the percentile on the acceptance/search criterion so the winner's curse does not bias it.

The map onto `var(d)`: the collapse is a within-set artifact, but the *search* over many candidate matches adds a second, distinct bias --- you report the percentile of the *best* match you found, which is optimistic. Conditioning the reference distribution on "this is the match my algorithm selected" corrects that. HEURISTIC: selective inference does not, by itself, break the *single-design* collapse (cells A/B), because for one fixed match there is no selection event to condition on. What it fixes is the *search* problem on top of whichever single-design test you choose. COST: you must specify the selection event (the caliper, the acceptance rule), which you usually know. VALUE: it is the correct treatment for the screening use case in HANDOFF #8 ("search/screen matched designs"), and it composes with cell (D) --- condition the pool-reference percentile on the search criterion. RECOMMENDED as the *companion* to (a)/(d), not as a substitute.

---

## 6. Recommendation

**The honest synthesis: "re-inflate by the large-pool harmonic-weighted `var(d)`" is both a fix and a reframing, and the two halves attach to two different deliverables.**

1. **As a magnitude (the fix): ship it.** Report `M = d' Sigma_0^+ d` (and/or the per-covariate standardized differences) computed on the observed assignment, with `Sigma_0` the full-sample / large-pool harmonic-weighted `var(d)`. It is precision-invariant, design-respecting, requires no superpopulation and no tolerance, and is nearly free on the existing engine (`randomization_cov_d` evaluated on the pool scatter rather than the realized scatter). This is HANDOFF #6, and the section-3 simulation confirms it does not shrink. It directly answers "is the imbalance large on a fixed, interpretable scale?" --- which is the part of your question that a number can settle. DEFENSIBILITY high, DESIGN-BASED yes, IMPLEMENTABILITY trivial.

2. **As a test (the reframing): this IS the pool reference / coarsest rung of the resolution profile.** Re-inflation cannot rescue the within-set *test* --- that is proved, not conjectured (section 3). The only non-degenerate test in the family is cell (D), which is a different null hypothesis: "balance no worse than a large-pool block-randomized experiment." You already have the better version of this: the resolution profile, which reports the *whole ladder* from the within-set null (cell A) to the pool null (cell D) and reads off `l*`, rather than committing to the single coarsest rung. So do not build a standalone "re-inflated test" --- it would just be the coarsest point of the profile. Keep the profile.

**Of the overfitting fixes, the prototype (`rescaling-prototype.R`) has now sorted them:**

- **Empirical-Bayes / Stein pooling of per-stratum scatter toward `S_pool` (fix (d)) --- a magnitude, not a test.** I had this as the first-choice prototype on the conjecture that it breaks the collapse. The prototype refutes that for the *test* while confirming it for the *statistic*: pooling each set's scatter toward the pool keeps `T` bounded as the match tightens, but leaves the within-set p-value collapsing (section 5, fix (d); the cause is the section-3 percentile-invariance theorem, since the pooled scatter is `z`-constant). So fix (d) belongs with recommendation #1 as a second, more design-based way to compute a precision-invariant *magnitude*: it plugs into the existing `V_d` assembly (replace `S_{xb}` with `(1-alpha_b) S_{xb} + alpha_b S_pool` before the SVD) and can be made automatic via `alpha_b` from the sampling variance of `S_{xb}`. Report it as a magnitude beside the profile; do not present it as a test. The disclosure is the same mild within-set superpopulation flavor, with no assumption on outcomes.

- **Selective-inference conditioning (fix (f), Branson) for the SEARCH case only.** When you screen many matches (HANDOFF #8), condition whichever percentile you report on the search criterion. This is orthogonal to the collapse and composes with the resolution profile (or with the pool-null test, cell D). Cite Branson (2021) and Fithian-Sun-Taylor (2014) as the anchors.

**What to drop:** df corrections (fix b) solve the rank-deficiency problem, not the collapse; statistic-only shrinkage under the within-set null (fix c without changing the null) is percentile-invariant and so does nothing; sample splitting (fix e) is dominated by using the pool directly.

**Bottom line for the paper's framing.** The clean message is the one your title implies: *the collapse is a property of the within-set null, not of the statistic, so you cannot rescale your way out of the test --- but you can (i) report a re-inflated magnitude that does not collapse (standardize by a fixed pool variance, or equivalently for this purpose by the Stein-pooled per-set scatter --- both bound the statistic), and (ii) to get a non-collapsing test, change the null deliberately by coarsening it (the resolution profile / pool reference), which trades the strict within-set null for a design-respecting block-randomized one, with the trade stated in the open.* Note which side of that line each rescue falls on: Stein pooling and the shrinkage dial are magnitude-side (they leave the within-set null in place and so leave the p-value collapsing); only changing the reference is test-side. Re-inflation is not a loophole in the impossibility theorem; it is the theorem pointing you at the null.

---

## 7. Status ledger

- COLLAPSE of cell (A): PROOF (scalar), SIMULATION (multivariate, both this memo and `resolution-profile-memo.qmd`).
- Cell (B) = cell (A) identically (fixed-metric re-inflation does not change the within-set p-value): PROOF (percentile-invariance under `z`-constant rescaling) + SIMULATION (Spearman rank cor = 1.0).
- `M` magnitude does not shrink: PROOF (scalar limit), SIMULATION (multivariate flat at ~1.65).
- Cell (D) (fixed `Sigma_0` + pool null) does not collapse: SIMULATION (p ~ 0.20 flat); it is a different null hypothesis (PROOF, by construction of the reference randomization).
- Impossibility of escaping the collapse while keeping the within-set null: PROOF within the family of `z`-constant standardizations of `d`; this sharpens the team's existing impossibility claim.
- Empirical-Bayes Stein pooling (fix d): the earlier CONJECTURE that it breaks the collapse is now CORRECTED by SIMULATION (`rescaling-prototype.R`). It bounds the *statistic* (`T` flattens instead of diverging) but does NOT break the within-set *p-value* collapse, because the pooled scatter is `z`-constant and so percentile-invariant under the within-set null (the section-3 theorem). It is a second route to the bounded MAGNITUDE, not a test rescue.
- Shrinkage dial (fix c) equals the resolution-profile ladder in covariance form: SIMULATION (`rescaling-prototype.R`) confirms it --- under the within-set null the dial moves the statistic (`T`: 18 -> 1.65 across `lambda`) but not the p-value (0 at every `lambda`); it becomes a non-collapsing test only when the null is changed too, i.e. it is the profile's caliper in covariance form.
- df correction (fix b) does not address the collapse: HEURISTIC.

*Citations used are limited to ones I am confident are real: Ledoit & Wolf (2004, J. Multivariate Anal.); James & Stein (1961); Efron & Morris (1975, JASA); Fithian, Sun & Taylor (2014, arXiv:1410.2597); Lee, Sun, Sun & Taylor (2016, Ann. Statist.); Branson (2021, Observational Studies 7(2):1-36); Hansen & Bowers (2008, Statist. Sci. 23(2):219-236). Verify volume/page details against the HANDOFF's citation-hygiene note before any of this reaches a manuscript.*
