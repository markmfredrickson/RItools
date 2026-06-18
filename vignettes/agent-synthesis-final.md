# Agent synthesis: balance tests and the scale-invariance pathology

Final report, two agents + numerical verification. 2026-04-10.

Neither agent could run R (sandbox restrictions). Every analytic claim
below has been independently verified numerically where feasible. The one
empirical question the agents could not close (energy statistic) I verified
myself and report below.

## Where the two agents agree (high confidence)

### 1. The impossibility theorem is correct

Both agents independently wrote clean theorem-proof versions of the
scale-invariance claim. The core statement is the same in both:

> **Theorem.** Let T(X, z) depend on X only through the within-stratum
> centered deviations tilde-X, and let T be positive-homogeneous of
> degree k in tilde-X (i.e. T(delta * tilde-X, z) = delta^k * T(tilde-X, z)
> for every z). Then the rank of T(X_delta, z_obs) inside the permutation
> distribution {T(X_delta, z') : z' in Z} is constant in delta.
>
> **Corollary.** The same holds for any strictly monotone transformation
> of a positive-homogeneous statistic (e.g. tanh of a distance).

Both agents stress-tested the theorem against exceptions and found none.
The only escape routes involve importing external scale information (fixed
kernel bandwidth, fixed bin edges, fixed threshold), which confirms rather
than refutes the memo's taxonomy.

### 2. None of the six candidate statistics escapes

Both agents checked KS, Anderson-Darling, Wasserstein, energy, entropy,
HSIC. Both conclude:

- **KS and AD**: MORE degenerate than HB08 (vertical CDF gap = 1
  regardless of delta, no information about delta at all).
- **Wasserstein**: 1-homogeneous in tilde-X. Falls to the theorem.
- **Energy**: Agent 1 claims it is 1-homogeneous; Agent 2 notes it is
  NOT formally homogeneous because it mixes within- and between-stratum
  distances. **My numerical check settles the question**: rank is constant
  at 1/16 (within-stratum) and 13/28 (complete randomization) across every
  delta. The energy statistic does NOT escape, regardless of which
  theoretical account is correct.
- **HSIC with median-heuristic bandwidth**: Falls to the theorem (bandwidth
  scales with data).
- **HSIC with FIXED bandwidth**: ESCAPES --- but only because the bandwidth
  is external information.
- **Multinomial entropy**: Escapes only with fixed (external) bin edges.
- **Rank-based tests (Wilcoxon, Stephenson)**: Even stronger invariance ---
  they are invariant to ANY monotone transformation of X, not just scaling.

### 3. The pre-matching pool reference: caliper-INVARIANT vs caliper-BLIND

**This is the single most important new finding from this round.** Both
agents independently arrived at the same sharp distinction:

- **The pre-matching pool reference with the statistic computed on the pool
  (what I was calling "Rosenbaum on the pool")** is caliper-invariant --- but
  it is also **caliper-blind**. It tells you about the raw observational
  study, not about the quality of the match. Tightening or loosening the
  caliper does not change the p-value, because the p-value never sees
  the caliper. This is a test of the pool, not of the match.

- **The pre-matching pool reference with the statistic computed on the
  matched sample** is **ill-defined** as a permutation test, because the
  statistic and the reference distribution are computed on different
  sample spaces (the matched N1 units vs. the full N0 units). You cannot
  compute a rank because the observed value and the reference values are
  not commensurable without some embedding convention.

- **The USEFUL combination is: pool-derived Sigma_x as the standardizer,
  matched-sample d as the numerator, and a fixed asymptotic reference
  (chi^2_p) as the comparison distribution.** This is the external-Sigma_x
  path. The numerator moves with the caliper (better match -> smaller d).
  The denominator stays fixed (pool Sigma_x). The p-value is caliper-
  monotone in the right direction. This is a proper test of match quality.

This means: the user's fullmatch-on-pool idea, the Rosenbaum chapter 6
complete-randomization idea, and my earlier "pre-matching-pool reference"
recommendation are all either caliper-blind (if done as permutation tests)
or reduce to the external-Sigma_x path (if done correctly). The memo
should be updated to make this distinction explicit.

### 4. Effect-size reporting is the right path

Both agents endorse effect-size reporting as the primary recommendation.
Agent 1 says: "Reporting the raw mean difference |d|, or any of its norms,
as a number with units --- not as a rank or a p-value --- breaks out of
the theorem because you are no longer comparing it to a permutation
distribution at all."

Agent 2 adds: "chi^2_p with pool-derived V* is just option (a) dressed
up differently. The memo should embrace it as the natural p-value
companion to effect-size reporting." I agree with this: if the user wants
a p-value alongside the effect sizes, the cleanest way to get one is
T = d' Sigma_pool^{-1} d compared to chi^2_p. This requires the user to
supply Sigma_pool (or a candidate_pool data frame from which we compute
it), which is the external information the theorem says is necessary.

### 5. API recommendation

Both agents converge on the same API:

- Keep HB08 as the default p-value for backward compatibility.
- Add a `candidate_pool` argument (a data frame with the pre-matching
  units). When supplied, compute pool-based Sigma_x from it.
- Add an `sigma_x_external` argument for users who want to supply
  their own matrix.
- Add an effect-size table as a new printed section. Columns:
  adjusted mean difference in original units, pool SD (if pool
  supplied), standardized mean difference in pool-SD units (if pool
  supplied).
- Do NOT expose "complete_random_pool" or "fullmatch_pool" as
  permutation-distribution references. They are either caliper-blind
  or ill-defined as permutation tests. The useful version is the
  external-Sigma_x path with a chi^2_p reference.

Agent 1 gives a sample printed output (see below) which is worth
adapting for the package.

## Where the agents differ (lower confidence)

### Energy statistic homogeneity

Agent 1 says the energy stat is 1-homogeneous in within-stratum
deviations. Agent 2 says it is NOT homogeneous because it mixes within-
and between-stratum pairwise distances. My numerical verification shows
the rank is constant, so the practical conclusion is the same. The
theoretical disagreement is about whether the between-stratum distances
matter for the homogeneity argument on this specific toy. Agent 2 is
more careful here: the between-stratum distances are approximately
constant in delta (dominated by the 8-unit gap), so the energy stat is
"approximately" homogeneous even though not formally so. This nuance
does not change the recommendation.

### How harsh to be about chi^2_p

Agent 2 explicitly disagrees with the memo's characterization of chi^2_p
with external V* as the "strongest assumption." Agent 2 argues it is
natural and should be embraced as the p-value companion to effect-size
reporting. Agent 1 does not directly address this but proposes the same
computation in the API (T = d' Sigma_pool^{-1} d compared to chi^2_p).
I agree with Agent 2: the memo should soften this language.

### Whether to expose complete-randomization-on-pool

Agent 1 is emphatic: do NOT expose it. Agent 2 is more nuanced: it could
be a separate diagnostic but should not be the primary path. I lean toward
Agent 1's position: offering a caliper-blind reference as a design-quality
measure is misleading, and the external-Sigma_x path is strictly better.

## The updated impossibility theorem (for the memo)

The sharpest version, synthesized from both agents:

> **Theorem (rank-invariance under within-stratum rescaling).** Let
> X(delta) be a family of covariate matrices with X(delta)_i =
> xbar_{s(i)}(delta) + delta * (X(1)_i - xbar_{s(i)}(1)) for delta > 0.
> Let T(X, z) be a test statistic that depends on X only through the
> within-stratum centered deviations tilde-X = X - stratum means, and let
> T be positive-homogeneous of some degree k >= 0 in tilde-X. Let Z be
> any set of candidate assignments that depends only on the strata
> structure, not on X. Then the rank-based p-value
>
>   p(X(delta), z_obs) = |{z in Z : T(X(delta), z) >= T(X(delta), z_obs)}| / |Z|
>
> is constant in delta.
>
> **Corollary.** The same holds for T = phi(T_0) where T_0 is positive-
> homogeneous and phi is strictly monotone. This covers tanh, log, sqrt,
> and any other monotone transformation.
>
> **Corollary.** Rank-based statistics (Wilcoxon, Stephenson, any
> function of within-stratum ranks) satisfy the stronger property that
> the p-value is invariant to ANY strictly monotone transformation of X
> within strata, not just positive scaling.

## Concrete next steps for the memo and the package

1. **Add the impossibility theorem** to the memo, replacing the informal
   sketch. Both agents provided proof-quality statements.

2. **Correct the "pre-matching pool" section** with the caliper-invariant
   vs. caliper-blind distinction. The pool reference is caliper-blind
   when done as a permutation test. The useful combination is pool-derived
   Sigma_x + matched-sample d + chi^2_p reference.

3. **Add a "recommended API" section** based on the agents' convergent
   proposals. Key elements: candidate_pool argument, sigma_x_external
   argument, effect-size table in the printed output.

4. **Soften the language** about chi^2_p with pool-derived V*. It is not
   the "strongest assumption" --- it is the natural p-value companion to
   effect-size reporting, and the pool provides the scale information the
   theorem says is needed.

5. **Add the energy-statistic numerical verification** to the memo's
   six-test-statistic diagnostic section, as a seventh entry. It is the
   one statistic where the theorem's hypothesis formally fails but the
   rank invariance empirically holds.

6. **For the package code**: the `sigma_x_test` infrastructure on the
   devel branch is still useful as the computation engine for
   T = d' Sigma_x^{-1} d, which is the statistic the external-Sigma_x
   path needs. What changes is the *framing*: it is no longer sold as a
   "better balance test" but as the machinery behind pool-calibrated
   effect-size reporting with an optional p-value.
