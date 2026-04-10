# HANDOFF.md

Handoff from the 2026-04-09/10 session on the `devel-sigma-x-omnibus` branch.
Nothing has been committed or pushed. Everything is in the working tree.

## 1. Key decisions made

### The original plan (implemented and working)

Ben Hansen proposed replacing HB08's omnibus chi-square statistic
`T = d' Cov(d)^{-1} d` with `T = d' Sigma_x^{-1} d`, where Sigma_x is
a covariance of the covariates rather than the permutation covariance.
This was motivated by a criticism from Paul Rosenbaum (relayed by Ben):
HB08 penalizes tight matches that have small residual imbalance, because
the permutation covariance shrinks with the match.

The new test was implemented with four null-distribution backends:
`satterthwaite_finite` (default), `satterthwaite_asymptotic`, `imhof`,
`simulate`. The finite-sample Satterthwaite backend uses second-order
moment formulas ported from the `i113-highermoments` branch (Mark
Fredrickson). All four backends pass 1561+ tests (3036 total suite).
`R CMD check` returns 0 errors, 0 warnings, 0 notes.

### The discovery that changed the framing

After implementing the new test, we discovered through a diagnostic toy
example that **no test built entirely from the matched sample --- neither
HB08 nor the new sigma_x test nor any of seven other test statistics
(Euclidean, L1, Wasserstein-type, Wilcoxon, tanh, max-distance, energy)
--- can distinguish a tight match from a loose match** when the
within-stratum spread is the only thing that changes. The p-value is
constant across all values of the tightness parameter delta, regardless
of which standardizer or reference distribution is used.

This is a structural property (scale invariance / studentization), not a
bug. Two independent AI agents proved the impossibility theorem
formally. A third and fourth agent verified it against seven statistics
and confirmed the impossibility holds for complete-randomization
references as well.

### The current recommendation (in the memo, awaiting team review)

1. **Keep HB08 intact** as a test of the strict randomization null.
2. **Add effect-size reporting** as a first-class output: adjusted mean
   differences in original units, plus standardized mean differences if
   the user supplies a pre-matching candidate pool or external SDs.
3. **Optionally offer a combined chi-square p-value** using
   pool-derived Sigma_x as the standardizer and chi^2_p as the
   reference. There is a reasonable disagreement among the agents about
   whether this p-value is helpful or misleading; the memo presents
   both views.
4. **Retain the sigma_x_test infrastructure** as the computation
   engine behind pool-calibrated effect-size reporting, not as a
   replacement for HB08.

### Key subsidiary findings

- In any pair-matched design, the new test with default Sigma_x equals
  HB08 exactly up to a constant 2/K. They are the same test.
- Rosenbaum's chapter 6 "complete randomization" benchmark, applied to
  the post-matching sample, is caliper-dependent: changing the caliper
  changes which units survive, which changes the reference distribution.
  Applied to the pre-matching pool, it is caliper-invariant but also
  caliper-blind (it doesn't see the match).
- The useful combination is: d from the matched sample (numerator),
  Sigma_x from the pre-matching pool (denominator), chi^2_p as the
  reference. This gives a p-value that tightens with better matching.
- The `satterthwaite_finite` backend genuinely outperforms asymptotic
  backends at small n (verified against exact enumeration on a 14-unit
  example; the gap is 0.044 in p-value). This finding is independent
  of the Rosenbaum criticism and is still useful.

## 2. Files changed and why

### Modified files

| File | Change |
|------|--------|
| `DESCRIPTION` | Version 0.3-5 -> 0.3-5.9000. Title flagged as devel. CompQuadForm in Suggests. RoxygenNote bumped to 7.3.3. |
| `NEWS.md` | New entry for 0.3-5.9000 describing the devel branch. |
| `.Rbuildignore` | Added CLAUDE.md, CLAUDE_CODING.md, .claude/ |
| `R/balanceTest.R` | Four new args (sigma_x_test, sigma_x, null, n_simulate). ~25 new lines for optional sigma_x integration after HB08. Full @param documentation. |
| `R/print.xbal.R` | Header tweak: prints "---Overall Tests (chi-square and sigma_x)---" when sigma_x columns present. |
| `R/utils.R` | subset.xbal forwards the new sigma_x_info attribute. |
| `R/xbal_tidiers.R` | glance.xbal docstring updated to mention new sigma_x columns. |
| `man/balanceTest.Rd` | Regenerated from roxygen. |
| `man/tidy.xbal.Rd` | Regenerated from roxygen. |

### New files (implementation)

| File | Purpose |
|------|---------|
| `R/sigma_x_test.R` (~300 lines) | All internal functions: `sigma_x_pvalue`, `default_sigma_x`, `randomization_cov_d`, `sigma_x_T_moments`, `sigma_x_test`, `draw_within_stratum_z`, `simulate_T_under_null`, `sigma_x_inferentials`. Ported moment machinery from i113-highermoments. |
| `man/sigma_x_pvalue.Rd` | Generated. Internal. |
| `man/default_sigma_x.Rd` | Generated. Internal. |
| `man/randomization_cov_d.Rd` | Generated. Internal. |
| `man/sigma_x_T_moments.Rd` | Generated. Internal. |
| `man/sigma_x_test.Rd` | Generated. Internal. |
| `man/sigma_x_inferentials.Rd` | Generated. Internal. |

### New files (tests)

| File | Purpose |
|------|---------|
| `tests/testthat/test.sigma_x_test.R` | 29 test_that blocks, 1561+ expectations. Covers: pvalue backends, closed-form V_d, default Sigma_x, test stat math (invariance, reduction to HB08, singular Sigma_x), finite-sample moment machinery (verified against exact enumeration), backend comparison (asymp-vs-finite gap, simulate-vs-enumeration, large-n agreement), end-to-end balanceTest integration, within-stratum permutation structural tests. |
| `tests/testthat/helper-sigma_x_test.R` | Exact-enumeration helpers: `enumerate_strata_assignments`, `exact_randomization_dist`, `closed_form_V_d_helper`, `within_stratum_pooled_cov_helper`, two fixtures. |

### New files (memos and analysis)

| File | Purpose |
|------|---------|
| `vignettes/sigma-x-rosenbaum-memo.qmd` | **THE MAIN DOCUMENT.** ~900 lines. The complete analysis of Rosenbaum's criticism, the toy example, the impossibility theorem, the caliper-dependence finding, the seven-statistic diagnostic, and the recommendation. Ready for team review (user will edit first). |
| `vignettes/devel-sigma-x-omnibus-memo.Rmd` | Earlier memo on the four backends (pre-studentization-discovery). Documents the asymp-vs-finite gap at small n. Partially superseded by the .qmd memo but the small-n findings are independently valid. |
| `vignettes/agent-synthesis-final.md` | Synthesis of the four agent reports. Reference material. |
| `vignettes/agent-synthesis-progress.md` | Earlier progress note. Superseded by the final synthesis. |
| `CLAUDE.md` | Codebase guide (from /init). Not committed. |

### Temporary files (not in the repo)

| File | Purpose |
|------|---------|
| `/tmp/rosenbaum-agent-prompt-v4.md` | The v4 agent prompt used for the final round. |
| `/tmp/rosenbaum-agent-prompt-v3.md` | Earlier version. |
| `/tmp/rosenbaum-agent-prompt-v2.md` | Earlier version. |

## 3. Current blockers and open questions

### Blockers

- **Nothing is committed or pushed.** The user explicitly said no push
  until the team reviews. The working tree is on the local-only branch
  `devel-sigma-x-omnibus`.
- **The memo needs human editing** before being shared with Ben and
  Mark. The user will review the writing, especially the impossibility
  theorem section and the recommendation.

### Open questions for the team

1. **Should balanceTest() offer a pool-calibrated chi^2_p p-value?**
   Two views in tension (see memo Recommendation section). Both agents
   think yes; the user is undecided.
2. **Should the default Sigma_x be changed from pooled-within to
   something else?** The current default is the within-stratum-pooled
   cov, which inherits the Rosenbaum criticism. Changing it to
   full-sample cov(X) doesn't fix the p-value (only the statistic).
   The "real" fix is the pool-calibrated path, which requires the user
   to supply a candidate_pool.
3. **Should the within-stratum-pooled sigma_x path be retained at all?**
   It is useful for the finite-sample Satterthwaite backend (which IS
   a genuine improvement over HB08's asymptotic chi-square at small n)
   but does not address the Rosenbaum criticism. It could be kept as a
   second-tier option.
4. **What should the API look like for candidate_pool?** A data frame
   argument? A named vector of SDs? Both? See the recommendation
   section of the memo for a sketch.
5. **Should we add a fifth null backend for complete randomization?**
   The memo argues it is caliper-blind when applied to the post-matching
   sample, but it could still be useful for the specific question "is
   this assignment unusual under any randomization of these units?"

## 4. Important context to preserve

### The toy example

Eight units, two strata of four, two covariates. Controls at (10, 0)
and (10, 8). Treated at (10+delta, 0) and (10, 8+delta). Parameter
delta controls tightness. cov(X) is full rank at every positive delta.
HB08 gives T = 6, p = 0.0498 at EVERY delta. Seven other statistics
also give constant ranks in the permutation distribution across delta.
This toy is the foundation of the impossibility argument.

### The impossibility theorem

Any test statistic that rescales by the same factor for every
hypothetical assignment when within-stratum deviations are shrunk gives
an invariant p-value. The proof is one line: multiplying every number in
a list by the same positive constant preserves rankings.

### The caliper-blind vs caliper-variant distinction

The agents' key new finding: pre-matching-pool references are
caliper-invariant only in a caliper-BLIND sense (they don't see the
match). The useful combination is pool for Sigma_x + matched sample for
d + chi^2_p as the reference. This gives a p-value that rewards better
matching.

### The pair-matched equivalence

In any pair-matched design (n_s = 2 for all strata), the new test with
default Sigma_x equals HB08 up to the constant 2/K. Verified at machine
precision.

### The finite-sample Satterthwaite finding

At small n (e.g., 14 units), the asymptotic chi-square reference is
off by 0.044 in p-value compared to the exact enumeration (0.338 vs
0.382). The finite-sample Satterthwaite backend gets the moments right
but the chi-square shape is still wrong at n = 14. Only simulate gives
the exact answer. This finding is independent of the Rosenbaum criticism
and is the main contribution of the sigma_x_test infrastructure at the
implementation level.

### The Rosenbaum chapter

`~/REVIEWS/rosenbaum_2025_chap6.pdf`, 13 pages. Chapter 6 of
*An Introduction to the Theory of Observational Studies* (2025).
Advocates comparing balance statistics to their empirical distribution
under repeated complete randomizations of the same units. The chapter's
approach is caliper-dependent when applied to the post-matching sample
(a finding from this session that does not appear to be acknowledged in
the chapter).

### Related work in block_test_power

`~/repos/block_test_power/` has a working paper on power analysis for
block-randomized tests. The test-statistic catalog there (six scores:
raw, rank, mean_dist, mean_rank_dist, max_dist, tanh) was verified to
also be scale-invariant on the toy. The connection: those tests work for
outcome analysis because the outcome scale is baked in by the
experiment's measurement conventions; on the balance-testing side the
analogous scale is the covariate scale, which studentization discards.

## 5. What's done vs. what remains

### Done

- [x] Implementation of sigma_x_test with four null backends
- [x] Integration into balanceTest() (opt-in via sigma_x_test = TRUE)
- [x] 29 test blocks, 1561+ expectations, all passing
- [x] R CMD check: 0/0/0
- [x] print, glance, subset methods all work with new columns
- [x] man pages generated for all internal functions
- [x] Discovery: the new test doesn't fix Rosenbaum's criticism
- [x] Impossibility theorem (two independent proofs)
- [x] Seven-statistic diagnostic (including energy stat)
- [x] Caliper-dependence finding
- [x] Caliper-blind vs caliper-variant distinction
- [x] Main memo (sigma-x-rosenbaum-memo.qmd) incorporating all findings
- [x] Agent synthesis (agent-synthesis-final.md)

### Remains

- [ ] User reviews memo prose before sending to team
- [ ] Team (Ben, Mark, Josh) reviews the recommendation
- [ ] Decision on chi^2_p p-value (offer it or not?)
- [ ] Decision on candidate_pool API shape
- [ ] Decision on whether to keep the pooled-within Sigma_x path
- [ ] Implementation of candidate_pool argument and effect-size table
- [ ] Implementation of pool-calibrated chi^2_p p-value (if approved)
- [ ] Regression tests for the impossibility (toy gives constant p)
- [ ] Regression test for pair-matched equivalence (T_HB / T_new = 2/K)
- [ ] Regression test for pool-calibrated delta-monotonicity
- [ ] Update CLAUDE.md to reflect the new code and findings
- [ ] Commit and push (after team review)
- [ ] Open GitHub issue for tracking
