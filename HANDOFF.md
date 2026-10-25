---
next: Decide with Ben Hansen the four design questions in section 6 of HANDOFF.md (the magnitude report's user-facing form; whether to print size and precision; the 0.25 default and its denominator; two sim-results calls), which need decisions rather than code.
kind: code
project_kind: R package
milestone: the four design questions for the sigma_x branch decided with Ben Hansen
sessions_to_milestone: 1
blocked_on:
estimated_by: session 2026-10-25
---

# HANDOFF: the two-denominators result, the impossibility pressure-test, and the package plan

## 0. Session of 2026-10-25: main merged, truncated Cauchy combination, riposte level tests

This section is the newest. Sections 1 to 8 below are from June 2026 and
still describe the research results and the open design questions. Section
3 describes a workflow that has since finished; its memo is
vignettes/impossibility-pressure-test-memo.md.

### Decisions Jake made this session

1. Jake chose to merge main into this branch. Ben Hansen's pull request #143
   on main fixed
   the svg margin of the base-graphics plot, trimmed the `tcov` attribute to
   the rows of `$results` and gave it dimnames, and made `balanceTest()`
   message rather than stop when every covariate is constant. Commit cfe6fc7
   merges it. Conflicts were resolved this way:
   - `R/Design.R`: the branch's `HB08()` was kept, because its keep/abstain
     block already does what main's all-constant guard does.
   - `R/balanceTest.R`: main's ordering was taken, except that the Cauchy
     combination (ACAT) block sits ABOVE the `p.adjust()` loop. Main's loop
     now overwrites `descriptives[, "p", ]`, and ACAT must combine the
     unadjusted p-values. On a 40-unit example, combining the Holm-adjusted
     ones moved `cauchy_comb_p` from 5.3e-05 to 1. Test 9 in
     `tests/testthat/test.cauchycomb.R` checks that `cauchy_comb_p` is the
     same under `p.adjust.method = "holm"` and `"none"`.
   - When a stratification has df 0 (no covariate left to test),
     `balanceTest()` prints one message containing main's phrase "perhaps all
     covariates are constant (within strata)?". Jake chose this over keeping
     the branch's two messages, so Ben's "Constant variables" test in
     `test.balanceTest.R` runs unchanged and the file matches main. Test 5 in
     `test.omnibus-degeneracy-screen.R` checks that exactly one message prints.
   - Version 0.3-5.9004, one above main's 0.3-5.9003.
2. `acat_pvalue()` uses the truncated Cauchy conversion of Gui, Jiang, and
   Wang (2025, Biometrika 112(4), asaf038) with truncation share t = 0.9, as
   riposte 0.0.0.9008 does. Commit 3f80d9d. Each p becomes
   tan((0.5 - 0.9 p) pi), and with n values summing to S the combined p-value
   is min(1, n P(C > S) / 0.9). Liu and Xie's conversion sent p = 1 to minus
   infinity, so one covariate with an adjusted difference of exactly 0 made
   `cauchy_comb_p` about 1 (0.99999999999999933 in a 40-unit example where
   the truncated version gives 4.6e-05). `truncation = 1` keeps Liu and Xie's
   version. The unused `weights` argument is gone.
   - Jake chose t = 0.9 after seeing the rejection rates under the null
     hypothesis in vignettes/truncated-cauchy-memo.md. At 0.05 the rate rises
     (0.0497 to 0.0626 with 40 independent covariates); at 0.01 it stays
     within 0.003 of 0.01. The memo records the reasoning: for a balance
     check, missing real imbalance is the worse error.
   - Jake approved three test changes. Liu and Xie's idempotence and
     uniformity tests now run with `truncation = 1`. The checks that results
     lie in (0, 1) now allow 1, because min(1, ...) can return exactly 1. The
     end-to-end null test in `test.notforCRAN.R` now requires the mean of
     `cauchy_comb_p` to lie between 0.5 and 0.7 (it was 0.615) instead of
     equal 0.5; Jake chose that bound over dropping the check.
3. At Jake's request ("Strengthen. More sims."), I strengthened riposte's
   tests that its combinations reject a true null hypothesis at no more than
   the nominal rate. Branch `stronger-level-tests` in ~/repos/riposte, commit
   ee1780c, pushed; no pull request opened yet. The large-sample test now
   reassigns treatment 10,000 times (it was 1,000) and also checks 0.01. The
   test of the hybrid (the Cauchy combination with the p-value of riposte's
   quadratic test added) with re-randomization p-values uses 5,000 data sets
   (it was 300). On riposte's own designs the truncated combination rejected
   a true null hypothesis at close to the nominal rate: 0.0487 (Cauchy) and
   0.0471 (hybrid) at 0.05 over 100,000 reassignments of the 36 Wallsten and
   Nteta outcomes, and 0.0507 and 0.0495 over 50,000 fresh data sets.
   `dev/level-truncated-cauchy.R` in riposte reruns these.
4. Jake set rules for how to write to him, after he found a 1,500-word reply
   too long. A new
   subsection, "Replies that need a decision from me", is in
   ~/repos/ai_workflow/CLAUDE.md, uncommitted in that repository: open with
   the one decision and a recommendation, put evidence in a memo, ask one
   decision at a time with the AskUserQuestion prompt, and put a memo's
   pending decision in this header's `next:` line, never in `blocked_on:`.

### Files changed this session

- `R/balanceTest.R`: merge resolution (ACAT above `p.adjust()`), the single
  df-0 message, and roxygen for `cauchy.combination` mentioning the truncation.
- `R/cauchycomb.R`: `acat_pvalue(p_values, truncation = 0.9)` and the new
  internal `truncated_cauchy_log()`, ported from riposte.
- `man/balanceTest.Rd`: regenerated for the `cauchy.combination` text.
- `DESCRIPTION`, `NEWS.md`: version 0.3-5.9004 and entries for the merge and
  the truncation.
- `tests/testthat/test.cauchycomb.R`, `test.omnibus-degeneracy-screen.R`,
  `test.notforCRAN.R`, `test.balanceTest.R` (Ben's version from main).
- `vignettes/truncated-cauchy-memo.md`: the evidence and the decision.

The branch `devel-sigma-x-omnibus` was pushed to
github.com/markmfredrickson/RItools, moving from 431401a to the commit that
adds this section. After commit 3f80d9d, `devtools::test()`
checked 3205 expectations with 0 failures and `devtools::check()` gave 0
errors, 0 warnings, 0 notes. The later commits change only the memo,
NEWS.md, and this file.

### Open questions and loose ends

- The four design questions in section 6 are unchanged and are the next step.
- riposte: open and merge a pull request for `stronger-level-tests`.
- ~/repos/ai_workflow/CLAUDE.md: commit the new subsection.
- Liu and Xie's combination, used with t = 1, already rejects too often at
  0.05 when the p-values are correlated: 0.058 with 6 covariates and 0.065
  with 40 when every pair of z-scores has correlation 0.5. The calibration
  test in `test.cauchycomb.R` draws independent p-values and cannot see this.
  Whether to document it in `?balanceTest` or add a test with correlated
  covariates is undecided.
- Local `main` in this checkout is 13 commits behind `origin/main`; run
  `git pull` on main before the next comparison with it.
- Not from this session and left alone: the modified
  `vignettes/sim-results-memo.pdf` and the untracked png, rds, and pdf files in
  `vignettes/`, and `revdep/`.

### Gotchas from this session

- Keep the ACAT block above the `p.adjust()` loop in `balanceTest()`; test 9
  in `test.cauchycomb.R` fails if it moves below.
- The local roxygen2 is 8.1.0. `devtools::document()` rewrites NAMESPACE into
  multi-line `importFrom()` blocks and bumps `Config/roxygen2/version` in
  DESCRIPTION. That churn was reverted this session; only Rd changes were kept.
- `devtools::document()` warns at `R/harmonic.R:13` that a `@param` lacks a
  description. The warning predates this session.

## The June 2026 handoff (sections 1 to 8)

Date: 2026-06-17. Branch: devel-sigma-x-omnibus. Audience: a fresh Claude
(possibly a Remote Control session driven from Jake's iPad while traveling),
plus Jake Bowers (PI; teaching with RItools soon) and Ben Hansen (co-author of
the d^2 test).

This supersedes the earlier handoffs. Their substance is preserved in the
auto-memory notes (read those first): balance-calibration-trilemma,
collapse-calibration-strengthened-impossibility, omnibus-acat-switch-and-screening,
resolution-profile-invention, branson-2021-randchecks,
balance-calibration-failure-mode-checklist, rct-vs-matched-balance-test,
sigma-x-null-backend-default.

Jake's writing rules in /Users/jwbowers/repos/ai_workflow/CLAUDE.md are mandatory
for every memo: plain words over jargon; motivate before method; name the actor,
the criterion, and the rejected alternative for any evaluative claim; no
decorative structural/industrial/security metaphors; theorem-vs-conjecture
honesty; ASCII only (no unicode anywhere). Every memo gets a writing-critic pass
before it is "done."

## 1. THE RESULT THIS SESSION SHARPENED -- it is the DENOMINATOR, not the omnibus

The collapse: the d^2 omnibus (and any percentile against within-stratum
re-randomization) is homogeneous of degree zero in the imbalance. As the match
tightens, the within-stratum randomization variance V_d -> 0, so a FIXED
substantive imbalance earns an arbitrarily small p. A well-matched design gets an
alarming omnibus p. Verified repeatedly.

The new, clearer statement (the "two denominators"): imbalance is a fraction. The
numerator is the gap dbar (treated-minus-control within sets), fixed in covariate
units. The collapse lives entirely in what you divide by:

- FIXED denominator (covariate units, or a pooled pre-match SD): a 0.2-unit gap
  is 0.2 units whether the strata are loose or tight. max|SMD| with a fixed
  pooled-SD denominator uses this. It does NOT collapse. Its only weakness: it has
  no built-in opinion about whether 0.2 is "a lot" -- you must supply that
  (a substantive tolerance, e.g. "5 years of age").
- RE-RANDOMIZATION denominator (the spread of the statistic under within-set
  re-randomization, = sqrt(V_d)): tight sets give z almost no room, so this spread
  -> 0 as the match tightens. The omnibus p uses this. So does "calibrate max|SMD|
  against re-randomizations of the matched sets" -- it is the SAME arithmetic as
  dividing the fixed gap by the vanishing sqrt(V_d). That is why Jake's
  re-randomized-max|SMD| idea inherits the collapse.

The aggregation (one number over all covariates) is NOT the problem. The
reference/denominator is. See vignettes/two-rulers-demo.R (and .rds): same fixed
imbalance, the within-SET re-randomization 95th-percentile of max|SMD| collapses
from 0.031 to 0.005 as the match tightens, so the observed percentile pins at 1.0
and the omnibus p falls to 4e-11, while a FIXED coarse reference (within-GROUP
re-randomization) gives a stable verdict near the 30th percentile.

## 2. THE ESCAPE -- a FIXED reference (Jake's whole-pool CRE idea), and the catch

A reference that does not shrink with the match does not collapse. Two forms,
same idea:
(a) Jake's substantive tolerances ("max age difference 5 years"): fixed, in real
    units, set at design. Non-collapsing. It gives one judgment per covariate,
    not one omnibus number.
(b) Re-randomize against a FIXED, COARSER design (complete randomization on the
    whole pool), not the tight matched sets.

Validated comparison (vignettes/pressure-test-helpers.R; run reproduced this
session):
- design-internal (d^2; max|SMD| vs within-set): COLLAPSE (alarm 0.66 -> 1.000 as
  the match tightens at fixed gap).
- whole-pool CRE, raw max|SMD| or fixed-Sigma Mahalanobis: do NOT collapse, but
  FLOOR-BLIND -- a full 1-SD residual gap barely registers (alarm ~0.16 / ~0.01),
  because beating a coin-flip on the UNMATCHED pool is trivially easy. "Better
  than complete randomization on the pool" is nearly always true and so tells you
  little.
- whole-pool CRE + COVARIANCE ADJUSTMENT (Jake's parenthetical): removes the
  trivial between-group structure from the reference too, so the comparison is
  about the residual you care about. Non-collapsing AND magnitude-sensitive
  (alarm 0.002 -> 0.93 as the gap grows 0 -> 1). This is the contender worth
  building. CAVEAT being tested by the running workflow: the working demo used the
  ORACLE coarse-group membership; the real test is adjustment on OBSERVED
  covariates (regression / estimated score), no oracle.

The three-way bind (the trilemma, in Jake's own terms): you can have at most two
of {a single calibrated omnibus number, a reference INTERNAL to the design you
actually built, no collapse}. The omnibus p takes the first two and pays with
collapse. The 5-year tolerance takes the last two and pays the single-number
convenience. The fixed-coarse-reference percentile takes the first and last and
pays by comparing to a design you did not build.

## 3. RUNNING NOW -- the impossibility pressure-test workflow

vignettes/impossibility-pressure-test.workflow.js (launched this session; may
still be running, or its result may be in hand). It (1) sends 9 inventors to find
a DESIGN-INTERNAL, single-number, non-collapsing, magnitude-sensitive calibration
(each implements its candidate as vignettes/cand-N.R and runs the shared collapse
+ magnitude harness), with 3 adversarial breakers per self-claimed escape, and
(2) compares the external references (raw whole-pool CRE, Mahalanobis, oracle and
ORACLE-FREE covariance-adjusted CRE, fixed-tolerance baseline). Output: a memo
(intended path vignettes/impossibility-pressure-test-memo.md) plus a return value
with impossibility_supported (true if zero inventors survived), survivors, and
the comparison table. The dimensional argument it tests: any design-internal
calibration is a function of (dbar, V_d); a scale-free function depends on
dbar/sqrt(V_d) and blows up as V_d -> 0 (collapse); the only V_d-bounded ones
ignore V_d and are bare magnitudes -- so no single scalar is design-internal,
magnitude-sensitive, AND non-collapsing at once. The collapse is a THEOREM;
"no design-internal escape" is a STRONGLY SUPPORTED CONJECTURE, not a proof.

WHEN IT FINISHES: write its finalMemo to vignettes/impossibility-pressure-test-memo.md,
read it, fold the verdict (esp. the oracle-free covariance-adjusted CRE result)
into the recommendation, and commit it. The cand-N.R and eval-*.R files are agent
scratch -- keep or delete after curating the memo; they are not canonical.

## 4. A CORRECTION (do not lose)

Two different knobs, often conflated:
- ORGANIZE fixed units into more strata (cut a score finer): BENIGN. The summary
  precision drifts DOWN (sample-size/weighting effects), so p RISES as balance
  improves; no collapse. The non-monotonicity Jake saw is partition jitter
  (~8-11 sign reversals along the ladder; see vignettes/surface-regimes.R/.out).
- TIGHTEN the match toward exact matching (units genuinely more similar within
  sets, V_d -> 0): this is the collapse (vignettes/two-rulers-demo.R).
So "fixed units, vary stratification" is the SAFER framing; the collapse is
specifically the V_d -> 0 operation. (An earlier claim that "refining raises P,
it is a race" was imprecise for the organize-fixed-units knob -- P fell there.)

## 5. FILES MADE THIS SESSION (vignettes/ is in .Rbuildignore, not shipped)

- two-rulers-demo.R (+ .rds) -- the denominator/collapse demonstration (sec 1).
- surface-prototype.R (+ .png/.rds) -- first 2D-surface sweep (organize fixed
  units; size x precision decomposition per rung).
- surface-regimes.R (+ .png/.rds/.out) -- regimes (strong/weak/null proxy) x
  ladders (non-nested cut vs nested) x 30 seeds; the correction in sec 4.
- pressure-test-helpers.R -- the validated shared harness (DGP with fixed gap +
  match-tightening; collapse_curve; magnitude_curve; named contenders). Base R,
  self-contained (no package load).
- impossibility-pressure-test.workflow.js -- the adversarial workflow (sec 3).
- cand-N.R, eval-*.R -- workflow agent scratch (not canonical).

## 6. DONE (through 2026-06-18) and WHAT REMAINS

(As of 2026-10-25 the version is 0.3-5.9004 and the suite checks 3205
expectations; see section 0.) ALL the implementation is DONE and pushed to
devel-sigma-x-omnibus. Version 0.3-5.9002. Full testthat suite 3148 pass / 0 fail / 0 error (NOT_CRAN, Monte
Carlo included); devtools::check() Status OK (0 errors, 0 warnings, 0 notes).

DONE:
- SCREEN + ACAT: relative within-stratum variance screen in HB08/HB08_2016 reading
  a new VarFraction slot (computed in alignDesignsByStrata, assignment-invariant);
  balanceTest() messages dropped covariates, lowers df, ABSTAINS gracefully (df 0 /
  NA p) instead of erroring. check_for_degenerate() points to cauchy.combination.
  Tests: test.omnibus-degeneracy-screen.R (14), test.acat-omnibus.R (5).
- MAGNITUDE REPORT, both forms (R/magnitude.R, test.magnitude-report.R, 19 tests):
  (a) balanceMagnitude() EXPORTED -- per-covariate std.diff vs a settable threshold
  (default 0.25, documented as a convention) + global max|SMD|, print method, reads
  $results only, never pre-filters the omnibus (decision D, CONFIRMED by Jake).
  (b) poolCRE_adjusted_percentile() INTERNAL/EXPERIMENTAL -- structure-adjusted
  whole-pool CRE percentile; k chosen data-driven by the gap statistic (k=NULL
  default); centering = c("conservative","uniform").
- BUG FIXES: scale.DesignOptions was unregistered (missing @export) -> base::scale
  dispatched; fixed. tidy.xbal/glance.xbal registered via S3method(generics::tidy/
  glance, xbal) (kept the function exports); broom/generics -> Suggests; Date test
  skip_if_not_installed. HANDOFF.md -> .Rbuildignore (cleared the check NOTE).
- MEMOS: answer memo (4 edits + writing pass); impossibility memo (writing pass);
  sim-results memo RESTRUCTURED around size x precision (A4/B5/B6/provenance cut,
  re-rendered clean).

REMAINS -- design calls for Jake and Ben (not coding):
1. The magnitude report's final USER-FACING form: keep poolCRE_adjusted_percentile
   internal/experimental, or promote it? Wire balanceMagnitude into the default
   print.xbal/summary, or keep it a standalone reader?
2. Whether to print the SIZE (M) and PRECISION (P) numbers in balanceTest output or
   keep them teaching-only; how to report P across several covariate directions.
3. The 0.25-vs-substantive default and the denominator (pooled vs treated-group SD);
   define the "pool" when units are dropped.
4. Two sim-results judgment calls left as-is (commit e888af1): A1's tightest t=0.001
   row now abstains (kept as a teaching point -- trim the sweep to 0.01 if you'd
   rather show only the pure collapse); B3-B7/B9 still use the base-R prototype
   screen from sim-helpers.R rather than exercising the shipped balanceTest screen.

## 7. REMOTE / TRAVEL (set up before leaving)

Continue this work from an iPad/iPhone with NO SSH via Remote Control (Claude Code
v2.1.181 installed; needs a subscription login, not an API key):
- Keep the Mac awake: lid open, plugged in, `caffeinate -dimsu` in a spare terminal.
- From the repo dir: `claude remote-control --name "RItools sigma-x omnibus"`.
  It serves FRESH sessions in this directory with full LOCAL file access (sees all
  uncommitted work and the running workflow's outputs). Connect from the iPad by
  the URL/QR it prints, or find it by name at claude.ai/code, or scan into the
  Claude mobile app.
- A remote session does NOT inherit this conversation's in-memory context -- it
  picks up from THIS file + the memory notes. So a fresh remote session should
  start by reading HANDOFF.md and the memory.
- Backup if the laptop drops off: the branch is pushed to GitHub, so fall back to
  Claude Code on the web (cloud, GitHub-only, no local files) on
  devel-sigma-x-omnibus.

## 8. GOTCHAS

- Use ABSOLUTE paths in Bash (CWD has drifted into vignettes/ before).
- Do NOT attribute platinum/gold to Tukey; do NOT call a small omnibus p
  "imbalance" when SMDs are small; do NOT teach P = 1/(1-eta^2) (wrong: omits the
  sample-size factor).
- The two reference points are real: exact-matching standard (M=0) and
  block-randomized standard (the omnibus reference). Plain names only.
- Comparability on observed X is NOT ignorability; unobserved confounding is a
  separate sensitivity analysis (Rosenbaum Gamma) and the real next threat.
- revdep/ and CRAN-SUBMISSION: leave alone unless preparing a release.
