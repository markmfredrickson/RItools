# HANDOFF: the two-denominators result, the impossibility pressure-test, and the package plan

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
    units, set at design. Non-collapsing. Cost: per-covariate, not one omnibus
    number.
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

## 6. STILL OPEN from the prior plan (NOT done this session)

1. Four edits to comparing-design-to-randomized-standard-memo.md still pending:
   (a) DELETE the "## A correction about Tukey" section (Tukey is a red herring;
       drop platinum/gold language); (b) "not a gap in our cleverness" ->
       theorem-vs-conjecture wording; (c) "settled, not open" -> same; (d) ensure
       "direction by direction with several (no single scalar P)". Then a
       writing-critic pass.
2. Rewrite sim-results-memo.qmd around the size-x-precision answer (drop old
   flag/dilute/directions/platinum/gold language). Still old framing.
3. Implement (tests-first suites are RED and define the contract):
   test.omnibus-degeneracy-screen.R -> relative within-stratum variance screen +
   "too little variance" MESSAGE + graceful abstention (R/Design.R ~991, R/utils.R
   ~381); test.acat-omnibus.R is mostly green (ACAT shipped in commit a53dc97) but
   1 failure waits on the screen. Then the magnitude report (per-covariate SMD vs
   a settable reference; global max|SMD|). Then make document + check.
4. Get Jake's explicit yes on decision D (keep the magnitude screen OUT of the
   omnibus; never pre-filter -- B9 showed pre-filtering inflates omnibus size).
5. Package questions: ship M and P or teaching-only? define the "pool" when units
   are dropped? 0.25 vs 0.1 plus the denominator choice? Jake dislikes 0.25 as a
   universal rule -- prefers substantive tolerances at design, or a fixed-coarse
   re-randomization reference (sec 2).

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
