# Collapse-method invention: result (claude)

WINNING LENS: graft (A feasible-matchings + E efficiency-ratio + an absolute
fixed-scale magnitude gate), realized most completely in Codex's Protocol-Replay
Feasible Balance Calibration. But the headline finding is a STRENGTHENED
IMPOSSIBILITY: no lens (A-F, in this run or Codex's) yields a single number that
is at once non-collapsing, floor-sensitive (actually answers "is this design
balanced"), and free of a forbidden reduction. The honest deliverable is a
two-number report, not one calibration.

---

## 0. What this file is

A fresh multi-agent search was run on the challenge in
`collapse-method-invention-challenge.md`: 6 inventor agents (lenses A-F) in
parallel, 3 adversarial skeptics per proposal, then a synthesis. This file
reports the result, cross-checked against four standalone R simulations I ran
directly and against an independent Codex proposal
(`protocol-replay-feasible-balance-calibration-codex.md`).

The short version: every proposal broke. They did not break at random. They
broke on a small, recurring set of obstructions that, taken together, sharpen
the challenge's stated trilemma into a stronger claim. I give that claim in
Section 2, the evidence in Sections 3-6, and the constructive thing a researcher
can still do in Section 7.

---

## 1. The collapse, reproduced

Before assessing cures I reproduced the disease, method-agnostically, holding the
observed mean imbalance fixed at about 0.20 and shrinking the within-pair
idiosyncratic spread sigma (a tighter match). Code: a 200-pair design,
within-pair difference drawn at mean 0.20 (the residual) and SD sqrt(2)*sigma.

```
 sigma within_var_X obs_meandiff_d null_var_Vd        T        p
  1.00     0.843609       0.295895    0.008436  10.3785 0.001275
  0.50     0.281646       0.217627    0.002816  16.8161 0.000041
  0.25     0.084250       0.206770    0.000843  50.7460 0.000000
  0.10     0.030068       0.207363    0.000301 143.0092 0.000000
  0.05     0.022041       0.196471    0.000220 175.1318 0.000000
  0.02     0.020649       0.201062    0.000206 195.7750 0.000000
  0.01     0.019863       0.198868    0.000199 199.1074 0.000000
```

The observed imbalance is essentially constant; the null variance V_d shrinks; T
rises toward S = 200 and p collapses to 0. A tighter match with the SAME residual
imbalance is reported as far worse balance. This is the pathology the challenge
names, and it is the yardstick every proposed cure must clear: as the match
tightens with residual imbalance held fixed, the reported number must not
degrade.

---

## 2. The result: a strengthened impossibility

The challenge frames a trilemma --- you cannot keep all of (1) the within-set
standard, (2) a non-collapsing calibration, and (3) no comparison to unbuilt
designs --- and asks for a method that keeps corner (2) by trading (1) and/or (3)
toward analogous designs. The six skeptic panels, working independently, kept
finding that the trade does not buy what it promises. Their breaks reduce to one
structural argument, which I state as a claim and then support.

CLAIM (strengthened obstruction). A calibration that does not collapse needs a
reference distribution whose spread stays bounded away from zero as the
within-set spread of X goes to zero. There are exactly three sources of a
bounded reference spread, and each is independently forbidden or collapsing:

  Source (I) --- the pooled / marginal covariate scale. Any reference built by
  injecting pooled-scale noise into X, by drawing covariates exchangeably across
  the pool, or by scoring fixed matched sets on held-out noise, has reference
  variance Var(dbar) = 2*sigma_pool^2 / S. For 1:1 matching that equals
  sigma_pool^2 * (1/n1 + 1/n0): the complete-randomization-of-the-pool null,
  exactly. This is the comparison the challenge fences off. Two skeptics verified
  the variance ratio numerically at 1.000 (PIC, HOPB), and a third reproduced it
  at 1.0018 with B = 2e5 for the optimism-correction proposal. Source (I) is
  MODE 1.

  Source (II) --- the spread across feasible matchings (the assignment freedom).
  This is design-respecting, but it fails in one of two ways. (a) Matching is an
  optimization, so the analyst reports the most-balanced feasible matching; the
  percentile-among-feasible-matchings then pins at about 1.0 by construction
  (MODE 3). (b) To keep the reference from collapsing, these methods cancel the
  irreducible imbalance I_floor = delta' Sigma^{-1} delta shared by every
  feasible matching; that cancellation makes the number blind to I_floor, which
  IS the residual confounding the balance question is about. What survives the
  cancellation is the removable assignment part --- the matching-as-pruning
  content the challenge also forbids (MODE 4).

  Source (III) --- no reference distribution at all (a sensitivity margin). With
  the optimism correction driven to its limit the margin is a deterministic
  monotone function of the bare fixed-scale magnitude d' Sigma^{-1} d (Spearman
  0.99993, exactly the forbidden "magnitude wearing a percentile costume",
  MODE 2), or it saturates at a data-independent ceiling tau/abar and stops
  discriminating among good designs (MODE 5).

The three sources exhaust the ways to get a non-vanishing spread. Each lands on a
forbidden or collapsing mode. So the balance question about the MATCHED
covariates --- "is the residual confounding on the covariates I matched on small,
calibrated against a reference" --- has no non-collapsing, non-forbidden
calibration. Keeping corner (2) for that question is not possible under the
challenge's own constraints. That is the strengthened impossibility.

THE DICHOTOMY LEMMA (closes the held-out / sample-splitting family). The
sample-splitting and held-out-covariate proposals (lenses B and C, and SOG
reading (ii)) deserve their own argument, because they look like the way out and
are not. Let the held-out reference score balance on an evaluation covariate E.
Exactly one of two things holds. Either E depends on the matching (it is a matched
covariate or a function of them), so the matcher shrank E's within-set spread and
the reference variance vanishes as the match tightens --- the collapse, relocated
into the reference (MODE 5; a simulation drove the in-sample reference variance /
(2/S) from 0.11 to 0.03 across a tightening sequence). Or E is independent of the
matching (genuinely held out), so the realized pairs are independent of E, the
reference variance equals the complete-randomization null (MODE 1), AND the
held-out statistic is a function of E alone and therefore blind to confounding on
the matched covariates (MODE 4; the same simulation: as true matched-X confounding
climbed 0 -> 1.0, the matched-covariate residual and the in-sample d^2 tracked it,
the held-out gap did not). There is no third case: non-collapse of a held-out
reference requires independence of E from the matching, and that independence forces
both the forbidden scale and matched-X blindness. The held-out idea cannot be at
once non-collapsing and an answer to the matched-covariate balance question.

One honesty correction, surfaced by an adversarial check of the proof itself. The
equality "held-out reference variance = 2*sigma^2/S exactly" holds for iid
pool-scale noise (verified ratio 1.000, and 1.0018 at B = 2e5); for a held-out
RE-PAIRING reference over a FINITE control reservoir it is only approximate (the
ratio runs about 0.02 to 0.46 as the reservoir grows, a finite-population
correction). The robust core of the impossibility is therefore the MODE-4
blindness, not the exact variance: an independent held-out covariate cannot see
confounding on the matched covariates however its reference variance is computed.
The variance identity is the clean special case, not the load-bearing step.

The cleanest single sentence comes from the skeptic who broke the
optimism-correction proposal: "the reference that does not collapse is the
reference that is forbidden." Non-collapse and the forbidden comparison are, for
the matched-covariate balance question, the same property.

The impossibility is not a counsel of despair. It tells you which questions DO
have a non-collapsing answer (Section 7): the design-search question and the
absolute-magnitude question. It also tells you to stop spending effort on the one
that does not.

---

## 3. Phase-1 proposals (A-F) and why each broke

One paragraph each. All six broke; the mode each hit is named.

LENS A --- Feasible-Matching Permutation Calibration (FMPC). Re-randomize which
eligible control attaches to which treated unit, score balance across feasible
matchings, place the realized matching as a percentile; evaluate on held-out
covariates to blunt the optimizer. Broke 3/3. The caliper rho that defines
feasibility is a continuous dial: at loose rho the feasible-matching reference
reproduces the complete-randomization reference to Monte Carlo error (MODE 1, a
skeptic measured p_FM = 0.797 vs p_cr = 0.789), and holding ONE realized design
fixed while varying only the reference caliper swung the verdict across the 0.5
line. The held-out evaluation did not stop cherry-picking: optimizing the held-out
statistic directly within the feasible set pushed the percentile to 0.88-0.99
across eight seeds (MODE 3).

LENS B --- Selection-Optimism Gap (SOG). Score held-out-covariate balance against
the balance the same matcher achieves on signal-free covariates. (Its three
first-round skeptics crashed on transient server errors; I re-ran the panel ---
see Section 4.) Broke 3/3 on re-run. MODE 1 under both readings, MODE 2/5 dilemma,
and MODE 3 cherry-pick via the analyst-chosen eval/match split.

LENS C --- Held-Out Placebo-Covariate Balance Index (HOPB). Score the held-out
imbalance against a placebo reference drawn at pool scale. Broke 3/3. The placebo
reference variance is 2*sigma_pool^2/S = the complete-randomization null
(MODE 1, verified ratio 1.000); the ratio is dominated by S so a fixed held-out
residual is reported as imbalanced at S = 400 but "better than chance" at S = 12
(MODE 5, mirror-image collapse); and a fresh noise reference is independent of
the selection that shrank the numerator, so cherry-picking is not cancelled
(MODE 3).

LENS D --- Design-Sensitivity Balance Tolerance (DSBT). Report Gamma*, the
multiple of a fixed per-set bias allowance at which worst-case displaced imbalance
breaches an absolute tolerance tau. Broke 3/3. As the match tightens with a
growing control reservoir the optimism term vanishes and Gamma* = tau -
sqrt(d' cov^{-1} d), a strictly monotone relabel of the forbidden magnitude
(MODE 2, Spearman -1); it saturates at tau/abar so a 20-fold quality difference
among good designs collapses to a 0.02 gap (MODE 5); and Gamma* moves with the
analyst's choice of which controls enter the pool covariance (gameable through
the metric).

LENS E --- Efficiency-Balance Equivalent (EBE). EBE = (I_null - I_real) /
(I_null - I_opt), the fraction of closable imbalance closed. Broke 2/2. The same
cancellation that delivers non-collapse cancels I_floor, so EBE = (R_null -
R_real)/R_null is blind to the irreducible confounding: a catastrophically
confounded design (I_real = 3.23) scored EBE = 1.00 and percentile 1.00, ranking
ABOVE a clean design (MODE 4). And because the analyst reports the optimum,
EBE = 1.000 and percentile = 1.000 by construction (MODE 3).

LENS F --- Persistence-of-Imbalance Calibration (PIC). Resample the covariate
channel at pooled dispersion, recompute d, report the percentile of |d|. Broke
3/3. Its non-collapsing form has reference variance 2*sigma_pool^2/S =
complete-randomization null (MODE 1, verified ratio 1.000000); the deliverable
also specified a set-centered bootstrap whose dispersion is exactly the within-set
spread that goes to zero, so the SPECIFIED method collapses while the simulated
closed form does not (the non-collapse evidence was for a different method than
the one written down).

---

## 4. Proposal B, judged (the repair run)

B (the optimism-correction lens) is the one structurally able to dodge MODE 3 and
MODE 4, so its untested status mattered. I re-ran its 3-skeptic panel. Verdict:
broken 3/0, with concrete numerics.

  Skeptic 1 (MODE 1 + MODE 4). The held-out reference IS the complete-randomization
  null on the held-out covariate, under both readings of "noise-optimism
  reference." Reading (i), fixed sets scored on iid pool-scale noise, has
  Var(dbar_W) = 2*sigma_pool^2/S exactly (ratio 1.0018 at B = 2e5). Reading (ii),
  re-run the matcher to optimize on noise columns and read column-1 balance, does
  NOT escape: matching minimizes within-pair Mahalanobis DISTANCE, which is not
  the across-pair MEAN of signed differences, so the optimizer does not shrink the
  per-column mean difference --- the claimed "optimism floor below the random-split
  center" does not exist for a mean-difference statistic (mean|md| = 0.091 equals
  the un-optimized paired-noise level 0.092; var ratio 0.96). SOG's percentile
  tracks the complete-randomization p-value almost exactly (confound 0: 0.9919 vs
  0.9913; confound 0.3: 0.0988 vs 0.0956). Favorable note, reported honestly: SOG
  is NOT floor-blind --- it detects genuine held-out confounding (percentile
  crashes 0.52 -> 0 as confounding rises). But it detects it only by scoring
  against the forbidden complete-randomization spread.

  Skeptic 2 (MODE 2 + MODE 5 dilemma, plus scope). If the reference is held fixed,
  the percentile is a deterministic monotone transform of the scalar held-out
  magnitude (Spearman 0.99993, zero decreases): MODE 2. If the reference is
  recomputed under the same tightening caliper, the optimizer-on-noise floor
  collapses while the realized held-out residual does not, so the percentile pins
  at 1.000 throughout (MODE 5). No third setting exists. Separately, SOG is silent
  on the matched covariates: a design with 0.58 SD residual on every matched
  covariate scored 0.967, indistinguishable from a clean design's 1.000.

  Skeptic 3 (MODE 3 cherry-pick via the split). Under reading (ii) the reference
  depends only on (n, Z, number of match columns) --- not on X, not on the
  eval/match split, not on the realized matching. So the analyst routes the
  confounded covariates into X_match (the matcher "absorbs" them and they never
  reach X_eval) and reports, as X_eval, the single clean held-out covariate whose
  realized imbalance drew lowest. Reported percentile 1.000 on a design carrying
  d^2 = 28.6 of unaddressed confounding.

The three breaks are independent and each is numerically reproduced. B is broken.
Its one favorable property --- genuine sensitivity to held-out confounding --- is
exactly what Section 2 predicts: that sensitivity comes from scoring against the
pooled scale, which is Source (I), which is MODE 1.

---

## 5. Codex's proposal, incorporated and assessed

Codex's `protocol-replay-feasible-balance-calibration-codex.md` is the strongest
single constructive proposal from either engine. It has two parts:

  - A tolerance GATE: PASS if L_obs = d' Sigma0^+ d <= c, where Sigma0 is a fixed
    pool-scale covariance and c = r*delta^2 for a pre-specified scientific
    threshold delta.
  - A calibration PRFP: the percentile of L_obs among replays of the FULL
    predeclared design protocol --- including the caliper grid, distances,
    candidate algorithms, tuning, and final selection rule --- under a documented
    probability law Q_P over matched-design search paths. An optional efficiency
    ratio SER = (mu_P - L_obs)/(mu_P - L_best_P) accompanies it.

What Codex gets right, and it is a genuine advance over every Claude inventor:
the MODE 3 defense. PRFP's reference replays the SAME search the analyst ran, so
the realized design is not special relative to the reference --- both are
"best-of-the-same-search." That neutralizes the cherry-pick that killed FMPC and
EBE. The price, which Codex states plainly, is that the analyst must have RECORDED
the entire search; an unrecorded specification search makes PRFP a sensitivity
analysis over protocol scope, not a clean calibration.

Where Codex does not escape the impossibility --- it embodies it. PRFP ranks
L_obs among other matchings of the SAME labeled pool, so the irreducible
confounding I_floor, shared by every feasible matching, cancels in the ranking.
PRFP is therefore floor-blind (MODE 4): on its own it cannot tell a clean design
from a confounded one. What restores floor-sensitivity is the GATE --- and the
gate is L_obs = d' Sigma0^+ d against an absolute threshold, exactly the bare
fixed-metric magnitude the challenge forbids "calling a test." Codex does not call
it a test; he calls it a gate and is explicit that "the gate prevents an optimized
design from receiving a favorable percentile while still having scientifically
unacceptable absolute imbalance." That is Codex conceding, in his own design, that
the percentile is floor-blind and that floor-sensitivity has to come from the
magnitude.

I verified the split numerically (vary the irreducible floor delta on a held-in
covariate; Q_P = random feasible matchings within a caliper):

```
 Realized = TYPICAL feasible draw          Realized = OPTIMUM on X_match
 delta gate_Lobs   PRFP                     delta gate_Lobs   PRFP
  0.00    0.0481 0.0623                       0.00    0.0001 0.9900
  0.25    0.0419 0.4115                       0.25    0.0596 0.2943
  0.50    0.2694 0.6608                       0.50    0.4107 0.0673
  1.00    0.7688 0.8753                       1.00    0.7010 0.7681
  2.00    2.4937 0.8354                       2.00    2.5515 0.6833
```

gate_Lobs tracks the true floor monotonically (0.05 -> 2.49) in both columns:
the magnitude is the floor-sensitive number. PRFP does not track it reliably ---
in the optimum column it is non-monotone in the true confounding (0.99, 0.29,
0.07, 0.77, 0.68). So the floor-sensitivity in Codex's report lives entirely in
the gate, i.e. in the forbidden magnitude; the non-collapsing percentile is the
floor-blind half. This is the Section-2 graft made concrete: a usable report needs
BOTH a forbidden magnitude (for floor-sensitivity) AND a feasible-matching
percentile (for non-collapse), because neither half is both.

(The exact floor-blindness of the efficiency RATIO is algebraic, not just
empirical: with I = I_floor + R, the SER ratio (mu - L_obs)/(mu - L_best) has
I_floor in every term and cancels it, leaving a function of the removable
assignment part R alone --- the EBE skeptic's proof, which transfers to SER
verbatim.)

---

## 6. The four simulations, in one place

All standalone R, deterministic seeds. The first establishes the disease; the rest
test the cures against the failure modes.

  SIM 1 (collapse). Section 1. The old within-set test punishes a tighter,
  equally-imbalanced match (p 0.0013 -> 0). PASS = the disease is real.

  SIM 2 (cherry-pick, MODE 3). Pool with the realized match = the optimizer's
  output, reference = uniform feasible matchings. The feasible-matching percentile
  and the efficiency ratio sat at exactly 1.000 across every match tightness ---
  uninformative. This is what Codex's protocol-replay (reference = the same search)
  is built to fix, and does fix when the search is recorded.

  SIM 3 + 4 (optimism gap, MODE 1). Held-out permutation-importance gap. The gap
  detects true held-out confounding (floor-sensitive, escapes MODE 4: as b_eval
  rises 0 -> 0.40 the percentile rises 0.16 -> 1.00). But its reference variance
  equals the complete-randomization variance (var_ratio about 1) except when the
  held-out covariate is strongly predicted by the matched covariates --- the one
  regime where the held-out check is redundant. The decisive head-to-head:

```
   R rho b_eval d_real p_gap   p_cr abs_diff sd_ratio
   4 0.3    0.2 0.2536 0.982 0.9940   0.0120   1.1702
  16 0.3    0.2 0.0215 0.166 0.2075   0.0415   1.3277
   4 0.6    0.2 0.0346 0.287 0.2900   0.0030   0.9888
  16 0.6    0.2 0.4056 1.000 1.0000   0.0000   1.1841
   4 0.9    0.2 0.1264 0.990 0.8365   0.1535   0.5425
  16 0.9    0.2 0.2239 1.000 0.9940   0.0060   0.6182
```

  At rho = 0.3 (the held-out covariate genuinely informative) p_gap and p_cr agree
  to 0.01-0.04 and sd_ratio is about 1: the optimism gap IS the complete-
  randomization comparison (MODE 1). The escape (sd_ratio about 0.55) shows up only
  at rho = 0.9, where matching already balanced the held-out covariate.

  SIM 5 (PRFP floor, MODE 4). Section 5. The gate magnitude is floor-sensitive;
  the percentile is an erratic, unreliable floor indicator.

The simulations and the skeptic panels agree. Nothing escapes.

---

## 6a. The impossibility, adversarially attacked

The impossibility was not asserted and left alone. A dedicated skeptic was told to
break it by constructing a survivor --- a reference that is non-collapsing, not the
complete-randomization scale, and able to track matched-covariate confounding. It
tried five distinct constructions: within-caliper alternate-control
re-randomization; the sign-flip null; a pair bootstrap; a matched-vs-discarded-
control selection reference; and a "floor-aware" held-out-but-confounded covariate.
Every one landed in one of the three sources and hit a named mode. The hardest, the
floor-aware held-out covariate, looked like a survivor --- non-collapsing AND
tracking the confound --- until its percentile was shown to be analytically
1 - pnorm(gap / sqrt(2/S)), which is the complete-randomization z-test of that
covariate's marginal gap: MODE 1 relabeled onto a held-out covariate. The skeptic
could not break the impossibility, and said so.

A second skeptic broke an over-reach in one draft of the recommendation, and the
break is worth keeping. That draft added a design-sensitivity number
Gamma_b = sqrt(M_tol / M), a "how much unobserved imbalance would overturn the
verdict" margin. But Gamma_b is a strictly monotone transform of the magnitude M
(Spearman -1), exactly the defect that killed lens D's Gamma* --- a bare-magnitude
relabel, MODE 2. The lesson is folded into Section 7: do not report a
design-sensitivity margin of that form, because it is the forbidden magnitude
wearing a sensitivity label. The two-number report below avoids it (its Number 2 is
a genuine reference distribution, not a transform of the magnitude).

## 7. What to report instead (the constructive part)

The impossibility forbids a single non-collapsing, floor-sensitive, non-forbidden
calibration. It does not forbid answering the underlying questions. They are
separate questions, and they have separate, honest answers. Report two numbers (a
third is optional), and never let one masquerade as the other.

NUMBER 1 --- Absolute balance, as a descriptive magnitude with a pre-specified
tolerance. Report the per-covariate standardized mean differences on a FIXED
pool-scale metric (and, if a scalar is wanted, L = d' Sigma0^+ d), against a
scientific tolerance delta declared before scoring. This is floor-sensitive and
non-collapsing (Sigma0 is fixed, so L -> d0' Sigma0^+ d0 < infinity as the match
tightens), and it is what actually answers "is the residual confounding small."
State it as what it is: a descriptive magnitude, not a calibrated percentile. The
challenge calls this an already-available number and forbids dressing it up as a
test --- so do not dress it up. Report it bare, with the tolerance, and own that
its non-collapse is bought by being the magnitude rather than a reference
distribution. Trades: this is the corner-(2)-for-this-question concession --- you
do not get a percentile here, because none exists.

NUMBER 2 --- Design-search calibration (Codex's PRFP, kept intact). The percentile
of L_obs among replays of the FULL predeclared protocol, INCLUDING the
specification search, under a documented law Q_P over search paths; optionally
SER. This is non-collapsing (the reference is over analogous matched designs of
the same labeled pool, whose spread does not shrink to zero), design-respecting,
and --- given a recorded search --- not gameable by cherry-picking. It is
floor-blind by construction, so read it only as the answer to its own question:
"did this design extract the balance the protocol could, and is its favorable
balance robust to the design search I actually ran." Run it AFTER Number 1 passes
its gate, exactly as Codex specifies, so a floor-blind percentile can never
certify a design that the magnitude already failed. Trades: sacrifices corner (1)
(no within-set treatment re-randomization) and corner (3) toward analogous designs
(the protocol-feasible matchings) --- both permitted.

NUMBER 3 (optional) --- Held-out transportability diagnostic. The held-out
optimism gap of lens B, reported as a DIAGNOSTIC, never as a calibration, with the
eval/match split PRE-REGISTERED before seeing data and the eval covariates not
selectable. State the caveat the simulations force: when the held-out covariate is
informative, this diagnostic's reference is the complete-randomization comparison
on that covariate, so it is the forbidden comparison restricted to held-out
covariates --- useful as a flag for over-fitting and poor generalization, not as a
certificate of balance.

Why this is not l* and not complete randomization. Number 2 randomizes the
matched-design SEARCH PATH, holding treatment labels fixed; it never permutes Z
over the pool and never merges matched sets or sweeps a coarsening ladder, so it
is neither the complete-randomization comparison nor l*. Number 1 is a fixed
magnitude with a fixed tolerance, not a re-randomization at any resolution.
Number 3 is labeled, honestly, as the forbidden comparison on held-out
covariates, which is why it is demoted to a diagnostic and not offered as the
answer.

---

## 8. Surviving objections I cannot fully answer

  - The descriptive magnitude (Number 1) does not satisfy the challenge's goal 3
    (a percentile or efficiency ratio a collaborator can act on). My answer is the
    impossibility: for the matched-covariate balance question there is no
    non-forbidden percentile, so the magnitude with a pre-specified tolerance is
    the most a researcher can honestly report. A reader who insists on a percentile
    for THAT question is asking for something Section 2 shows cannot exist.

  - PRFP (Number 2) is only as good as the recorded search and the sampler Q_P. An
    analyst who explored specifications without logging them cannot replay them,
    and a poorly mixed Q_P gives a misleading percentile. Codex's own hardest
    objection. The honest fallback is a sensitivity analysis over protocol scope
    (final recipe only; plus caliper search; plus distance and discard search),
    which is weaker than a single calibrated number.

  - The whole construction assumes the balance question is about OBSERVED
    covariates. None of it touches hidden bias. A design that passes Number 1's
    gate and scores well on Number 2 can still be confounded on an unmeasured
    covariate; that is a separate analysis (design sensitivity for the treatment
    effect, not for balance).

  - The strengthened impossibility (Section 2) is an argument, not a published
    theorem. It enumerates three sources of a bounded reference spread and shows
    each is forbidden or collapsing; I have not proved the enumeration is
    exhaustive in full generality, only that every concrete proposal across two
    engines and six lenses, plus five further constructions a dedicated skeptic
    tried (Section 6a), fell into one of the three. Two specific seams remain. (a)
    The clean cases cover statistics that are LINEAR in X within sets (the
    mean-difference family); the dichotomy lemma is airtight for any held-out
    linear-in-E statistic, but a genuinely nonlinear interaction statistic
    evaluated in-sample is covered only by appeal to Fact B (the collapse is
    statistic-agnostic), not by an independent derivation. A determined inventor
    should be aimed there. (b) The variance identity "held-out reference variance =
    2*sigma^2/S" is exact only for iid pool-scale noise; for finite re-pairing
    references it carries a finite-population correction, so the load-bearing part
    of the lemma is the MODE-4 blindness, which is correction-free. The honest
    status is a strong, simulation-backed conjecture with no counterexample yet,
    and an open invitation to find a fourth source.

---

## 9. Which lens won (one line, for cross-run comparison)

graft (A feasible-matchings + E efficiency-ratio + an absolute fixed-scale
magnitude gate), best realized in Codex's PRFP; but no lens survives as a
standalone non-collapsing, floor-sensitive, non-forbidden calibration --- the
result is a strengthened impossibility plus a two-number report (descriptive
magnitude + gate; protocol-replay design-search percentile).

claude
