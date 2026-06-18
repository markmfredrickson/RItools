# The magnitude (SMD) screen: a non-collapsing balance measure

Date: 2026-06-17 (rewritten). Design memo for RItools. Companion to
sim-results-memo.qmd (B7 = the SMD does not collapse; B8 = the small-sample floor).

## The point: the SMD prevents the collapse

The within-set omnibus (d^2, and ACAT) is homogeneous of degree zero in the
imbalance: it reads the systematic DIRECTION of imbalance relative to within-set
re-randomization and is blind to its SIZE. So a tight match with substantively
tiny imbalance still gets a small omnibus p ("imbalanced!"), and a large matched
sample rejects on trivial bias. That is the collapse (A1, A3 in the qmd).

The standardized mean difference

    SMD = std.diff = adj.diff / pooled.sd,  pooled.sd = sqrt((var_treated + var_control)/2)

does NOT collapse, because its denominator is a FIXED scale (the covariate's
pooled spread, computed once from the data and dominated by between-stratum
variation), not the within-set null variance that shrinks as the match tightens.
So the SMD tracks the real size of imbalance: as a match tightens the SMD falls
(B7, where SMD(x1) goes 0.87 -> 0.05 while the omnibus p stays ~0).

Therefore: judging balance by the SMD is how we prevent the collapse. We do not
"fix" the omnibus p -- the impossibility result forbids a non-collapsing
within-set test. We stop leaning on the collapsing p and report a magnitude that
does not collapse instead. This is the constructive form of "lead with
magnitude," and it is exactly what a Love plot already does.

Two facts make this cheap:
- balanceTest() ALREADY computes and reports std.diff (R/Design.R:652), with
  adj.diff and pooled.sd. So the statistic exists; this is a reporting/decision
  layer, not a new estimator.
- covariate.scales (R/Design.R:646) lets a user fix the denominator (e.g. to
  pre-matching SDs), which guarantees the scale never shifts.

The SMD comes in two forms, and we want BOTH:
- per-covariate SMDs (the Love plot), and
- a single GLOBAL summary, max|SMD| across covariates -- the one-number
  "how balanced is the worst covariate, in SD units" measure (used in B8).

## Three things called "screen" -- keep them distinct

This bundling caused confusion; only one prevents the collapse.

1. DESIGN screen (w_j / singular-value): drops covariates or directions with
   ~zero WITHIN-STRATUM variance. Purpose: numerical stability + honest
   abstention (decisions 1, 2 elsewhere). It does NOT prevent the collapse (B3:
   surviving-covariate p is still degree-zero). Numerical hygiene, not the fix.
2. MAGNITUDE criterion (the SMD): the subject of this memo. DOES prevent the
   collapse. Per-covariate and/or global max|SMD|.
3. Pre-filter question (decision D below): a narrow technical caveat about whether
   the magnitude criterion may choose the omnibus's inputs. Not the main point.

"Dropping covariates with collapsing within-stratum variance" is (1), the design
screen -- a separate, numerical thing. The collapse fix is (2).

## What to build

- A. CUTOFF (APPROVED). Default margin 0.25 SD, changeable via an argument
  (e.g. balance.margin = 0.25); a covariate is marked when |std.diff| >= margin.
  Caveats on the value and denominator: see Provenance below -- the default value
  (0.25 vs 0.1) and the denominator (pooled vs treated-group SD) both deserve an
  explicit, documented choice; keep the margin settable.

- B. DISPLAY (APPROVED, expanded). Surface the magnitude in every place a user
  sees a result:
  1. a per-covariate marking (a column in the printed table / results object for
     |std.diff| >= margin);
  2. a +/- margin reference band on the Love plot (plot method);
  3. a one-line summary (e.g. "3 of 12 covariates exceed the 0.25-SD margin");
  4. a GLOBAL max|SMD| balance number (the worst standardized imbalance).

- C. SCOPE NOW (APPROVED). Build the marking (A + B) plus small-sample handling
  (next section). Defer a formal equivalence test (TOST / Hartman & Hidalgo 2018).

- D. SEPARATION (OPEN -- still to decide). Keep the magnitude criterion as
  something the user READS, with NO effect on which covariates the omnibus is
  computed on. The alternative (a "pre-filter": run the omnibus only on covariates
  with large SMD) is invalid as a test -- selecting covariates because their
  observed imbalance is large and then testing those same covariates rejects far
  too often under the null (selection on the outcome being tested). Contrast: the
  design screen (1) MAY filter the omnibus, because it decides using the
  within-stratum variance, which never looks at the observed treated-control
  difference. The magnitude criterion does look at it, so it must not filter.
  EVIDENCE (B9, added 2026-06-17): on NULL data (200 pairs, k=10/20, no degeneracy),
  selecting covariates with |SMD| >= tau and testing only those makes the omnibus
  reject 0.98-1.00 of the time CONDITIONAL on selection (vs nominal 0.05) once
  tau >= 0.15, and inflates the OVERALL size to ~0.12-0.19 at tau=0.10-0.15. The
  0.25 case looks "safe" in overall rate (~0.001) only because selection rarely
  fires at this n -- its conditional size is still 1.0. So pre-filtering is invalid.
  Recommendation: yes, keep separate. (Still awaiting Jake's explicit call.)

## The small-sample problem (REMEMBER) and the planned permutation approach

This is the piece to carry forward; we want a permutation-based approach here.

The problem (B8). The SMD is a point estimate with sampling noise. Under TRUE
balance (treatment as-if-random within strata) the observed SMD is not zero; its
null spread is about

    rho_match * sqrt(2 ln k) / sqrt(ESS),   rho_match = rms within-stratum diff / pooled SD,

for the family-wise max over k covariates. So a fixed 0.25 cutoff can flag a
perfectly balanced design by chance when ESS is small and/or k is large (B8:
ESS=25 flags a balanced design ~42% of the time at 0.25; ESS=100 ~0%). The
threshold should be the LARGER of the substantive margin (0.25, an effect-size
judgment, sample-size-independent) and this statistical floor.

Planned permutation approach (FUTURE -- Jake wants this). Instead of the analytic
floor above, compute the SMD's small-sample null distribution EXACTLY for the data
at hand: permute the treatment assignment WITHIN strata many times, recompute the
SMD (and max|SMD|) each time, and read off a data-specific NOISE BAND (e.g. the
0.95 quantile of |SMD| under balance). This needs no normal approximation and
respects the actual stratum structure and covariate scales.

CRITICAL CAVEAT -- do not turn the permutation band into a test, or the collapse
comes back. The within-stratum permutation null of the SMD SHRINKS as the match
tightens (its numerator's null spread shrinks with the within-stratum scatter,
while the fixed pooled-SD denominator does not). So if you FLAG whenever the
observed SMD exceeds its permutation band, you re-create the collapse exactly: a
tight match has a tiny band, so even a substantively negligible SMD lands outside
it. The permutation band is therefore an UNCERTAINTY display, not a reject rule.

The honest use is a TWO-NUMBER report, kept distinct:
- MAGNITUDE: is |std.diff| (or max|SMD|) above the substantive margin (0.25)?
  This does not collapse and answers "is the imbalance big enough to matter?"
- PRECISION: where does the observed SMD sit relative to its within-stratum
  permutation band? This answers "is the point estimate even distinguishable from
  sampling noise at this ESS?" -- a caveat on how much to trust the point estimate,
  NOT a balance verdict.
At small ESS the band is wide, telling the user not to over-read a noisy SMD; at a
tight match the band is narrow but the SMD is also small, so magnitude (not the
band) carries the verdict. (A complete-randomization / pooled permutation band
would NOT shrink with the match, but it ignores the strata -- the same
trilemma corner Rosenbaum's marginal benchmark sits in -- so it is not the
design-respecting choice. Use the within-stratum band, as uncertainty only.)

This two-number framing matches the project's standing recommendation
(descriptive magnitude + design percentile, kept separate).

## Provenance of the 0.25 rule (summary; full version in sim-results-memo.qmd)

- "0.25, Rosenbaum & Rubin (1985)" is a MISATTRIBUTION: that 0.25 is a
  propensity-score CALIPER, not a balance cutoff.
- The 0.25 BALANCE threshold is Rubin (2001) as restated by Stuart (2010);
  Rubin's own words are "half a standard deviation," for the propensity score.
- Two camps on the NUMBER: matching tradition 0.25 (Stuart/MatchIt); medical
  tradition 0.1 (Normand 2001; Austin 2009/2011).
- Two camps on the DENOMINATOR: treated-group SD (Stuart/MatchIt) vs pooled SD
  (Austin/medical). RItools uses pooled SD -- the lineage that pairs with 0.1 --
  so 0.25 + pooled SD mixes conventions.
- Ho, Imai, King & Stuart (2007): balance is in-sample and p-values should not be
  used as balance criteria -- prior art for leading with magnitude over the p.
- Confidence: Rubin 2001, Stuart 2010, Austin 2009/2011 quotes are verbatim; R&R
  1985 caliper and Normand 2001 are from secondary quotation -- confirm those two
  primaries before citing in a paper.

## Open decisions / TODO

- D (separation): keep the magnitude screen out of the omnibus computation? (Rec:
  yes.) Awaiting Jake.
- Default value + denominator: 0.25 + pooled SD (mixed), 0.1 + pooled SD
  (consistent medical), or 0.25 + treated-group SD (consistent matching)? Keep
  settable; pick a documented default.
- PERMUTATION noise band (future): implement the within-stratum permutation null
  of the SMD / max|SMD| as an uncertainty display -- NOT a reject threshold (see
  the critical caveat). This is the small-sample piece to remember.
- Formal equivalence test (TOST / Hartman & Hidalgo 2018): later.
