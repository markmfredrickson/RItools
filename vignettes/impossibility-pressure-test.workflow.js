export const meta = {
  name: 'impossibility-pressure-test',
  description: 'Adversarially search for a design-internal single-number balance calibration that does not collapse; compare external (whole-pool CRE) references',
  phases: [
    { title: 'Generate', detail: 'inventors propose design-internal non-collapsing calibrations; each runs the harness' },
    { title: 'Refute',   detail: 'breakers try to break each self-claimed escape with their own sims' },
    { title: 'Compare',  detail: 'evaluate the external whole-pool-CRE references, incl. oracle-free covariance adjustment' },
    { title: 'Synthesize', detail: 'adjudicate the impossibility, rank references, stats + writing critics, write memo' },
  ],
}

// Shared context every agent needs. The harness is self-contained base R (no
// package load). cwd is the repo root.
const HARNESS = `
SHARED HARNESS (already on disk): vignettes/pressure-test-helpers.R -- self-contained base R, NO devtools/package needed.
It defines: make_matched(shrink, g, seed) [hierarchical matched data, FIXED per-covariate gap g, match tightness = shrink, returns list with X, z, set, group, pool, sd_pool];
low-level adj_diff_sum/Vperm/d2_stat/eig_rank/perm_within/perm_pool/smd_vec/pooled_sd/ref_percentile;
collapse_curve(calib) [gap FIXED, shrink 0.6->0.012; returns $curve(shrink,alarm), $ref, $collapsed];
magnitude_curve(calib) [match FIXED, gap 0->1; alarm should RISE with gap];
named contenders cal_d2_internal, cal_maxsmd_internal (both DESIGN-INTERNAL, both COLLAPSE),
cal_maxsmd_poolCRE, cal_maha_poolCRE (EXTERNAL, do not collapse but FLOOR-BLIND),
cal_maxsmd_poolCRE_adj (EXTERNAL + covariance adjustment on ORACLE group; non-collapsing AND magnitude-sensitive),
cal_maxsmd_bare (bare magnitude, ref="none").

A calib is: function(dat) -> list(alarm = <numeric, higher = design looks more imbalanced>, ref = "internal"|"external"|"none").
Run R as: R_LIBS=.local Rscript -e '<code>'  OR  R_LIBS=.local Rscript path/to/file.R . Always source("vignettes/pressure-test-helpers.R") first.

DEFINITIONS (be strict):
- DESIGN-INTERNAL: the reference/null that defines "how alarmed" is generated from the REALIZED strata (within-set re-randomization, or any function of the within-stratum randomization covariance V_d). NOT a fixed coarser design, NOT a pooled pre-match scale, NOT a fixed substantive threshold.
- COLLAPSE: with the imbalance held FIXED (gap g constant) and the match tightening (shrink -> 0, V_d -> 0), alarm climbs toward its maximum. A beautifully matched design is reported as ever-more-imbalanced for a gap that never changed.
- A genuine escape must be ALL of: (i) single number over all covariates, (ii) ref="internal" honestly, (iii) collapse_curve NOT collapsing, (iv) magnitude_curve RISING with the gap (not floor-blind).
DIMENSIONAL ARGUMENT to confirm or break: any design-internal calibration is a function of (dbar, V_d). A scale-free (degree-zero) function of (dbar, V_d) depends only on dbar/sqrt(V_d), which -> infinity as V_d->0 (collapse). The only functions bounded as V_d->0 ignore V_d -- i.e. they are bare magnitudes (Ruler A), not calibrations. The claim: no single scalar can be design-internal, magnitude-sensitive, AND non-collapsing at once.
`;

const INVENTOR_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['idx','name','ref_claim','escapes','argument','collapse_alarm','magnitude_alarm'],
  properties: {
    idx: { type: 'integer' },
    name: { type: 'string' },
    strategy: { type: 'string' },
    ref_claim: { type: 'string', enum: ['internal','external','none'] },
    escapes: { type: 'boolean', description: 'true only if design-internal AND non-collapsing AND magnitude-sensitive' },
    collapse_alarm: { type: 'array', items: { type: 'number' }, description: 'alarm across shrink 0.6->0.012 at fixed gap' },
    magnitude_alarm: { type: 'array', items: { type: 'number' }, description: 'alarm across gap 0->1 at fixed match' },
    file: { type: 'string', description: 'path to the saved candidate .R file' },
    argument: { type: 'string', description: 'honest verdict: why it escapes, or why it collapsed/was secretly external' },
  },
};

const BREAKER_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['candidate_idx','verdict','failure_mode','evidence'],
  properties: {
    candidate_idx: { type: 'integer' },
    verdict: { type: 'string', enum: ['survives','broken'] },
    failure_mode: { type: 'string', enum: ['collapses_on_harder_dgp','secretly_external_or_bare','floor_blind','invalid_calibration','none','other'] },
    evidence: { type: 'string', description: 'concrete sim numbers supporting the verdict' },
  },
};

const COMPARE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['reference','collapses','magnitude_sensitive','valid_calibration','recommendation'],
  properties: {
    reference: { type: 'string' },
    collapses: { type: 'boolean' },
    magnitude_sensitive: { type: 'boolean' },
    floor_blind: { type: 'boolean' },
    valid_calibration: { type: 'string', enum: ['yes','no','unknown'], description: 'is alarm ~Uniform under the true reference assignment mechanism' },
    oracle_free_works: { type: 'string', enum: ['yes','no','na','unknown'], description: 'does it still work when covariance adjustment uses observed covariates (regression), not oracle group' },
    unbuilt_design_cost: { type: 'boolean', description: 'does it compare to a design that was not actually built' },
    numbers: { type: 'string', description: 'key simulated numbers' },
    recommendation: { type: 'string' },
  },
};

// ---------------------------------------------------------------- Phase 1
phase('Generate')
const STRATEGIES = [
  { idx: 1, hint: 'A RATIO of two design-internal quantities chosen so V_d cancels -- e.g. observed d^2 divided by its own within-set null mean/quantile. Check whether canceling V_d just reproduces the degree-zero collapse.' },
  { idx: 2, hint: 'The generalized eigenvalues of within-stratum covariance W versus total covariance Tot (the fraction of variance in each direction surviving stratification). Build ONE number from that spectrum that you claim tracks imbalance magnitude.' },
  { idx: 3, hint: 'A SELF-NORMALIZED / studentized omnibus: divide the gap by a robust internal spread estimate engineered NOT to vanish as the match tightens.' },
  { idx: 4, hint: 'A Bayesian / shrinkage posterior probability that the standardized within-set imbalance exceeds a FIXED substantive threshold (note: a fixed threshold may make it external -- judge honestly).' },
  { idx: 5, hint: 'An e-value / betting martingale accumulated across the matched sets.' },
  { idx: 6, hint: 'A rank- or sign-based (Wilcoxon-like) design-internal statistic at the matched sets, calibrated by within-set re-randomization.' },
  { idx: 7, hint: 'An R^2 / information-criterion of treatment-on-covariates within strata, corrected for the K/(n-S) chance floor.' },
  { idx: 8, hint: 'An explicit tightness PENALTY: take the collapsing alarm and subtract a function of the match tightness, trying to cancel the collapse. Check whether the penalty needs information you only get from a fixed external reference.' },
  { idx: 9, hint: 'WILDCARD: any other design-internal single-number calibration you believe escapes. Be creative and genuinely adversarial toward the impossibility claim.' },
];
const inventors = (await parallel(STRATEGIES.map(s => () =>
  agent(
`${HARNESS}

You are an INVENTOR. Goal: find a DESIGN-INTERNAL, single-number balance calibration that does NOT collapse as the match tightens, yet still RISES with real imbalance. If you succeed you refute a standing impossibility claim, so try hard and be honest.

Your assigned strategy (idx ${s.idx}): ${s.hint}

Steps:
1. Implement your candidate as an R function cal_candidate(dat) returning list(alarm=<numeric, higher=worse>, ref="internal"|"external"|"none"). Be ruthlessly honest about ref: if your "reference" is a fixed coarser design, a pooled pre-match scale, or a fixed threshold, it is "external", NOT "internal".
2. Save it to vignettes/cand-${s.idx}.R (source the harness at the top of that file).
3. Run BOTH: collapse_curve(cal_candidate) and magnitude_curve(cal_candidate). Print the alarm vectors.
4. Decide escapes = TRUE only if ref="internal" AND collapse_curve did NOT collapse AND magnitude_curve rises with the gap.

Return the structured result. In argument, state plainly what happened and, if it collapsed or turned out external, WHY (connect to the dimensional argument).`,
    { label: `invent:${s.idx}`, phase: 'Generate', schema: INVENTOR_SCHEMA }
  )
))).filter(Boolean);

log(`inventors returned ${inventors.length}; self-claimed escapes: ${inventors.filter(i => i.escapes).length}`)

// ---------------------------------------------------------------- Phase 2
phase('Refute')
const claimed = inventors.filter(i => i && i.escapes);
let breakerResults = [];
if (claimed.length === 0) {
  log('no candidate self-claimed an escape; impossibility holds at the generate step')
} else {
  breakerResults = await pipeline(
    claimed,
    cand => parallel(['HARDER-DGP: re-run collapse_curve with smaller gap (g=0.05), more covariates, and shrink down to 0.005; does alarm still stay bounded?',
                      'AUDIT-REFERENCE: prove whether the reference is truly internal (must move when ONLY the match tightness changes at fixed gap) or secretly external/bare; mislabeled ref="internal" => broken.',
                      'FLOOR-BLIND + VALIDITY: re-run magnitude_curve and check alarm is ~Uniform at gap=0 under the true within-set null; flat alarm or invalid calibration => broken.']
      .map(lens => () =>
        agent(
`${HARNESS}

You are a BREAKER. Default assumption: the candidate is BROKEN unless it clearly survives your lens. Source the candidate: source("vignettes/cand-${cand.idx}.R") (it sources the harness). Its function is cal_candidate(dat).

Candidate idx ${cand.idx} ("${cand.name}") claims to be a design-internal, non-collapsing, magnitude-sensitive calibration. Its inventor reported collapse_alarm=[${(cand.collapse_alarm||[]).join(', ')}] and magnitude_alarm=[${(cand.magnitude_alarm||[]).join(', ')}].

YOUR LENS: ${lens}

Run the relevant R yourself, get numbers, and return a verdict. Be specific and adversarial.`,
          { label: `break:${cand.idx}:${lens.slice(0,12)}`, phase: 'Refute', schema: BREAKER_SCHEMA }
        )
      )
    )
  );
  breakerResults = breakerResults.filter(Boolean).flat().filter(Boolean);
}

// a candidate SURVIVES only if no breaker broke it
const survivors = claimed.filter(c =>
  !breakerResults.some(b => b && b.candidate_idx === c.idx && b.verdict === 'broken')
);
log(`survivors after refutation: ${survivors.length}`)

// ---------------------------------------------------------------- Phase 3
phase('Compare')
const REFERENCES = [
  { name: 'whole-pool CRE, max|SMD| (cal_maxsmd_poolCRE)', focus: 'Confirm non-collapse and quantify the FLOOR-BLINDNESS: how large a residual gap before alarm is informative? Is it a usable test or always "passes"?' },
  { name: 'whole-pool CRE, fixed-Sigma Mahalanobis (cal_maha_poolCRE)', focus: 'Why is it even more floor-blind than max|SMD|? Is the fixed-Sigma multivariate metric salvageable, or dominated by between-group structure?' },
  { name: 'whole-pool CRE + covariance adjustment, ORACLE group (cal_maxsmd_poolCRE_adj)', focus: 'Confirm it both does-not-collapse and is magnitude-sensitive. Is the alarm a VALID calibration (uniform under the true CRE at gap=0)? State the unbuilt-design cost plainly.' },
  { name: 'whole-pool CRE + covariance adjustment, ORACLE-FREE (regression on observed covariates)', focus: 'THE KEY TEST. Re-implement adjustment by residualizing X on a model of OBSERVED covariates / an estimated propensity or prognostic score (NO oracle group). Does it still escape collapse AND stay magnitude-sensitive AND valid? Implement and run it.' },
  { name: 'fixed substantive tolerance on max|SMD| (Ruler A baseline)', focus: 'The non-calibrated baseline: a fixed cutoff c on max|SMD|. Non-collapsing by construction; document its cost (no built-in reference; per-covariate; you must choose c) and contrast with the adjusted-CRE.' },
];
const comparisons = (await parallel(REFERENCES.map(r => () =>
  agent(
`${HARNESS}

You are an EVALUATOR comparing EXTERNAL fixed-reference balance calibrations (the escape from the design-internal collapse). Evaluate this one rigorously with your own R runs:

REFERENCE: ${r.name}
FOCUS: ${r.focus}

For your reference, report: does it collapse (collapse_curve)? is it magnitude-sensitive (magnitude_curve)? is it floor-blind? is it a VALID calibration (alarm ~Uniform under its own reference mechanism at gap=0 -- test this by simulating the reference assignments and checking the alarm distribution)? does it have the unbuilt-design cost? For the ORACLE-FREE row you MUST implement regression-based adjustment on observed covariates and run it. Return the structured result with concrete numbers.`,
    { label: `compare:${r.name.slice(0,18)}`, phase: 'Compare', schema: COMPARE_SCHEMA }
  )
))).filter(Boolean);

// ---------------------------------------------------------------- Phase 4
phase('Synthesize')
const ctx = JSON.stringify({
  inventors: inventors.map(i => ({ idx:i.idx, name:i.name, ref_claim:i.ref_claim, escapes:i.escapes,
                                   collapse_alarm:i.collapse_alarm, magnitude_alarm:i.magnitude_alarm, argument:i.argument })),
  breakers: breakerResults,
  survivors: survivors.map(s => s.idx),
  comparisons,
}, null, 1);

const memo = await agent(
`${HARNESS}

Write a memo (GitHub-flavored markdown, ASCII ONLY -- no unicode, use --- for em dash, -> for arrows, straight quotes) titled "Pressure-testing the balance-calibration impossibility". Audience: Jake Bowers (applied statistician; randomization inference) and Ben Hansen. Jake's writing rules are strict: plain words over jargon; motivate before method; name the actor, criterion, and rejected alternative for any evaluative claim; NO decorative metaphors (no machinery/scaffolding/firewall/load-bearing); make claims directly (no "it is important to note"); say what is a THEOREM vs a CONJECTURE.

The result data (inventors, breakers, survivors, external-reference comparisons):
${ctx}

Structure:
1. The question, in one paragraph: is there a design-internal, single-number balance calibration that does not collapse as the match tightens? Motivate with the matched-design-gets-an-alarming-p problem.
2. The dimensional argument (why design-internal calibrations must collapse): a design-internal calibration is a function of (dbar, V_d); scale-free combinations depend on dbar/sqrt(V_d) and blow up as V_d->0; the only V_d-bounded ones ignore V_d and are bare magnitudes. Present as the conjecture the search tested.
3. The adversarial search: how many inventors, what strategies, how many claimed escape, how many SURVIVED refutation. Report honestly. If zero survived, that is strong support (not proof) for the impossibility; name the closest near-miss and exactly how it broke.
4. The escape is external references -- compare Jake's ideas on one table: collapse? magnitude-sensitive? floor-blind? valid? unbuilt-design cost? oracle-free? Lead to the finding that raw whole-pool CRE is non-collapsing but floor-blind, and covariance adjustment is what makes it both stable and informative -- with the oracle-free caveat resolved or flagged.
5. Recommendation and honest scope: the three-way bind (single calibrated omnibus number / design-internal reference / non-collapsing -- pick two); what to actually report; that this is all about OBSERVED covariates and is not ignorability.

Return ONLY the memo markdown.`,
  { label: 'synthesize', phase: 'Synthesize' }
);

const statsCritic = await agent(
`${HARNESS}

You are a STATISTICS CRITIC. Re-derive and, where cheap, RE-RUN with your own R the key quantitative claims in this memo. Flag any claim that is wrong, overstated, or unsupported by the harness behavior (e.g., a "does not collapse" that actually collapses on a harder DGP, a "valid calibration" not checked for uniformity, a theorem/conjecture mislabel). Be concrete.

MEMO:
${memo}

Return a bullet list of required corrections (or "no corrections" if genuinely clean), each with the fix.`,
  { label: 'stats-critic', phase: 'Synthesize' }
);

const writingCritic = await agent(
`You are a WRITING CRITIC enforcing Jake Bowers's rules (plain words; motivate before method; name actor/criterion/rejected-alternative for evaluative words like "appropriate/robust/valid"; NO decorative structural/industrial/security metaphors; no "it is important to note"; no nominalizations hiding the actor; ASCII only; theorem-vs-conjecture honesty). Quote each offending phrase and give the rewrite.

MEMO:
${memo}

Return a bullet list of specific edits (quote -> rewrite), or "no edits" if genuinely clean.`,
  { label: 'writing-critic', phase: 'Synthesize' }
);

const finalMemo = await agent(
`Revise the memo below by APPLYING the stats corrections and writing edits. Keep all correct content; do not introduce new claims. ASCII only. Return ONLY the final memo markdown.

STATS CORRECTIONS:
${statsCritic}

WRITING EDITS:
${writingCritic}

MEMO:
${memo}`,
  { label: 'finalize', phase: 'Synthesize' }
);

return {
  impossibility_supported: survivors.length === 0,
  inventors_total: inventors.length,
  claimed_escapes: claimed.length,
  survivors: survivors.map(s => ({ idx: s.idx, name: s.name })),
  comparisons,
  memo: finalMemo,
};
