# Challenge prompt: invent a non-collapsing, design-respecting balance calibration

Hand this prompt to a fresh Claude or Codex session to launch a multi-agent
method-invention search. It is self-contained: the receiving AI does not need this
repository, the HANDOFF, or the memory store.

---

ROLE

You are a mathematical statistician working on randomization-based inference for
covariate balance in matched observational studies. Your job in this session is to
INVENT a new method, not to summarize the literature. Be adversarial with your own
ideas: state assumptions, give the strongest objection to each proposal, and discard
what does not survive. Plain ASCII only (use --- for em dash, -- for en dash, -> for
arrows, straight quotes). No decorative metaphors, no "is appropriate / robust /
natural" without a named criterion. For a hypothesis TEST there is no "estimand" --
do not use that word.

BACKGROUND (self-contained; you do not have the repository)

A matched design pairs or groups treated and control units so that, within each matched
set s, the units are similar on observed covariates X. To assess balance, the
Hansen-Bowers (2008) d^2 omnibus test computes, per covariate, the within-set adjusted
treated-minus-control mean difference, stacks these into a vector d, and forms

    T = d' V_d^+ d,

where V_d is the covariance of d UNDER THE WITHIN-SET RANDOMIZATION NULL: the
distribution obtained by re-randomizing the treatment indicator within the realized
matched sets, holding the sets fixed. V_d^+ is the Moore-Penrose inverse. Under that
null T is approximately chi-square, and the reported p-value is the upper-tail
probability, equivalently the percentile of the observed balance among within-set
re-randomizations of THIS design. Call this the "within-set standard."

THE PROBLEM TO SOLVE: THE COLLAPSE

V_d scales with the within-set spread of X. The whole point of matching is to shrink
that spread. So as a match tightens:

    within-set var(X) -> 0   ==>   V_d -> 0   ==>   for fixed residual d,  T -> infinity,
    p -> 0.

A BETTER match (smaller within-set spread, same residual mean imbalance) is reported as
WORSE balance. The test rewards loose matches and punishes tight ones. This is the
central pathology. It is not a bug in the code; it is what "surprising relative to
within-set re-randomization" means when the within-set reference has almost no spread
left.

Two facts already established, which constrain your search:

  (A) METRIC RESCALING DOES NOT HELP. Replacing V_d with any FIXED covariance Sigma (not
      recomputed under the null) cancels in the within-block percentile: the rank of the
      observed statistic among within-set re-randomizations is invariant to a fixed
      positive-definite metric. So you CANNOT fix the collapse by choosing a cleverer
      fixed Mahalanobis metric. The reference DISTRIBUTION, not the metric, is what must
      change.

  (B) THE COLLAPSE IS STATISTIC-AGNOSTIC. It is not special to d^2. The F test, the
      logistic-regression likelihood-ratio test, and the Westfall-Young max-T all
      collapse under the same within-set standard for the same reason. So a different
      test statistic, by itself, will not save you.

THE SUBSTANTIVE GOAL (what a good method must deliver)

A calibrated, interpretable assessment of whether a realized matched design is
well-balanced that:
  1. does NOT punish tighter matching (monotone-sane: improving balance should not
     worsen the reported quality);
  2. respects the matched structure (it is an assessment of THIS design, not of an
     unrelated design);
  3. yields a number a researcher can act on with collaborators -- ideally a percentile
     or efficiency ratio of the form "this design is better balanced than X% of
     [reference class]" or "this design achieves a fraction X of the achievable
     balance," not merely a raw magnitude.

THE TRILEMMA (the formal obstruction)

Against the within-set standard, ANY calibration collapses as the match tightens. So one
cannot simultaneously have all three of:

    (1) the within-set standard (re-randomize treatment within the realized sets);
    (2) a non-collapsing calibration;
    (3) no comparison to designs the analyst did not actually build.

You must sacrifice at least one corner. The ONE non-negotiable is corner (2): your method
must keep a NON-COLLAPSING CALIBRATION. You may sacrifice corner (3) (compare to analogous
matched designs you did not build) OR corner (1) (abandon within-set treatment
re-randomization in favor of a different reference -- e.g., sample-splitting, held-out
covariates, or a design-sensitivity-style analysis), or both. What you may NOT do is the
two specific things fenced off under HARD CONSTRAINTS below: comparison to NON-analogous
designs (complete randomization), and the l* resolution sweep. State plainly which
corner(s) you trade and why.

HARD CONSTRAINTS ON THE SOLUTION

  - IF you relax corner (3), do it ONLY toward ANALOGOUS designs. The comparison class
    must be other MATCHED designs constructed from the same data (the same treated/control
    pool, the same covariates), or a penalized/anti-overfit version of the realized
    design. Examples of permitted reference classes: "all feasible matchings of this
    pool," "the most-balanced feasible matching consistent with the data," "matchings
    under permuted or noise covariates," "held-out-covariate balance of this matching."

  - IF you relax corner (1) instead (or in addition), permitted moves include
    sample-splitting (choose the match on one part of the data, assess balance on a
    held-out part), held-out-covariate balance (match on some covariates, evaluate on
    covariates withheld from the matching), and design-sensitivity-style analyses (report
    how large an unobserved within-set imbalance, or how much overfitting optimism, would
    have to be present to overturn the balance verdict). These abandon within-set
    treatment re-randomization, which is allowed -- but you must then say what your new
    reference is and why it does not collapse.

  - FORBIDDEN: comparison to a design that is not analogous to the matched design --
    above all, comparison to complete randomization of the whole pool. That comparison
    answers the matching-as-pruning question, which is rejected here. Do not propose it,
    and do not propose anything that reduces to it.

  - FORBIDDEN: reinventing the "resolution profile" / l*. For your awareness so you do
    not rediscover it: that idea sweeps a LADDER of coarsenings of the matched sets (from
    the realized sets up toward complete randomization) and reports the coarsening level
    l* at which the within-set percentile crosses 0.5. It keeps corner (1) and trades on
    a within-set re-randomization at varying resolution. Your method must NOT be a
    resolution/coarsening sweep and must not report an l*-like crossing point. If your
    idea starts to look like "re-randomize within progressively merged sets," stop and
    change direction.

  - FORBIDDEN: reporting a raw fixed-metric magnitude and calling it a test. A magnitude
    like d' cov(X_pool)^{-1} d is already on the table as a descriptive number; it is NOT
    a calibration and does not satisfy goal 3. You must produce an actual reference
    distribution or a principled normalization, not a bare distance.

PERMITTED AND ENCOURAGED DIRECTIONS (suggestive, not prescriptive -- invent your own)

  - A reference distribution over FEASIBLE MATCHINGS: define a distribution over valid
    matchings of the same pool (e.g., randomize which control attaches to which treated
    subject within caliper/feasibility constraints), compute balance for each, and place
    the realized matching's balance in that distribution. Note this re-randomizes the
    MATCHING ASSIGNMENT, not the treatment label -- a different null from corner (1).

  - ANTI-OVERFITTING / OPTIMISM CORRECTION: matching is an optimization that SELECTS a
    balanced configuration, so good realized balance is partly selection optimism.
    Quantify and subtract that optimism. Candidate moves: match on PERMUTED or NOISE
    covariates to get the balance achievable by chance optimization (a permutation-
    importance analogue for matching); an effective-degrees-of-freedom / complexity
    penalty for the matching procedure; sample-splitting (choose the match on one part,
    assess on another); CROSS-VALIDATED or HELD-OUT-COVARIATE balance (match on some
    covariates, evaluate balance on covariates withheld from the matching). Out-of-sample
    balance does not collapse the way in-sample within-set balance does, and it detects
    overfitting directly.

  - EFFICIENCY-RATIO normalization against the achievable optimum: normalize realized
    balance by the best feasible matching (lower bound) and a null/expected matching
    (upper bound) to get a bounded "fraction of achievable balance" that improves as the
    match tightens.

  Combine these if it helps. The strongest answer may pair a feasible-matching reference
  class WITH an overfitting correction so the percentile cannot be gamed by cherry-picking
  the single most-balanced matching.

PRIOR ART YOU MUST POSITION AGAINST

  - Branson (2021), "Randomization Tests to Assess Covariate Balance When Designing and
    Analyzing Matched Datasets," Observational Studies 7(2); R package randChecks. He
    runs a design-as-null test (H0: treatment ~ the assignment mechanism implied by a
    chosen design) and places several designs on one univariate fixed-covariance
    Mahalanobis scale. He does NOT sweep resolution and his single-resolution fixed-metric
    framework does not exhibit the collapse. State clearly how your method differs from
    his design-as-null test and from his fixed-Sigma scale.

DELIVERABLE

Produce a methods proposal with:

  1. A precise definition of the new reference distribution and/or normalization: what is
     randomized or held out, what is held fixed, and the exact statistic and the exact
     calibration number a user would report.
  2. A PROOF or tight argument that it does NOT collapse: show formally what happens to
     your reported number as the match tightens (within-set var(X) -> 0) with residual
     imbalance held fixed. Demonstrate it does not drive the calibration to a degenerate
     limit.
  3. An explicit statement of WHICH trilemma corner(s) you sacrifice and why that
     sacrifice is acceptable here. Keeping corner (2) is mandatory; trading corner (3)
     (only toward analogous designs) and/or corner (1) (sample-splitting / held-out /
     design-sensitivity) is your choice. Name the trade.
  4. A demonstration that it is NOT l* and NOT a complete-randomization comparison:
     name the concrete difference.
  5. If overfitting is a threat to your reference class (e.g., the best feasible matching
     is selected to look good), show how your method neutralizes it.
  6. Computational feasibility: how to sample feasible matchings or compute the held-out
     quantity at realistic n; cost relative to one d^2 evaluation.
  7. A FALSIFICATION / VERIFICATION plan: a small simulation that would expose the method
     if it secretly collapses or secretly reduces to a forbidden comparison. Specify the
     data-generating setup, the sequence of progressively tighter matches, and the
     expected behavior of your reported number across that sequence (it should improve or
     stay stable as balance improves, never degrade).
  8. Honest limitations: assumptions it depends on, where it could mislead, and the one
     objection you find hardest to answer.

Lead with the single best proposal. If you have a strong runner-up that sacrifices a
different aspect, sketch it briefly at the end and say why you ranked it second. Show your
reasoning; do not just assert that the method works.

ORCHESTRATION (how to run this as a multi-agent search)

Do not answer with one pass of your own reasoning. Run a fan-out / judge / synthesize
search so that several independent inventions compete and the survivors are
adversarially stress-tested. Concretely:

  PHASE 1 -- DIVERGE (independent invention). Spawn k >= 5 inventor agents IN PARALLEL,
  each given this entire prompt and NOTHING about the others' work. To force diversity,
  assign each inventor a different starting lens so they do not converge on one idea:
    - inventor A: reference distribution over FEASIBLE MATCHINGS of the same pool
      (re-randomize the matching assignment, not the treatment label);
    - inventor B: ANTI-OVERFITTING / optimism correction (match on permuted or noise
      covariates; effective-degrees-of-freedom penalty for the matching procedure);
    - inventor C: SAMPLE-SPLITTING / CROSS-VALIDATED held-out-covariate balance;
    - inventor D: DESIGN-SENSITIVITY style (how large an unobserved within-set imbalance
      or how much selection optimism would overturn the verdict);
    - inventor E: EFFICIENCY-RATIO normalization against the achievable optimum (best
      feasible matching as lower bound, null/expected matching as upper bound);
    - inventors F+ (if k > 5): free lens -- invent something not on this list.
  Each inventor returns the full DELIVERABLE (items 1-8) for its proposal.

  PHASE 2 -- JUDGE (adversarial verification). For each surviving proposal, spawn 2-3
  independent SKEPTIC agents whose job is to BREAK it, not to praise it. Each skeptic
  must try to show one of: (a) it secretly collapses (re-derive item 2 and look for the
  degenerate limit the inventor missed); (b) it secretly reduces to a complete-
  randomization / matching-as-pruning comparison; (c) it is l* in disguise; (d) it can be
  gamed by cherry-picking the most-balanced feasible matching; (e) item 7's falsification
  simulation would actually fail. A proposal survives only if a majority of its skeptics
  cannot break it. Default the skeptic to "broken" when genuinely uncertain -- it is
  cheaper to discard a salvageable idea than to ship a collapsing one.

  PHASE 3 -- SYNTHESIZE. From the survivors, write one consolidated recommendation: the
  single strongest method, plus any genuinely distinct survivor as a runner-up. Where two
  survivors are complementary (e.g., a feasible-matching reference class PLUS an
  overfitting correction that closes skeptic objection (d)), graft them into one method
  and say so. The synthesis must carry forward, intact, the non-collapse argument (item 2)
  and the falsification plan (item 7) of whatever it recommends -- do not let the
  hardest-to-verify parts get lost in the merge.

  RUN THE FALSIFICATION SIM IF YOU CAN. If you have code execution, actually run item 7's
  simulation for the synthesized method (R preferred; the target package is RItools, but a
  standalone implementation is fine) and report whether the reported number stays stable
  or improves across progressively tighter matches. If you cannot execute code, write the
  simulation as runnable code and state exactly what output would falsify the method.

  REPORT FORMAT. Return: (1) the synthesized recommendation as a full deliverable; (2) a
  one-paragraph-each summary of every Phase-1 proposal and why it did or did not survive;
  (3) the surviving skeptic objections you could NOT fully answer, stated honestly;
  (4) WHICH LENS WON -- name the Phase-1 lens (A-F) the synthesized recommendation came
  from, or "graft" plus the contributing lenses if it merges several, in one line at the
  top so the result is comparable across separate runs.

  WRITE THE RESULT TO DISK. In addition to returning the report, write the full report to
  a Markdown file in the working directory named
  "collapse-method-invention-RESULT-<engine>.md" where <engine> is "claude" or "codex"
  depending on which system you are. Put the winning-lens line (item 4) as the first line
  of the file after the title, so two runs can be diffed at a glance. ASCII only, same
  style rules as this prompt.
