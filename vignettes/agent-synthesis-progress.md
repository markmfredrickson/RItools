# Agent synthesis progress (started 2026-04-09 11pm CT)

This file tracks the progress of the agent-based mathematical-statistics
analysis dispatched at 11:04pm CT on April 9, 2026. Claude will update it
as agents return. Check back in the morning.

## Agent status

- **Agent 2 (ac4d29972fc0d3586)**: COMPLETED. Full report received.
- **Agent 1 (aea526da0249fab48)**: Still running as of 11:55pm CT.

## Key findings so far (from Agent 2's report)

### Task 5: Impossibility proof

Agent 2 wrote a clean theorem-statement version of the scale-invariance
result:

> **Theorem.** Let T(X, z, s) be any test statistic satisfying
> T(X_delta, z, s) = g(delta) * T(X_1, z, s) for all z in Z, where
> g is strictly positive and does not depend on z. Then the rank-based
> p-value is invariant to delta.

The key insight: what matters is not that T is homogeneous in X, but that
the SAME g(delta) applies to every z. This makes every element of the
permutation distribution rescale by the same factor, preserving ranks.

### Task 1: Hunting for escaping statistics

Agent 2 checked KS, Anderson-Darling, Wasserstein, energy, entropy, HSIC.
Findings:

- **KS and Anderson-Darling**: MORE degenerate than HB08. They measure
  vertical CDF gaps, which are exactly 1 on this toy regardless of delta
  (until delta = 0 exactly). They don't see the perturbation at all.

- **Wasserstein**: 1-homogeneous. Falls to the theorem. Confirmed.

- **Energy statistic (Szekely-Rizzo)**: The theorem does NOT formally
  cover this because the energy stat mixes within-stratum (delta-sensitive)
  and between-stratum (constant) pairwise distances, and different z's
  produce different functional dependencies on delta. But Agent 2
  conjectured the rank probably doesn't move because between-stratum
  dominates. **I verified this numerically**: rank is constant at 13/28
  under complete randomization and 1/16 under within-stratum, across
  every delta from 1 to 0.001. The energy statistic does NOT escape.

- **HSIC with median-heuristic bandwidth**: Falls to the theorem because
  the bandwidth scales with the data.

- **HSIC with FIXED (external) bandwidth**: ESCAPES. But only because the
  bandwidth is external information --- confirming the memo's taxonomy
  that external information is needed.

- **Multinomial entropy with fixed bins**: Would escape but only because
  the bin edges are external. Trades one pathology for another (arbitrary
  bin choice).

**Bottom line for Task 1**: No test statistic built entirely from the
matched sample escapes the scale invariance on this toy. The only
candidates that escape do so by importing external scale information
(fixed bandwidth, fixed bin edges), which confirms rather than refutes
the memo's argument.

### Task 2: Pre-matching pool as reference

Agent 2 argues that caliper-invariance requires BOTH the standardization
AND the reference distribution to come from the pool. Supplying only
Sigma_x from the pool while keeping the within-stratum randomization null
leaves the scale-invariance. Agent 2 also argues the user's fullmatch-pool
idea, correctly implemented, amounts to the same maneuver as supplying an
external Sigma_x from the pool and using a fixed reference distribution.

### Disagreement with the memo

Agent 2 disagrees that option 3 (chi^2_p with external V*) is the
"strongest assumption" among the three external-reference options. Agent 2
argues that chi^2_p with pool-derived V* is really just option (a) (user-
supplied pool-based Sigma_x) dressed up differently, and the memo should
embrace it as the natural p-value companion to effect-size reporting.

## Energy statistic verification (done by me, not the agents)

```
--- Complete randomization (28 assignments, 2 of 8 treated) ---
delta      T_obs_energy   rank           p
1          -4.32426       13/28          0.4643
0.5        -4.55689       13/28          0.4643
0.1        -4.7503        13/28          0.4643
0.01       -4.795         13/28          0.4643
0.001      -4.7995        13/28          0.4643

--- Within-stratum randomization (16 assignments, 1 per stratum) ---
delta      T_obs_energy   rank           p
1          -4.32426       1/16           0.0625
0.5        -4.55689       1/16           0.0625
0.1        -4.7503        1/16           0.0625
0.01       -4.795         1/16           0.0625
0.001      -4.7995        1/16           0.0625
```

Constant ranks across all deltas under both reference distributions.
Energy statistic does NOT escape the scale invariance on this toy, despite
the formal theorem not covering it.

## Still waiting for

Agent 1's report. Will update this file when it arrives and write the
full synthesis below.
