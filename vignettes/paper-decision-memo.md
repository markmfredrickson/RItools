# Should the sigma_x / scale-invariance work become a short paper?

A discussion memo for the RItools team (Ben, Mark, Josh).
Jake Bowers, with Claude. 2026-06-04. Branch: `devel-sigma-x-omnibus`.

## The question

We have a memo (`vignettes/sigma-x-rosenbaum-memo.qmd`) and a package
update that grew out of Rosenbaum's criticism of the Hansen-Bowers (2008)
omnibus balance test. The memo's findings: HB08, and every balance test
built only from the matched sample, is scale-invariant. It cannot tell a
tight match from a loose one when only the within-stratum spread changes,
because studentization divides the absolute scale away. The fix is to
report standardized effect sizes against an external (pre-matching)
yardstick. The question on the table: is there a short paper here, or is
this documentation?

## How I got these answers

I asked two AI agents to role-play expert reviewers and judge whether the
work merits a paper. One simulated Paul Rosenbaum, reasoning from his
published positions on observational-study design and balance assessment.
One simulated Don Green, reasoning from his published positions on field
experiments and applied practice. **Both are simulations, not the real
people.** They are grounded in each scholar's books and papers, but
nothing here should be read as a statement from Rosenbaum or Green. I
treat their output the way I would treat a sharp graduate student's
steelman: useful for finding weaknesses while we can still fix them.

Each agent read the full memo and the project handoff before responding.

## Bottom line

Neither persona would write this as a stand-alone methods paper. Both
would put the material in the package documentation (a vignette) and, if
we want something citable, a short expository note or a software paper.
The two reviews were produced independently and converged on the same
verdict and most of the same reasons. That convergence is the strongest
signal in this memo.

## Where the two reviews agreed

Five points came up in both reviews without any coordination between
them.

1. **Drop "impossibility theorem."** The mathematics is one line: a
   statistic homogeneous in the within-stratum deviations, compared to a
   reference built from those same deviations, is scale-invariant,
   because multiplying every value in a list by a positive constant
   preserves ranks. Both reviewers said that calling this an
   impossibility theorem oversells a corollary of studentization and
   costs us credibility on the parts that are genuinely useful. State it
   plainly and claim it modestly.

2. **The core recommendation is received wisdom, and we do not cite it.**
   Standardized differences against a fixed pre-matching yardstick is the
   Cochran-Rubin-Stuart tradition and, centrally, the "balance test
   fallacy" of Imai, King, and Stuart (2008, JRSS-A). The memo cites none
   of them in its body. Both reviewers flagged this as a real problem:
   without those citations the memo reads as reinventing a settled idea
   and claiming credit for it. Any paper has to foreground that
   literature and contribute the one new thing on top of it.

3. **Cut the chi^2_p combined p-value.** This is the sharpest agreement,
   and it answers one of our open team questions. A single pool-calibrated
   p-value reintroduces exactly the misleading scalar the whole document
   exists to retire. The Rosenbaum persona: it re-imports the confusion we
   just spent ten pages dispelling. The Green persona: a practitioner will
   fixate on the new number for the same reason she fixated on the old
   one, and declining to pick a side reads as unresolved. If we want one
   summary number, report the Mahalanobis distance d' Sigma_pool^{-1} d in
   pool-SD units as a distance, with no p-value attached.

4. **The eight-unit toy cannot be the only evidence.** Three coincident
   controls plus one offset treated unit per stratum is engineered so that
   the only thing varying with delta is a pure rescaling. A scale-
   invariant test is then trivially flat on it, which a referee will read
   as assuming the conclusion. Both reviewers said: keep the toy as
   pedagogy, but the paper needs a real matched study where we tighten the
   match and show the omnibus p-value sitting still while the standardized
   differences shrink.

5. **The genuinely new contribution is narrow.** It is not that HB08
   fails, since HB08 was never offered as a measure of design quality. It
   is that *no* choice of balance statistic computed on the matched sample
   escapes the invariance: rank-based, distance-based, energy, all of it.
   That cataloguing, together with the point that the invariance is
   *exact* rather than the asymptotic underpowering that Imai-King-Stuart
   describe, is the only thing worth building a note around.

## What each review added on its own

**The Rosenbaum persona: the Chapter 6 critique is a strawman.** The memo
attacks Rosenbaum's complete-randomization benchmark for a workflow he
never proposed, namely recomputing it as the caliper varies and reading
the movement as a quality score. The persona's objection was blunt: the
memo commits, against his chapter, the same error it indicts in others
(using a tool for a job it was not built for). The demonstration is also
weak on its own terms, moving a p-value from 0.90 to 0.89 in an example
the memo admits was not tuned for drama. Either back the caliper claim
with a consequential example or cut it. The persona also noted that he
and the memo agree on the destination, so framing the memo as a correction
of Rosenbaum misreads the relationship.

**The Green persona: separate the two audiences.** `balanceTest()` serves
both experimenters and observational researchers, and the result lands
differently on each. For an experimenter, a balance test is a check that
randomization was implemented correctly. Scale-invariance does not impair
that use and is arguably a feature, since she does not want her integrity
check to depend on the units in which age happens to be recorded. The
result bites only the observational researcher who reads the balance
p-value as a quality grade, and for her it reinforces advice the matching
tutorials already give. The persona's distinct warning: the "reward small
absolute imbalance" framing, if applied to an experiment, could nudge
someone toward conditioning on realized balance, which is the error that
corrupts experimental inference. A paper needs an explicit paragraph
saying this critique is about matching, not a license for experimenters to
grade or select among realized randomizations.

## My read, including the counterweight

I think the reviews are right but slightly too quick to dismiss one thing.
The exact-versus-asymptotic distinction is sharper than Imai-King-Stuart
and is the quotable core: the feature that makes HB08 a valid
randomization test, studentization by the design's own variance, is the
same feature that disqualifies it as a measure of match quality.
Imai-King-Stuart say that significance conflates imbalance with sample
size. Our point is different and cleaner: the randomization-calibrated
balance test is blind to absolute imbalance by construction, to the last
decimal. That, plus "no statistic escapes," is a real contribution for the
randomization-inference balance-testing lineage that HB08 and RItools
represent, which the regression-based framing of Imai-King-Stuart does not
cover.

So the honest contribution survives. It is modest, and it survives only if
we make all of the changes above: drop the impossibility language, credit
the prior literature up front, fix or cut the Chapter 6 caliper claim, cut
the chi^2_p p-value, lead with real data, and separate the two audiences.
That is a long list of conditions on a thin result, which is why both
reviews land on "vignette plus short note," not "paper."

## Decisions this memo asks the team to make

1. **Paper, note, or vignette?** My recommendation, following both
   reviews: build the vignette regardless, since that is the strongest
   home for the material, and decide separately whether to also write a
   short note for a citable artifact.

2. **Kill the chi^2_p combined p-value?** Both reviews say yes. This was
   open question 2 from the handoff. I lean toward cutting it from the
   user-facing output and reporting the pool-standardized distance
   instead. This is the most consequential code decision in front of us.

3. **Venue, if we write the note.** Observational Studies (short note,
   Rosenbaum's own journal, engages Chapter 6 directly) or a software
   paper in the R Journal or JSS whose contribution is that RItools now
   reports pool-calibrated effect sizes and documents why the omnibus
   p-value should not be read as match quality. The Green persona ranked
   the software paper highest on the grounds that it is what gets cited.
   Neither persona would send it to a theory venue, and the Green persona
   would not send it to Political Analysis.

4. **The Chapter 6 caliper claim.** Back it with a consequential example
   or demote it to a footnote. As written it is both a strawman against
   Rosenbaum and numerically negligible.

## Suggested next step

If we want the citable artifact, the natural move is a reframed outline
for a roughly five-page Observational Studies note built on this advice:
lead with the exact-invariance and no-statistic-escapes contribution,
credit Imai-King-Stuart and the standardized-difference tradition up
front, use a real matched study, separate the two audiences, and drop the
chi^2_p p-value. I can draft that outline, or draft the vignette version,
or push either persona further on a specific point (for example, what a
consequential caliper example would have to look like).
