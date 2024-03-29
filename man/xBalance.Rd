% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xBalance.R
\name{xBalance}
\alias{xBalance}
\title{Standardized Differences for Stratified Comparisons}
\usage{
xBalance(
  fmla,
  strata = list(unstrat = NULL),
  data,
  report = c("std.diffs", "z.scores", "adj.means", "adj.mean.diffs",
    "adj.mean.diffs.null.sd", "chisquare.test", "p.values", "all")[1:2],
  stratum.weights = harmonic,
  na.rm = FALSE,
  covariate.scaling = NULL,
  normalize.weights = TRUE,
  impfn = median,
  post.alignment.transform = NULL,
  pseudoinversion_tol = .Machine$double.eps
)
}
\arguments{
\item{fmla}{A formula containing an indicator of treatment
assignment on the left hand side and covariates at right.}

\item{strata}{A list of right-hand-side-only formulas containing
the factor(s) identifying the strata, with \code{NULL} entries
interpreted as no stratification; or a factor with length equal
to the number of rows in data; or a data frame of such
factors. See below for examples.}

\item{data}{A data frame in which \code{fmla} and \code{strata}
are to be evaluated.}

\item{report}{Character vector listing measures to report for each
stratification; a subset of \code{c("adj.means",
"adj.mean.diffs", "adj.mean.diffs.null.sd", "chisquare.test",
"std.diffs", "z.scores", "p.values", "all")}. P-values reported
are two-sided for the null-hypothesis of no effect. The option
"all" requests all measures.}

\item{stratum.weights}{Weights to be applied when aggregating
across strata specified by \code{strata}, defaulting to weights
proportional to the harmonic mean of treatment and control group
sizes within strata.  This can be either a function used to
calculate the weights or the weights themselves; if
\code{strata} is a data frame, then it can be such a function, a
list of such functions, or a data frame of stratum weighting
schemes corresponding to the different stratifying factors of
\code{strata}.  See details.}

\item{na.rm}{Whether to remove rows with NAs on any variables
mentioned on the RHS of \code{fmla} (i.e. listwise deletion).
Defaults to \code{FALSE}, wherein rows aren't deleted but for
each variable with \code{NA}s a missing-data indicator variable
is added to the variables on which balance is calculated and
medians are imputed for the variable with missing data (in
RItools versions 0.1-9 and before the default imputation was the
mean, in RItools versions 0.1-11 and henceforth the default is
the median). See the example below.}

\item{covariate.scaling}{A scale factor to apply to covariates in
calculating \code{std.diffs}.  If \code{NULL}, \code{xBalance}
pools standard deviations of each variable in the treatment and
control group (defining these groups according to whether the
LHS of \code{formula} is greater than or equal to 0).  Also, see
details.}

\item{normalize.weights}{If \code{TRUE}, then stratum weights are
normalized so as to sum to 1.  Defaults to \code{TRUE}.}

\item{impfn}{A function to impute missing values when
\code{na.rm=FALSE}. Currently \code{\link{median}}. To impute
means use \code{\link{mean.default}}.}

\item{post.alignment.transform}{Optional transformation applied to
covariates just after their stratum means are subtracted off.}

\item{pseudoinversion_tol}{The function uses a singular value decomposition
to invert a covariance matrix. Singular values less than this tolerance
will be treated as zero.}
}
\value{
An object of class \code{c("xbal", "list")}.  There are
  \code{plot}, \code{print}, and \code{xtable} methods for class
  \code{"xbal"}; the \code{print} method is demonstrated in the
  examples.
}
\description{
Given covariates, a treatment variable, and a stratifying factor,
calculates standardized mean differences along each covariate,
with and without the stratification and tests for conditional
independence of the treatment variable and the covariates within
strata.
}
\details{
Note: the newer \code{\link{balanceTest}} function provides the same
functionality as \code{xBalance} with additional support for clustered
designs. While there are no plans to deprecate \code{xBalance}, users are
encouraged to use \code{balanceTest} going forward.

In the unstratified case, the standardized difference of covariate
means is the mean in the treatment group minus the mean in the
control group, divided by the S.D. (standard deviation) in the
same variable estimated by pooling treatment and control group
S.D.s on the same variable.  In the stratified case, the
denominator of the standardized difference remains the same but
the numerator is a weighted average of within-stratum differences
in means on the covariate.  By default, each stratum is weighted
in proportion to the harmonic mean \eqn{1/[(1/a +
1/b)/2]=2*a*b/(a+b)} of the number of treated units (a) and
control units (b) in the stratum; this weighting is optimal under
certain modeling assumptions (discussed in Kalton 1968, Hansen and
Bowers 2008).  This weighting can be modified using the
\code{stratum.weights} argument; see below.

When the treatment variable, the variable specified by the
left-hand side of \code{fmla}, is not binary, \code{xBalance}
calculates the covariates' regressions on the treatment variable,
in the stratified case pooling these regressions across strata
using weights that default to the stratum-wise sum of squared
deviations of the treatment variable from its stratum mean.
(Applied to binary treatment variables, this recipe gives the same
result as the one given above.)  In the numerator of the
standardized difference, we get a ``pooled S.D.'' from separating
units into two groups, one in which the treatment variable is 0 or
less and another in which it is positive.  If \code{report}
includes "adj.means", covariate means for the former of these
groups are reported, along with the sums of these means and the
covariates' regressions on either the treatment variable, in the
unstratified (``pre'') case, or the treatment variable and the
strata, in the stratified (``post'') case.

\code{stratum.weights} can be either a function or a numeric
vector of weights.  If it is a numeric vector, it should be
non-negative and it should have stratum names as its names. (i.e.,
its names should be equal to the levels of the factor specified by
\code{strata}.) If it is a function, it should accept one
argument, a data frame containing the variables in \code{data} and
additionally \code{Tx.grp} and \code{stratum.code}, and return a
vector of non-negative weights with stratum codes as names; for an
example, do \code{getFromNamespace("harmonic", "RItools")}.

If \code{covariate.scaling} is not \code{NULL}, no scaling is
applied. This behavior is likely to change in future versions.
(If you want no scaling, set \code{covariate.scaling=1}, as this
is likely to retain this meaning in the future.)

\code{adj.mean.diffs.null.sd} returns the standard deviation of
the Normal approximated randomization distribution of the
strata-adjusted difference of means under the strict null of no
effect.
}
\note{
Evidence pertaining to the hypothesis that a treatment
  variable is not associated with differences in covariate values
  is assessed by comparing the differences of means (or regression
  coefficients), without standardization, to their distributions
  under hypothetical shuffles of the treatment variable, a
  permutation or randomization distribution.  For the unstratified
  comparison, this reference distribution consists of differences
  (more generally, regression coefficients) when the treatment
  variable is permuted without regard to strata.  For the
  stratified comparison, the reference distribution is determined
  by randomly permuting the treatment variable within strata, then
  re-calculating the treatment-control differences (regressions of
  each covariate on the permuted treatment variable). Significance
  assessments are based on the large-sample Normal approximation
  to these reference distributions.
}
\examples{
data(nuclearplants)
##No strata, default output
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         data=nuclearplants)

##No strata, all output
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         data=nuclearplants,
         report=c("all"))

##Stratified, all output
xBalance(pr~.-cost-pt, strata=factor(nuclearplants$pt),
         data=nuclearplants,
         report=c("adj.means", "adj.mean.diffs",
                  "adj.mean.diffs.null.sd",
                  "chisquare.test", "std.diffs",
                  "z.scores", "p.values"))

##Comparing unstratified to stratified, just adjusted means and
#omnibus test
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=list(unstrat=NULL, pt=~pt),
         data=nuclearplants,
         report=c("adj.means", "chisquare.test"))

##Comparing unstratified to stratified, just adjusted means and
#omnibus test
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=data.frame(unstrat=factor('none'),
           pt=factor(nuclearplants$pt)),
         data=nuclearplants,
         report=c("adj.means", "chisquare.test"))

##Missing data handling.
testdata<-nuclearplants
testdata$date[testdata$date<68]<-NA

##na.rm=FALSE by default
xBalance(pr ~ date, data = testdata, report="all")
xBalance(pr ~ date, data = testdata, na.rm = TRUE,report="all")

##To match versions of RItools 0.1-9 and older, impute means
#rather than medians.
##Not run, impfn option is not implemented in the most recent version
\dontrun{xBalance(pr ~ date, data = testdata, na.rm = FALSE,
           report="all", impfn=mean.default)}

##Comparing unstratified to stratified, just one-by-one wilcoxon
#rank sum tests and omnibus test of multivariate differences on
#rank scale.
xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=data.frame(unstrat=factor('none'),
           pt=factor(nuclearplants$pt)),
         data=nuclearplants,
         report=c("adj.means", "chisquare.test"),
	 post.alignment.transform=rank)
}
\references{
Hansen, B.B. and Bowers, J. (2008), ``Covariate
  Balance in Simple, Stratified and Clustered Comparative
  Studies,'' \emph{Statistical Science} \bold{23}.

  Kalton, G. (1968), ``Standardization: A technique to control for
  extraneous variables,'' \emph{Applied Statistics} \bold{17},
  118--136.
}
\seealso{
\code{\link{balanceTest}}
}
\author{
Ben Hansen and Jake Bowers and Mark Fredrickson
}
\keyword{design}
\keyword{nonparametric}
