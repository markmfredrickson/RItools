##' Given covariates, a treatment variable and, optionally, a stratifying factor,
##' calculates standardized mean differences along each covariate,
##' with and without the stratification and tests for conditional
##' independence of the treatment variable and the covariates within
##' strata.
##'
##' In the unstratified case, the standardized difference of covariate
##' means is the mean in the treatment group minus the mean in the
##' control group, divided by the S.D. (standard deviation) in the
##' same variable estimated by pooling treatment and control group
##' S.D.s on the same variable.  In the stratified case, the
##' denominator of the standardized difference remains the same but
##' the numerator is a weighted average of within-stratum differences
##' in means on the covariate.  By default, each stratum is weighted
##' in proportion to the harmonic mean \eqn{1/[(1/a +
##' 1/b)/2]=2*a*b/(a+b)} of the number of treated units (a) and
##' control units (b) in the stratum; this weighting is optimal under
##' certain modeling assumptions (discussed in Kalton 1968, Hansen and
##' Bowers 2008).  This weighting can be modified using the
##' \code{stratum.weights} argument; see below. 
##'
##' When the treatment variable, the variable specified by the
##' left-hand side of \code{fmla}, is not binary, \code{xBalance}
##' calculates the covariates' regressions on the treatment variable,
##' in the stratified case pooling these regressions across strata
##' using weights that default to the stratum-wise sum of squared
##' deviations of the treatment variable from its stratum mean.
##' (Applied to binary treatment variables, this recipe gives the same
##' result as the one given above.)  In the numerator of the
##' standardized difference, we get a ``pooled S.D.'' from separating
##' units into two groups, one in which the treatment variable is 0 or
##' less and another in which it is positive.  If \code{report}
##' includes "adj.means", covariate means for the former of these
##' groups are reported, along with the sums of these means and the
##' covariates' regressions on either the treatment variable, in the
##' unstratified (``pre'') case, or the treatment variable and the
##' strata, in the stratified (``post'') case.
##'
##' \code{stratum.weights} can be either a function or a numeric
##' vector of weights.  If it is a numeric vector, it should be
##' non-negative and it should have stratum names as its names. (i.e.,
##' its names should be equal to the levels of the factor specified by
##' \code{strata}.) If it is a function, it should accept one
##' argument, a data frame containing the variables in \code{data} and
##' additionally \code{Tx.grp} and \code{stratum.code}, and return a
##' vector of non-negative weights with stratum codes as names; for an
##' example, do \code{getFromNamespace("harmonic", "RItools")}.
##'
##' If the stratifying factor has NAs, these cases are dropped.  On the other
##' hand, if NAs in a covariate are found then the default behavior is to
##' mean-impute while adding a dummy variable for whether NAs are found (and
##' checking balance for it as well).
##' 
##' If \code{covariate.scaling} is not \code{NULL}, no scaling is
##' applied. This behavior is likely to change in future versions.
##' (If you want no scaling, set \code{covariate.scaling=1}, as this
##' is likely to retain this meaning in the future.)
##'
##' @title Standardized Differences for Stratified Comparisons
##' @param fmla A formula containing an indicator of treatment
##'   assignment on the left hand side and covariates at right.
##' @param data A data frame in which \code{fmla} and \code{strata}
##'   are to be evaluated.
##' @param strata A list of right-hand-side-only formulas containing
##'   the factor(s) identifying the strata, with \code{NULL} entries
##'   interpreted as no stratification; or a factor with length equal
##'   to the number of rows in data; or a data frame of such
##'   factors. See below for examples.
##' @param report Character vector listing measures to report for each
##'   stratification; a subset of \code{c("adj.means",
##'   "adj.mean.diffs", "chisquare.test", "std.diffs", "z.scores",
##'   "p.values", "all")}. P-values reported are two-sided for the
##'   null-hypothesis of no effect. The option "all" requests all
##'   measures.
##' @param p.adjust.method Method of p-value adjustment.
##' @param element.weights Per-element weights. If there are clusters, the cluster weight is the sum of weights of elements within the cluster.  Within each stratum, cluster and element weights will be normalized to sum to 1.
##' @param stratum.weights Weights to be applied when aggregating
##'   across strata specified by \code{strata}, defaulting to weights
##'   proportional to the harmonic mean of treatment and control group
##'   sizes, in numbers of clusters if clusters are present, within strata.
##'   This can be either a function used to
##'   calculate the weights or the weights themselves; if
##'   \code{strata} is a data frame, then it can be such a function, a
##'   list of such functions, or a data frame of stratum weighting
##'   schemes corresponding to the different stratifying factors of
##'   \code{strata}.  See details.
##' @param subset Optional vector specifying a subset of observations to be used.
##' @param na.rm Whether to remove rows with NAs on any variables
##'   mentioned on the RHS of \code{fmla} (i.e. listwise deletion).
##'   Defaults to \code{FALSE}, wherein rows aren't deleted but for
##'   each variable with \code{NA}s a missing-data indicator variable
##'   is added to the variables on which balance is calculated and
##'   medians are imputed for the variable with missing data (in
##'   RItools versions 0.1-9 and before the default imputation was the
##'   mean, in RItools versions 0.1-11 and henceforth the default is
##'   the median). See the example below.
##' @param covariate.scaling A scale factor to apply to covariates in
##'   calculating \code{std.diffs}.  If \code{NULL}, \code{xBalance}
##'   pools standard deviations of each variable in the treatment and
##'   control group (defining these groups according to whether the
##'   LHS of \code{formula} is greater than or equal to 0).  Also, see
##'   details.
##' @param impfn A function to impute missing values when
##'   \code{na.rm=FALSE}. Currently \code{\link{median}}. To impute
##'   means use \code{\link{mean.default}}.
##' @param include.NA.flags Should NA flags be included?
##' @param post.alignment.transform Optional transformation applied to
##'   covariates just after their stratum means are subtracted off.
##' @return An object of class \code{c("xbal", "list")}.  There are
##'   \code{plot}, \code{print}, and \code{xtable} methods for class
##'   \code{"xbal"}; the \code{print} method is demonstrated in the
##'   examples.
##' @note Evidence pertaining to the hypothesis that a treatment
##'   variable is not associated with differences in covariate values
##'   is assessed by comparing the differences of means (or regression
##'   coefficients), without standardization, to their distributions
##'   under hypothetical shuffles of the treatment variable, a
##'   permutation or randomization distribution.  For the unstratified
##'   comparison, this reference distribution consists of differences
##'   (more generally, regression coefficients) when the treatment
##'   variable is permuted without regard to strata.  For the
##'   stratified comparison, the reference distribution is determined
##'   by randomly permuting the treatment variable within strata, then
##'   re-calculating the treatment-control differences (regressions of
##'   each covariate on the permuted treatment variable). Significance
##'   assessments are based on the large-sample Normal approximation
##'   to these reference distributions.
##' @export
##' @references Hansen, B.B. and Bowers, J. (2008), ``Covariate
##'   Balance in Simple, Stratified and Clustered Comparative
##'   Studies,'' \emph{Statistical Science} \bold{23}.
##'
##'   Kalton, G. (1968), ``Standardization: A technique to control for
##'   extraneous variables,'' \emph{Applied Statistics} \bold{17},
##'   118--136.
##' @author Ben Hansen and Jake Bowers and Mark Fredrickson
##' @keywords design nonparametric
##' @import SparseM svd
##' @examples
##' data(nuclearplants)
##' ##No strata, default output
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants)
##'
##' ##No strata, all output
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants,
##'          report=c("all"))
##'
##' ##Stratified, all output
##' xBalance(pr~.-cost-pt + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "adj.mean.diffs",
##'                   "chisquare.test", "std.diffs",
##'                   "z.scores", "p.values"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Missing data handling.
##' testdata<-nuclearplants
##' testdata$date[testdata$date<68]<-NA
##'
##' ##na.rm=FALSE by default
##' xBalance(pr ~ date, data = testdata, report="all")
##' xBalance(pr ~ date, data = testdata, na.rm = TRUE,report="all")
##'
##' ##To match versions of RItools 0.1-9 and older, impute means
##' #rather than medians.
##' ##Not run, impfn option is not implemented in the most recent version
##' \dontrun{xBalance(pr ~ date, data = testdata, na.rm = FALSE,
##'            report="all", impfn=mean.default)}
##'
##' ##Comparing unstratified to stratified, just one-by-one wilcoxon
##' #rank sum tests and omnibus test of multivariate differences on
##' #rank scale.
##' xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"),
##' 	 post.alignment.transform=rank)
xBalance <- function(fmla,
                     data,
                     strata = NULL,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs",
                         "chisquare.test","p.values", "all")[1:2],
                     #                     include.means=FALSE, chisquare.test=FALSE,
                     element.weights,
                     stratum.weights = harmonic,
                     subset,
                     na.rm = FALSE,
                     impfn = median,
                     include.NA.flags = TRUE,
                     covariate.scaling = NULL,
                     post.alignment.transform = NULL,
                     p.adjust.method = "holm") {

  if (!is.null(strata)) {
    stop("The strata argument has been deprecated. Use 'z ~ x1 + x2 + strata(s)' instead. See ?xBalance for details.")
  }

  stopifnot(is.null(post.alignment.transform) || is.function(post.alignment.transform))

  if (missing(data)) 
     data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "element.weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  if (cwpos <- match("element.weights", names(mf), nomatch=0))
      names(mf)[cwpos] <- "weights"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$na.action <- quote(stats::na.pass)
  data <- eval(mf, parent.frame())
  if (!cwpos)
      data$'(weights)' <- 1
  
  # Using charmatch instead of pmatch to distinguish between no match and ambiguous match. It reports
  # -1 for no match, and 0 for ambiguous (multiple) matches.
  valid.for.report <- c("adj.means","adj.mean.diffs","chisquare.test",
                                     "std.diffs","z.scores","p.values","all")
  report.good <- charmatch(report, valid.for.report, -1)

  if (any(report.good == -1)) {
    stop(paste("Invalid option(s) for report:", paste(report[report.good == -1], collapse=", ")))
  }
  if (any(report.good == 0)) {
    stop(paste("Option(s) for report match multiple possible values:", paste(report[report.good == 0], collapse=", ")))
  }

  # Now that we've found the partial matches, get their proper names
  report <- valid.for.report[report.good]

  if("all" %in% report)
    report <- c("adj.means","adj.mean.diffs","chisquare.test", "std.diffs","z.scores","p.values")

  design          <- makeDesign(fmla, data, imputefn = impfn, na.rm = na.rm, include.NA.flags = include.NA.flags)
  aggDesign       <- aggregateDesign(design)

  design <- as(design, "StratumWeightedDesign")
  design@Sweights <- DesignWeights(aggDesign, # Have to aggregate 1st to figure stratum weights properly
                                   effectOfTreatmentOnTreated) #For now we override any user-provided stratum.weights
  descriptives    <- designToDescriptives(design, covariate.scaling)

  # going forward, we use the user's weights, not ETT always
  aggDesign.weighted <- as(aggDesign, "StratumWeightedDesign")
  aggDesign.weighted@Sweights <-
      DesignWeights(aggDesign, stratum.weights)

  strataAligned <- alignDesignByStrata(aggDesign.weighted, post.alignment.transform)

  tmp <- lapply(strataAligned, function(i) { do.call(alignedToInferentials, i) })
  names(tmp) <- names(aggDesign@StrataMatrices)

  ans <- list()

  # append the z and p to the "descriptives" array (making it somewhat misnamed)
  tmp.z <- as.data.frame(lapply(tmp, function(tt) { tt$z }))
  tmp.p <- as.data.frame(lapply(tmp, function(tt) { tt$p }))
  nstats.previous <- dim(descriptives)[2]
  descriptives <- abind(descriptives, along = 2, tmp.z, tmp.p, use.first.dimnames = TRUE)
  names(dimnames(descriptives)) <- c("vars", "stat", "strata")

  dimnames(descriptives)[[2]][nstats.previous + 1:2] <- c("z", "p")

  inferentials <- do.call(rbind, lapply(tmp, function(s) {
    c(s$csq, s$DF, pchisq(s$csq, df = s$DF, lower.tail = FALSE))
  }))
  colnames(inferentials) <- c("chisquare", "df", "p.value")

  # the meat of our xbal object
  ans$overall <- inferentials
  ans$results <- descriptives

  # do p.value adjustment
  ans$results[, "p", ] <- p.adjust(ans$results[, "p", ], method = p.adjust.method)
  ans$overall[, "p.value"] <- p.adjust(ans$overall[, "p.value"], method = p.adjust.method)


  attr(ans$results, "originals") <- design@OriginalVariables
  attr(ans$overall, "tcov") <- lapply(tmp, function(r) {
    r$tcov
  })
  attr(ans, "fmla") <- formula(fmla)
  attr(ans, "report") <- report # hinting to our summary method later

  class(ans) <- c("xbal", "list")
  ans
}

xBalance.make.stratum.mean.matrix <- function(ss, mm) {
  stopifnot(inherits(ss, "factor")) # just in case a numeric variable is passed in.

  post.nobs <- dim(mm)[1]
  nlev <- nlevels(ss)

  # for this matrix, finding the indices of the rows is easy, as there is only one
  # item per row, and there post.nobs number of rows.
  tR <- new("matrix.csr",
            ja = as.integer(as.integer(ss)),
            ia = as.integer(1:(post.nobs+ 1)),
            ra = unsplit(1/tapply(ss,ss,length),ss),
            dimension = c(post.nobs,nlev))

  # With many items per row, we need to break the list of strata
  # down into indices of where each row starts and ends
  # e.g. Say ss = 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 the row indices would be
  # 1 5 12 16 (where 16 is the start of the non existant 4th row)

  L <- new("matrix.csr", #ifelse(oldver,"tripletMatrix","dgTMatrix"),
           ia = as.integer(1:(post.nobs + 1)),
           ja = as.integer(as.integer(ss)),
           ra = rep(1,length(ss)),
           dimension = c(post.nobs,nlev))


  msmn <- t(tR) %*% mm
  msmn <- L %*% msmn

  msmn <- as.matrix(msmn)

  return(msmn)
}
