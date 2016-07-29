
##' Given an assignment variable, variables on which to compare
##' groups assigned to treatment and control conditions and,
##' optionally, a clustering variable and/or
##' one or more stratifying factors,
##' compares univariate and multivariate measures suitable for comparing
##' the two groups and for testing the proposition that assignment
##' was random or effectively random within levels of a stratifying factor. 
##'
##' The function assembles various univariate descriptive statistics
##' for the groups to be compared: (weighted) means of treatment and
##' control groups; differences of these (adjusted differences); and
##' adjusted differences as multiples of a pooled S.D. of the variable
##' in the treatment and control groups (standard differences). This
##' is done separately for each provided stratifying factor and, by
##' default, for the unstratified comparison, in each case reflecting
##' a standardization appropriate to the designated (post-)
##' stratification of the sample.  In the case without stratification
##' or clustering, the only weighting used to calculate treatment and
##' control group means is that provided by the user as an
##' \code{element.weight}; in the absence of such an argument, these
##' means are unweighted.  When there are strata, within-stratum means
##' of treatment or of control observations are calculated using
##' \code{element.weights}, if provided, and then these are combined
##' across strata according to a \sQuote{effect of treatment on
##' treated}-type weighting scheme. (The function's
##' \code{stratum.weights} argument figures in the function's
##' inferential calculations but not these descriptive calculations.)
##' To figure a stratum's effect of treatment on treated weight, the
##' sum of all \code{element.weights} associated with treatment or
##' control group observations within the stratum is multiplied by the
##' fraction of clusters in that stratum that are associated with the
##' treatment rather than the control condition.  (Unless this
##' fraction is 0 or 1, in which case the stratum is downweighted to
##' 0.)
##'
##' The function also calculates univariate and multivariate inferential
##' statistics, targeting the hypothesis that assignment was random within strata. These
##' calculations also pool \code{element.weight}-ed, within-stratum group means across strata,
##' but the default weighting of strata differs from that of the descriptive calculations, and is
##' determined by the \code{stratum.weights} argument. By default, each stratum is weighted
##' in proportion to the product of the stratum mean of \code{element.weight}s and
##' the harmonic mean \eqn{1/[(1/a +
##' 1/b)/2]=2*a*b/(a+b)} of the number of treated units (a) and
##' control units (b) in the stratum; this weighting is optimal under
##' certain modeling assumptions (discussed in Kalton 1968, Hansen and
##' Bowers 2008).  The multivariate assessment is based on a Mahalanobis-type
##' distance that combines each of the univariate mean differences while accounting
##' for correlations among them. It's similar to the Hotelling's T-squared statistic,
##' except standarized using a permutation covariance.  See Hansen and Bowers (2008).  
##'
##' In contrast to the earlier function \code{xBalance} that it is intended to replace,
##' \code{balanceTest} accepts only binary assignment variables (for now). 
##'
##' \code{stratum.weights} can be either a function or a numeric
##' vector of weights.  If it is a numeric vector, it should be
##' non-negative and it should have stratum names as its names. (i.e.,
##' its names should be equal to the levels of the factor specified by
##' \code{strata}.) If it is a function, it should accept one
##' argument, a data frame containing the variables in \code{data} and
##' additionally \code{Tx.grp} and \code{stratum.code}, and return a
##' vector of non-negative weighting factors with stratum codes as names.
##' (To see the function that's applied by default, 
##' do \code{getFromNamespace("harmonic", "RItools")}.)  These weighting factors
##' will be multipled by the stratum mean of \code{element.weights} to determine
##' the stratum weights used for inferential calculations.
##'
##' If the stratifying factor has NAs, these cases are dropped.  On the other
##' hand, if NAs in a covariate are found then those observations are dropped for descriptive
##' calculations and "imputed" to the stratum mean of the variable for inferential calculations. 
##' When covariate values are dropped due to missingness, proportions of observations not missing on
##' that variable are recorded and returned.  The printed output presents non-missing proportions alongside of
##' the variables themselves, distinguishing the former by placing them at the bottom of the list and enclosing the
##' variable's name in parentheses.  If a variable shares a missingness pattern with other another variable,
##' its missingness information may be labeled with the name of the other variable in the output. 
##' 
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
##' @param element.weights Per-element weight, or 0 if element does not meet condition specified by subset argument. If there are clusters, the cluster weight is the sum of weights of elements within the cluster.  Within each stratum, cluster and element weights will be normalized to sum to 1.
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
##' @param subset Optional; condition or vector specifying a subset of observations to be given positive element weights.
##' @param covariate.scaling A scale factor to apply to covariates in
##'   calculating \code{std.diffs} (currently ignored). 
##' @param include.NA.flags Present item missingness comparisons as well as covariates themselves?
##' @param post.alignment.transform Optional transformation applied to
##'   covariates just after their stratum means are subtracted off.
##' @return An object of class \code{c("xbal", "list")}.  There are
##'   \code{plot}, \code{print}, and \code{xtable} methods for class
##'   \code{"xbal"}; the \code{print} method is demonstrated in the
##'   examples.
##' @note Evidence pertaining to the hypothesis that a treatment
##'   variable is not associated with differences in covariate values
##'   is assessed by comparing the differences of means, without standardization, to their distributions
##'   under hypothetical shuffles of the treatment variable, a
##'   permutation or randomization distribution.  For the unstratified
##'   comparison, this reference distribution consists of differences
##'   as the treatment
##'   assignments of clusters are freely permuted.  For
##'   stratified comparisons, the reference distributions describes re-randomizations of
##'   this type performed separately in each stratum. Significance
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
##' @import svd stats
##' @examples
##' data(nuclearplants)
##' ##No strata, default output
##' balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants)
##'
##' ##No strata, all output
##' balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
##'          data=nuclearplants,
##'          report=c("all"))
##'
##' ##Stratified, all output
##' balanceTest(pr~.-cost-pt + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "adj.mean.diffs",
##'                   "chisquare.test", "std.diffs",
##'                   "z.scores", "p.values"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Comparing unstratified to stratified, just adjusted means and
##' #omnibus test
##' balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"))
##'
##' ##Missing data handling.
##' testdata<-nuclearplants
##' testdata$date[testdata$date<68]<-NA
##'
##'
##'
##' ##Comparing unstratified to stratified, just one-by-one wilcoxon
##' #rank sum tests and omnibus test of multivariate differences on
##' #rank scale.
##' balanceTest(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
##'          data=nuclearplants,
##'          report=c("adj.means", "chisquare.test"),
##' 	 post.alignment.transform=rank)
balanceTest <- function(fmla,
                     data,
                     strata = NULL,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs",
                         "chisquare.test","p.values", "all")[1:2],
                     #                     include.means=FALSE, chisquare.test=FALSE,
                     element.weights,
                     stratum.weights = harmonic,
                     subset,
                     include.NA.flags = TRUE,
                     covariate.scaling = NULL,
                     post.alignment.transform = NULL,
                     p.adjust.method = "holm") {
### API Assumptions:
### - no ... in the xBal formula
### (if this assumption ceases to be met then we have to add an explicit check that
### the user hasn't tried to specify an offset, given that we're repurposing that
### model.frame option)

  if (!is.null(strata)) {
    stop("The strata argument has been deprecated. Use 'z ~ x1 + x2 + strata(s)' instead. See ?balanceTest, examples.")
  }

  stopifnot(is.null(post.alignment.transform) || is.function(post.alignment.transform))

  if (missing(data)) 
     data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "element.weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  if (cwpos <- match("element.weights", names(mf), nomatch=0L))
      names(mf)[cwpos] <- "weights"
  ## Here's where we rely on assumption of no ... in the xBal formula
  ## it helps us avoid adding a second offset argument to the model frame call
  if (sspos <- match("subset", names(mf), nomatch=0L))
      names(mf)[sspos] <- "offset"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$na.action <- quote(stats::na.pass)
  data <- eval(mf, parent.frame())
  if (!cwpos)
      data$'(weights)' <- 1
  if (sspos)
      {
          ss <- data$'(offset)'
          ss <- as.logical(ss)
          if (any(is.na(ss)))
              {
                  ss[is.na(ss)] <- FALSE
                  warning("subset specification gave NAs; interpreting these as FALSE")
              }
          data$'(weights)' <- ifelse(ss, data$'(weights)', 0)
          data$'(offset)' <- NULL
      }

  if (any(is.na(data$'(weights)'))) stop("NAs detected in weights")
  
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

  design          <- makeDesigns(fmla, data)

  ## We need aggDesign to make descriptives because the stratum weights need to be calculated from the aggregated setup
  ## however we're not going to feed aggDesign itself to `designToDescriptives()`, so we don't risk "polluting"
  ## the descriptives calcs w/imputation that's being done along with the aggregation procedure
  aggDesign       <- aggregateDesigns(design)

  design <- as(design, "StratumWeightedDesignOptions")
  design@Sweights <- DesignWeights(aggDesign, # Have to aggregate 1st to figure stratum weights properly
                                   effectOfTreatmentOnTreated) #For now we override any user-provided stratum.weights
  descriptives    <- designToDescriptives(design, covariate.scaling)

  # going forward, we use the user's weights, not ETT always
  aggDesign.weighted <- as(aggDesign, "StratumWeightedDesignOptions")
  aggDesign.weighted@Sweights <-
      DesignWeights(aggDesign, stratum.weights)

  strataAligned <- alignDesignsByStrata(aggDesign.weighted, post.alignment.transform)
  origvars <- strataAligned[[1]]@OriginalVariables #to include NotMissing columns
  
  tmp <- lapply(strataAligned, alignedToInferentials)
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

  ## do p.value adjustment
  for (s in 1L:dim(descriptives)[3])
  ans$results[, "p", s] <- p.adjust(ans$results[, "p", s], method = p.adjust.method)
##  ans$overall[, "p.value"] <- p.adjust(ans$overall[, "p.value"], method = p.adjust.method)

  attr(ans$results, "originals") <- origvars
  attr(ans$results, "term.labels") <- design@TermLabels
  attr(ans$results, "include.NA.flags") <- include.NA.flags # hinting for print and plot methods

  attr(ans$overall, "tcov") <- lapply(tmp, function(r) {
    r$tcov
  })
  attr(ans, "fmla") <- formula(fmla)
  attr(ans, "report") <- report # hinting to our summary method later
  class(ans) <- c("xbal", "list")
  ans
}
