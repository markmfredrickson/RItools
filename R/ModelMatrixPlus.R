###############################################################################
## ModelMatrixPlus Objects: covariates with term-specific missingness info & indexing
################################################################################

setClassUnion("Contrasts", c("list", "NULL"))
##' ModelMatrixPlus S4 class
##'
##' If the Covariates matrix has an intercept, it will only be in the first column.
##'
##' More on NotMissing slot: It's matrix of numbers in [0,1].
##' First col is entirely \code{TRUE} or 1, like an intercept, unless corresponding
##' UnitWeight is 0, in which case it may also be 0 (see below). Subsequent cols
##' present only if there are missing covariate values, in which case these cols are
##' named for terms (of the original calling formula or data frame) that possess
##' missing values.  Terms with the same missing data pattern are mapped to a single
##' column of this matrix.  If the ModelMatrixPlus is representing elements, each column should
##' be all 1s and 0s, indicating which elements have non-missing values for the term
##' represented by that column.  If the ModelMatrixPlus as a whole represents clusters,
##' then there can be fractional values, but that situation should only arise in the
##' DesignOptions class extension of this class, so it's documented there.
##'
##' @slot OriginalVariables look-up table associating Covariates columns with terms of the originating model formula
##' @slot TermLabels labels of terms of the originating model formula
##' @slot contrasts Contrasts, a list of contrasts or NULL, as returned by `model.matrix.default`
##' @slot NotMissing Matrix of numbers in [0,1] with as many rows as the Covariates table but only one more col than there are distinct covariate missingness patterns (at least 1, nothing missing). First col is entirely T or 1, like an intercept.
##' @slot NM.Covariates integer look-up table mapping Covariates columns to columns of NotMissing.  (If nothing missing for that column, this is 0.)
##' @slot NM.terms integer look-up table mapping term labels to columns of NotMissing (0 means nothing missing in that column)
##' @slot UnitWeights vector of weights associated w/ rows of the ModelMatrixPlus
##' @keywords internal
setClass("ModelMatrixPlus",
         contains = "matrix",
         slots    = c(
             OriginalVariables = "integer",
             TermLabels        = "character",
             Contrasts         = "list",
             NotMissing        = "matrix",
             NM.Covariates     = "integer",
             NM.terms          = "integer",
             UnitWeights       = "numeric" )
         )

#' @method as.matrix ModelMatrixPlus
#' @export
as.matrix.ModelMatrixPlus <- function(x, ...)
    {
        ans <- x@.Data
        attr(ans, "assign") <- x@OriginalVariables
###        attr(ans, "term.labels") <- # a model.matrix wouldn't really
###            x@TermLabels # have this, although maybe it should
        attr(ans, "contrasts") <- x@Contrasts
        ans
    }

##' Grow a model matrix while at the same time compactly
##' encoding missingness patterns in RHS variables of a model frame.
##'
##'
##' @title Model matrices along with compact encodings of data availability/missingness
##' @param object Model formula or terms object (as in `model.matrix`)
##' @param data data.frame, as in `model.matrix()` but has to have \sQuote{\code{(weights)}} column
##' @param remove.intercept logical
##' @param ... passed to `model.matrix.default` (and further)
##' @return ModelMatrixPlus, i.e. model matrix enriched with missing data info
##' @author Ben B Hansen
##' @keywords internal
##' 
model_matrix <- function(object, data = environment(object), remove.intercept=TRUE, ...) {
  # mf <- model.frame(object, data, na.action = na.pass)
  tms <- terms(object)
  term.labels <- attr(tms, "term.labels")

  uweights <- as.vector(model.weights(data))
  if (is.null(uweights))
    stop("model_matrix() expects its data arg to be a model frame containing weights")
  stopifnot(is.numeric(uweights), all(!is.na(uweights)), all(uweights>=0))
  
  covariates <- model.matrix(object = object, data = data, ...) 

  assign <- attr(covariates, "assign")
  attr(covariates, "assign") <- NULL
  stopifnot(all(assign %in% 0:length(term.labels)))

  contrasts <- attr(covariates, "contrasts")
  attr(covariates, "contrasts") <- NULL

  if (remove.intercept &&
      (iloc <- match("(Intercept)", colnames(covariates), nomatch=0))
      ) {
    assign <- assign[-iloc]
    covariates <- covariates[,-iloc, drop=FALSE]
  }

    logicalsT <- (paste0(term.labels, "TRUE") %in% colnames(covariates)) &
        !(paste0(term.labels, "TRUE") %in% term.labels)
    logicalsF <- (paste0(term.labels, "FALSE") %in% colnames(covariates)) &
        !(paste0(term.labels, "FALSE") %in% term.labels)
    if (any(logicalsF & logicalsT))
    {
        rlocs <- colnames(covariates) %in% paste0(term.labels[logicalsF & logicalsT], "FALSE")
        assign <- assign[!rlocs]
        covariates <- covariates[,!rlocs, drop=FALSE]
    }
    if (any(logicalsT))
    {
        slocs <- colnames(covariates) %in% paste0(term.labels[logicalsT], "TRUE")
        colnames(covariates)[slocs] <-
            substr(colnames(covariates)[slocs], 1L,
                   nchar(colnames(covariates)[slocs]) - 4L)
        }

  cols.by.term <- lapply(1L:length(term.labels),
                         function(whichterm) (1:ncol(covariates))[assign==whichterm])
  ccs.by.term <-
    lapply(cols.by.term, function(whichcols) {
      complete.cases(covariates[,whichcols])
    })
  names(ccs.by.term) <- term.labels
  null.record <- rowSums(as.data.frame(ccs.by.term))==0
  if (any(null.record)) #in this case make sure it's NA in each covariate column
      null.record[null.record] <- apply(is.na(covariates[null.record,,drop=FALSE]), 1, all)

  terms.with.missings <- !sapply(ccs.by.term, all)

  ccs.by.term <- c(list('_non-null record_'=!null.record),
                   ccs.by.term)
  terms.with.missings <- c(TRUE, # in order always to have same leading entry
                           terms.with.missings)

  nm.covs <- integer(ncol(covariates))
  nm.terms <- integer(length(terms.with.missings))

  nmcols <- ccs.by.term[terms.with.missings]
  nmcols.dupes <- duplicated(nmcols, fromLast=FALSE)
  if (any(nmcols.dupes))
  {
      nm.terms[terms.with.missings][!nmcols.dupes] <- 1L:sum(!nmcols.dupes)
      nm.terms[terms.with.missings][nmcols.dupes] <-
        match(nmcols[nmcols.dupes], nmcols[!nmcols.dupes])
      ## At this point nm.terms is a look-up table with an entry for
      ## _non-null record_ and for each provided term. Entries are
      ## 0 if there are no null records, or no missings on the term in question;
      ## otherwise if there are null records and/or term missings, the
      ## entry gives the column of the matrix `notmissing`, formed a few lines
      ## down, that describes the relevant missingness pattern.
  } else nm.terms[terms.with.missings] <- 1L:length(nmcols)

  notmissing <- as.data.frame(nmcols[!nmcols.dupes])
  names(notmissing) <- names(ccs.by.term)[terms.with.missings][!nmcols.dupes]
  ## Now form the look-up table associating columns of the covariates matrix
  ## with columns of matrix `notmissing` describing relevant missingness patterns.
  ## Columns associated with terms on which there was no missingness get a 0 here.
  nm.terms <- nm.terms[-1L]
  nm.covs[assign>0] <- nm.terms[assign[assign>0]]

  notmissing <- as.matrix(notmissing)

  ## MMF: I'm getting a weird bug that won't let constrasts be "NULL"
  if (is.null(contrasts)) {
      contrasts <- list()
  }

  new("ModelMatrixPlus", covariates,
      OriginalVariables = assign,
      TermLabels        = term.labels,
      Contrasts         = contrasts,
      NotMissing        = notmissing,
      NM.Covariates     = nm.covs,
      NM.terms          = nm.terms,
      UnitWeights       = uweights )
}
