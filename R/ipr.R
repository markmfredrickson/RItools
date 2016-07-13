##' Assemble inverse probability of assignment weights on the fly,
##' e.g. as part of evaluation of `weights` argument to a `glm` or similar
##' function.   Assumes a designs that is completely randomized, either by cluster or by
##' the unit of observation, potentially within strata or blocks.
##'
##' The function assumes complete random assignment, within (optional)
##' blocks.  This can be complete random assignment of clusters; in this situation
##' the user avoids having to manage aggregation to the cluster level in order
##' properly to infer treatment probabilities.  There need to be at least two conditions,
##' and any block without all conditions represented will be dropped, as are observations
##' that are NA on the blocking variable.
##'
##' @title Inverse probability weighting suitable for completely randomized designs
##' @param z categorical variable recording distinctions among assignments to treatment conditions
##' @param strat categorical variable recording blocks
##' @param clus optional categorical variable indicating cluster membership
##' @return vector of inverse probability of assignment weights
##' @author Bendek B. Hansen
ipr <- function(z,strata,clus=NULL)
  {
            stopifnot(length(z)==length(strata), is.null(clus) | length(clus)==length(strata),
                                        !all(z==z[1], na.rm=TRUE))
                    strata <- as.factor(strata)
                    if (is.ordered(z)) warning("I received an ordinal z. I'll treat it the same as any other factor (FYI).")

                    if (is.character(z)) z <- as.factor(z)
                    z <- as.numeric(z)

                    if (!is.null(clus))
                      {
                        ## input checking
                        ## each cluster should have 1 and only 1 type of unit (treated or control)
                        tbl <- table(z, clus)
                        isGood <- apply(tbl, 2, function(x) { sum(x != 0)  == 1})
                        if (!(all(isGood))) {
                          stop("The following clusters did not have identical treatment assignment: ",
                                                                                                            paste(collapse = ", ", colnames(tbl)[!isGood]))
                                                                            }
                        ## each cluster should be nested entirely within a single stratum
                        tbl <- table(strata, clus, useNA = "ifany")
                        isGood <- apply(tbl, 2, function(x) { sum(x != 0) == 1 })
                        if (!(all(isGood))) {
                          stop("In ", s, ", the following clusters were not nested within stratum levels: ",
                               paste(collapse = ", ", colnames(tbl)[!isGood]))
                        }

                        cluster.representatives <- !duplicated(clus)
                        z <- z[cluster.representatives]
                        strat.c <- strata[cluster.representatives,
                                         drop=FALSE] # Can there be levels that
                        ## don't get represented? Not sure but if so then this is safer
                      } else
            strat.c <- strata


            strat.by.z <- table(strat.c, z)
            strat.all <- apply(strat.by.z, 1, # are *all* levels of z represented in
                               all) # the stratum? If not we're going to exclude it
            strat.tot <- apply(strat.by.z, 1, sum)
            HT.by.strat <- ifelse(rep(strat.all, ncol(strat.by.z)),
                                  strat.tot/strat.by.z, 0)
            dim(HT.by.strat) <- dim(strat.by.z)
            
            ans <- numeric(length(strata))
            for (zz in 1L:ncol(strat.by.z))
              {
                zval <- as.numeric(colnames(strat.by.z)[zz])
                htvals <- HT.by.strat[, zz, drop=TRUE]
                toreplace <- z==zval
                toreplace[is.na(z)] <- FALSE
                ans[toreplace] <- htvals[strat[toreplace]]
              }
            ans
          }
