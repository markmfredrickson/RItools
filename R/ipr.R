##' Assemble inverse probability of assignment weights on the fly,
##' e.g. as part of evaluation of `weights` argument to a `glm` or similar
##' function.   Also calculates inverse odds of assignment weights.
##' Assumes a design that is completely randomized, either by cluster or by
##' the unit of observation, potentially within strata or blocks.
##'
##' For inverse probability of assignment weights, the default, set argument
##' \code{type} to \dQuote{inverse probability}.  Inverse odds of assignment
##' calls for specification of a reference level, with observations in that
##' condition receiving weights of 1. State this reference level within your
##' \code{type} argument: for example, if your treatment condition corresponds
##' to \code{z==1}, and you want inverse odds of treatment weighting, then pass
##' a \code{type} argument of \dQuote{odds vs 1} or \dQuote{odds against 1}. 
##'
##' The function assumes complete random assignment, within (optional)
##' blocks.  This can be complete random assignment of clusters; in this situation
##' the user avoids having to manage aggregation to the cluster level in order
##' properly to infer treatment probabilities.  There need to be at least two conditions,
##' and any block without all conditions represented will be dropped, ie given observation 
##' weights of 0, as are observations that are NA on the blocking variable.
##'
##' NAs in the assignment variable also translate to zero.  If they occur within a
##' non-NA block, probabilities or odds for that block are calculated as though the 
##' observations with NAs for \code{z} just weren't there. 
##'
##' @title Infer inverse probability/odds of assignment weights in randomized designs
##' @param z variable recording distinctions among assignments to treatment conditions
##' @param strata categorical variable recording blocks
##' @param clus optional categorical variable indicating cluster membership
##' @param type character string naming desired type of weight. See Details for options other than inverse probability
##' @return vector of weights
##' @export
##' @author Ben B. Hansen
ipr <- function(z,strata,clus=NULL, type="inverse probability")
  {
            stopifnot(length(z)==length(strata), is.null(clus) | length(clus)==length(strata),
                      !all(is.na(z)), !all(z==z[which.min(is.na(z))], na.rm=TRUE),
                      is.character(type), substr(type, 1, 4) %in% c('inve', 'odds'))

            strata <- as.factor(strata)
            if (is.ordered(z)) warning("I received an ordinal z. I'll treat it the same as any other factor (FYI).")

            z <- as.factor(z)
            
            if (tolower(substr(type, 1, 4))=='odds')
                {
            reflev <-
                getRefLev(type[1], levels(z),
                              errmsg=paste0('For odds please specify reference category, e.g.type="odds vs ',
                                  levels(z)[nlevels(z)], '"')
                              )
            type <- "odds"
        }
                    
            
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
            strat.keep <- apply(strat.by.z, 1, # are *all* levels of z represented in
                               all) # the stratum? If not we're going to exclude it
            wt.by.strat <- if (type=="odds")
                {
                    strat.reflev <- strat.by.z[,reflev,drop=TRUE]
                    ifelse(rep(strat.keep, ncol(strat.by.z)),
                           strat.reflev/strat.by.z, 0)
                } else 
                    {
            strat.tot <- apply(strat.by.z, 1, sum)
            ifelse(rep(strat.keep, ncol(strat.by.z)),
                                  strat.tot/strat.by.z, 0)
        }
            dim(wt.by.strat) <- dim(strat.by.z)

            ans <- numeric(length(strata))
            for (zz in 1L:nlevels(z))
              {
                zval <- levels(z)[zz]
                wtvals <- wt.by.strat[, zz, drop=TRUE]
                toreplace <- z==zval
                toreplace[is.na(z)] <- FALSE
                toreplace[is.na(strata)] <- FALSE
                ans[toreplace] <- wtvals[strata[toreplace]]
              }
            ans
          }

getRefLev <- function(typestring, tab, errmsg="I can't infer the desired reference level")
    {
        stopifnot(is.character(typestring), length(typestring)==1,
                  is.character(tab))
        if (tolower(substr(typestring, 1L, 4L))!="odds") stop(errmsg)
        nc <- nchar(typestring)
        if (nc==4) stop(errmsg)

        startfrom <- 5L
        if (substr(typestring,5L,6L)=='vs' & nc>6) startfrom <- 7L
        if (substr(typestring,5L,7L)==' vs' & nc>7) startfrom <- 8L
        if (substr(typestring,5L,7L)=='.vs' & nc>7) startfrom <- 8L
        if (substr(typestring,5L,8L)==' vs ' & nc>8) startfrom <- 9L
        if (substr(typestring,5L,8L)=='.vs ' & nc>8) startfrom <- 9L
        if (substr(typestring,5L,8L)=='.vs.' & nc>8) startfrom <- 9L

        reflevspec <- substr(typestring,startfrom,nchar(typestring))

                            
        strReverse <- function(x)
            sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")

        rlevs <- strReverse(tab)
        rtype <- strReverse(reflevspec)

        matches <- match(rlevs,rtype, nomatch=0)
        if (!any(matches))
            {
                matches <- charmatch(rlevs,rtype, nomatch=-1)
                if (all(matches<0)) stop(errmsg)
                if (!any(matches==1))
                    stop(paste("I can't determine which of these reference levels you want:\n",
                               paste(tab[matches==0], collapse=" ")) )
                matches[matches<0] <- 0
            }
        ## If we've come this far then there's at least one candidate for a match
        matches <- as.logical(matches)
        
        longestmatch <- max(nchar(rlevs)[matches])
        if (sum(nchar(rlevs)[matches]==longestmatch)>1)
            stop(paste("I can't determine which of these reference levels you want:\n",
                       paste(tab[matches], collapse=" "))
                 )
        tab[matches][nchar(rlevs)[matches]==longestmatch]
    }
