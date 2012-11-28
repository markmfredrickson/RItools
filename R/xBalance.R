xBalance <- function(fmla, strata=list(unstrat=NULL),
                     data,
                     report=c("all", "std.diffs","z.scores","adj.means","adj.mean.diffs","adj.mean.diffs.null.sd",
                       "chisquare.test","p.values")[1],
#                     include.means=FALSE, chisquare.test=FALSE,
                     stratum.weights=harmonic, na.rm=FALSE,
                     covariate.scaling=NULL, normalize.weights=TRUE,impfn=median,
                     groups = NULL)
{
  stopifnot(class(fmla)=="formula",
            all(report %in% c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test",
                              "std.diffs","z.scores","p.values","all")),
            is.null(strata) || is.factor(strata) || is.list(strata),
            !is.data.frame(strata) || !any(is.na(names(strata))),
            !is.data.frame(strata) || all(names(strata)!=""),
            !is.data.frame(strata) || all(sapply(strata, is.factor)),
            is.null(data) || is.data.frame(data)
            )
  if (is.null(strata)) warning("Passing NULL as a 'strata=' argument is depracated;\n for balance w/o stratification pass 'list(nostrat=NULL)' instead.\n (Or did you mean to pass a non-NULL 'strata=' argument? Then check for typos.)")
  if (is.list(strata) && !is.data.frame(strata) &&
      !all(sapply(strata, function(x) (is.null(x) | inherits(x,"formula"))))
      ) stop("For balance against multiple alternative stratifications,\n please make 'strata' either a data frame or a list containing formulas or NULL entries.")
  if("all" %in% report){report<-c("adj.means","adj.mean.diffs","adj.mean.diffs.null.sd","chisquare.test",
                              "std.diffs","z.scores","p.values")}

### NA Handling ##  
  if (na.rm==TRUE)
    {
      tfmla <- terms.formula(fmla,dat=data, keep.order=TRUE) 
    } else
  {
    data <- naImpute(fmla,data,impfn)
    tfmla <- attr(data, 'terms')
  }
### End NA handling ###
  

###Extract the treatment var
  if (!attr(tfmla, "response")>0) 
    stop("fmla must specify a treatment group variable")

  zz <- eval(tfmla[[2]], data, parent.frame()) # changed for v.93, see comment in log
  zzname<- deparse(tfmla[[2]])
  if (!is.numeric(zz) & !is.logical(zz))
    stop("LHS of fmla should be logical or numeric")
  if (any(is.na(zz))) stop('NAs on LHS of fmla not allowed.')
### End extract treatment var

## Investigate the grouping variables

all.variables <- attr(tfmla, "term.labels")

if (is.null(groups)) {
  groups <- list("All" = all.variables)
} else {

  if (any(0 == sapply(groups, length))) {
    stop("Empty groups not permitted.")
  }

  vars.in.groups <- unique(unlist(groups))
  unknown.vars <- !(vars.in.groups %in% all.variables)

  if (zzname %in% vars.in.groups) {
    stop("Treatment variable is not permitted in groups.")
  }

  if (any(unknown.vars)) {
    stop(paste("Unknown variable(s):", paste(vars.in.groups[unknown.vars], collapse = ", ")))
  }

  groups <- append(list("All" = all.variables), groups)
}
  
# NB: I've tried without the explicit model.frame call, but weird errors would pop up)
mf <- model.frame(tfmla, data = data, na.action = na.pass)
mm1 <- make_nice_model_matrix(tfmla, mf)

# the groups variable is a list of list of variables
# this next few lines converts the variable names into the column indices of the variables in the model matrix
# a categorial variable, for example, should expand into several indices, one for each dummy
# likewise for splines or other unusual expansions

group.indices <- lapply(groups, function(x) {
  # x should be a character vector.
  return(attr(mm1, "assignnames") %in% x)
})


### Prepare ss.df, data frame of strata
  if (is.null(strata)) ss.df <- data.frame(unstrat=factor(numeric(length(zz))))
      
  if (is.factor(strata) & length(strata)!=length(zz)) stop("length of strata doesn\'t match dim of data")

  if (is.factor(strata)) ss.df <- data.frame(strat=factor(strata))
  if (is.data.frame(strata)) ss.df <- as.data.frame(lapply(strata,factor))
  if (is.list(strata) & !is.data.frame(strata))
### In this case strata should be a list of formulas    
    {
      pfr <- parent.frame()
    ss.df <-
      lapply(strata,
             function(fmla) {
               if (is.null(fmla)) factor(numeric(length(zz))) else {
               ss <- eval(attr(terms(fmla), "variables"), data, 
                          pfr) 
               if (length(ss)-1) interaction(ss, drop=TRUE) else factor(ss[[1]])
             }
             }
             )
    ss.df <- as.data.frame(ss.df)
  }
### End prepare ss.df, data frame of strata

### Remove stratification variables without levels (e.g., all NAs), with a warning
if (any(ss.rm <- !sapply(ss.df, nlevels)))
  {
    if (length(ss.df)==1) stop("'strata=' variable contains no strata.  Perhaps it evaluates to NAs?")
    if (all(ss.rm)) stop("'strata=' variables contain no strata.  Perhaps they all evaluate to NAs?")
    ss.rm.nms <- if (is.null(names(ss.df))) which(ss.rm) else names(ss.df)[ss.rm]
    ss.rm.nms <- paste(ss.rm.nms, collapse=" ,")
    warning(paste("Removing the following strata entries, which contained no strata.\n(Perhaps they evaluate to NAs?)\n",
            ss.rm.nms)
            )
    ss.df <- ss.df[!ss.rm]
}
### End remove stratification variables without levels  
  gs.df <- xBalance.find.goodstrats(ss.df,zz,mm1)
  
  swt.ls <- xBalance.make.stratwts(stratum.weights,ss.df, gs.df, zz, data, normalize.weights)

  s.p <- if (is.null(covariate.scaling)) {xBalance.makepooledsd(zz,mm1,dim(mm1)[1])
                                        } else 1

### Call xBalanceEngine here.
  
  RES <- lapply(names(ss.df),
                function(nm) {
###                  workingswt.ls<-swt.ls[[nm]]  # shouldn't be neccessary after r216 change to xBalance.make.stratwts
###                  workingswt.ls[["wtratio"]]<-swt.ls[[nm]][["wtratio"]][gs.df[[nm]]]
                  xBalanceEngine(factor(ss.df[gs.df[[nm]],nm]),
                                            zz[gs.df[[nm]]],
                                            mm1[gs.df[[nm]],,drop=FALSE],
                                            report, swt.ls[[nm]], 
                                            s.p, normalize.weights, zzname, group.indices)
                                            
                            }
                )
  names(RES) <- names(ss.df)
  ##nms <- paste(rep(names(ss.df), rep(length(RES[[1]]$dfr),length(ss.df))),
  ##            names(RES[[1]]$dfr), sep=".")
  ans <- list() ##the overall function still returns a list because of the overall test info.
  ##results is an array of variables by balance statistics by stratification.
  ##here assuming that the variables and statistics are the same across stratifications (including unstratified).
  ans$results<-array(dim = c(vars = nrow(RES[[1]][["dfr"]]),
                             stat = ncol(RES[[1]][["dfr"]]),
                             strata = length(RES)),
                     dimnames = list(vars = rownames(RES[[1]][["dfr"]]),
                                     stat = colnames(RES[[1]][["dfr"]]),
                                     strata = names(RES)))

  # groups is also an array, this time the first index is the strata, the 2nd are the chisquare test info, 3rd are any groups
  ans$groups <- array(dim = c(length(RES), # number of strata
                               3, # chi.squared stat, deg. freedom, p.value
                               length(groups)),
                      dimnames = list(strata  = names(RES), 
                                      tests = c("chisquare", "df", "p.value"), 
                                      groups = names(groups)))

  # we populate both arrays by iterating through the strata
  for(i in names(RES)) {
    tmp <- RES[[i]]
    ans$results[,,i] <- as.matrix(tmp[["dfr"]])

    for (j in names(groups)) {
      xxx <- tmp$groups[[j]]
      ans$groups[i, 'chisquare', j] <- xxx$chisquare
      ans$groups[i, 'df', j]        <- xxx$df
      ans$groups[i, 'p.value', j]   <- pchisq(xxx$chisquare,
                                              df = xxx$df,
                                              lower.tail = FALSE)
    }

  }
  ans$overall <- as.data.frame(ans$groups[,,1, drop = F]) # for backwared compatiblity, will probably be removed in future versions
  colnames(ans$overall) <- c("chisquare", "df", "p.value")

  attr(ans, "fmla") <- formula(tfmla)

  class(ans) <- c("xbal", "list")
  ans
}

xBalance.make.stratum.mean.matrix <- function(ss, mm) {

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

# Creates a "nice" model matrix that has the following properties:
# - all logicals and 2 level factors are just convert to 0/1 numerics
# - all factors with 3 or more levels are given "X = a", "X = b" labels
# - no intercept is included.
# It otherwise returns a standard model matrix with an additional "assignnames" property 
# that maps the usual "assign" attribute to the actual lables of the variables, not just their
# indices.
#
# @param tf A "terms" object, as applied to a formula and some data
# @param mf A model.frame() object. Requiring it to be passed so that 
# @return A model.matrix like object.
make_nice_model_matrix <- function(tf, mf) {

  # do a little "pre-processing" of the data so that two level stuff just gets single item
  mf <- lapply(mf, function(x) {
    if (is.logical(x)) {
      return(as.numeric(x))
    }
    if (is.factor(x) && nlevels(x) == 2) {
      return(as.numeric(x))
    }
    
    # don't change it
    return(x)
  })

  mm <- model.matrix(tf, mf, contrasts.arg = lapply(Filter(is.factor, mf), make_nice_contrasts)) 
  oldassign <- attr(mm, "assign")

  mm <- mm[, -1, drop = F] # drop the intercept
  attr(mm, "assign") <- oldassign[-1] # also drop intercept
  attr(mm, "assignnames") <- attr(tf, "term.labels")[attr(mm, "assign")]

  return(mm)
}

make_nice_contrasts <- function(x) {
  structure(diag(nlevels(x)),
            dimnames = list(paste(" =", levels(x)),
                            paste(" =", levels(x))))
}
