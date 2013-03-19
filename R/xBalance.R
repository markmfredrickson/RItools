xBalance <- function(fmla, strata=list(unstrat=NULL),
                     data,
                     report=c("std.diffs","z.scores","adj.means","adj.mean.diffs","adj.mean.diffs.null.sd",
                       "chisquare.test","p.values", "all")[1:2],
#                     include.means=FALSE, chisquare.test=FALSE,
                     stratum.weights=harmonic, na.rm=FALSE,
                     covariate.scaling=NULL, normalize.weights=TRUE,impfn=median)
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
  
mm1 <- xBalance.makeMM(tfmla,data)

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
                                            s.p, normalize.weights,zzname)
                            }
                )
  names(RES) <- names(ss.df)
  ##nms <- paste(rep(names(ss.df), rep(length(RES[[1]]$dfr),length(ss.df))),
  ##            names(RES[[1]]$dfr), sep=".")
  ans <- list() ##the overall function still returns a list because of the overall test info.
  ##results is an array of variables by balance statistics by stratification.
  ##here assuming that the variables and statistics are the same across stratifications (including unstratified).
  ans$results<-array(dim=c(vars=nrow(RES[[1]][["dfr"]]),stat=ncol(RES[[1]][["dfr"]]),strata=length(RES)),
                     dimnames=list(vars=rownames(RES[[1]][["dfr"]]),stat=colnames(RES[[1]][["dfr"]]),strata=names(RES)))
  for(i in names(RES)){
    ##print(i);print(RES[[i]][["dfr"]])
    ans$results[,,i]<-as.matrix(RES[[i]][["dfr"]])
  }
  ##dimnames(ans)[["stat"]][grep("Tx",dimnames(ans)[["stat"]])]<-c("adj.mean.strata=0","adj.mean.strata=1")
  ##ans$by.variable <- do.call(cbind, lapply(RES, function(x) x[['dfr']]) )
  ##colnames(ans$by.variable) <- nms
  attr(ans, "fmla") <- formula(tfmla)

  if ("chisquare.test" %in% report)
  {
    ans$overall <- data.frame(chisquare=numeric(length(RES)),
                              df=numeric(length(RES)),
                              p.value=numeric(length(RES)),
                              row.names=names(RES))
  for (nn in names(RES))
  {
   ans$overall[nn,'chisquare'] <- RES[[nn]]$chisq['chisquare']
   ans$overall[nn,'df'] <- RES[[nn]]$chisq['df']
   ans$overall[nn,'p.value'] <- pchisq(RES[[nn]]$chisq['chisquare'],
                                       df=RES[[nn]]$chisq['df'],
                                       lower=FALSE)
  }
}
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



