
xtable.xbal <- function(x,caption = NULL, label = NULL, align =c("l",rep("r",ncol(xdf))),
                          digits = 2, display = NULL, col.labels = NULL,
                          type = "variables", ...)
  {##By default use decimal alignment, which will require the dcolumn package in latex and an appropriate column definition like:
    ##\newcolumntype{.}{D{.}{.}{2.2}}
    ##Here is an example which works
    ##xb1<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
    ##         strata=data.frame(unstrat=factor(character(32)),
    ##           pt=factor(nuclearplants$pt)),
    ##         data=nuclearplants,
    ##         report=c("adj.means","adj.mean.diffs",'std.diffs', 'z.scores', 'chisquare.test','p.values'))
    ##
    ##junk<-xtable(xb1)
    ##print(junk,add.to.row=attr(junk,"latex.add.to.row"),hline.after=c(0,nrow(junk)),sanitize.text.function=function(x){x},floating=TRUE,floating.environment="sidewaystable")

  if (!any(type %in% c("variables", "groups"))) {
    stop("Unknown type of table: ", type)
  }
  
  stopifnot(require(xtable))
 
  attrs <- list()

  if (type == "variables") {
  
    xprint <- flatten.xbalresult(x)
    numstrata<-dim(x$results)[3]
    latex.annotation <- attr(xprint, "latex.annotation")
    xdf<-xprint$vartable
  
    if(!is.null(col.labels)) names(xdf) <- col.labels
    
    attrs <- list(latex.add.to.row = 
                    list(pos = list(-1),
                         command = attr(xprint, "latex.annotation")),
                  hline.after = c(0, nrow(xdf)))
  }


  if (type == "groups") {

    tmp <- ftable(x$groups, 
                  row.vars = "groups", 
                  col.vars = c("strata", "tests"))

    xdf <- tmp[,] # a way to cast back to just a simple matrix

    rownames(xdf) <- dimnames(x$groups)$groups
    colnames(xdf) <- rep(dimnames(x$groups)$tests, 2)

  }

  xt <- xtable(xdf, 
               caption = caption, 
               label = label, 
               digits = digits,
               align = align,
               display = display,
               col.labels = col.labels,
               ...) 

  return(do.call(structure, append(list(xt), attrs)))
  
}


