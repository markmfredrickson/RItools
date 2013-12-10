
xtable.xbal <- function(x,caption = NULL, label = NULL, align =c("l",rep("r",ncol(xvardf))),
                        digits = 2, display = NULL, col.labels=NULL, ...) {
  ##By default use decimal alignment, which will require the dcolumn package in latex and an appropriate column definition like:
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

  stopifnot(require(xtable))
  xprint <- flatten.xbalresult(x)
  numstrata<-dim(x$results)[3]
  latex.annotation <- attr(xprint, "latex.annotation")
  xvardf<-xprint$vartable

  if (!is.null(col.labels))
    names(xvardf) <- col.labels


  ##call xtable on the resulting data.frame
  vartab <- xtable(xvardf,caption=caption, label=label, digits=digits,align=align,display=display,col.labels=col.labels,...) ##NextMethod("xtable",xvardf)
  structure(vartab,
            latex.add.to.row=list(pos=list(-1),command=latex.annotation),
            hline.after=c(0,nrow(xvardf)))

}
