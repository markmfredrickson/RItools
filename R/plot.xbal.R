#' Plot of balance across multiple stratification schemes
#'
#' The plot allows a quick visual comparison of the effect of different
#' stratification designs on the comparability of different
#' variables. This is not a replacement for the omnibus statistical test
#' reported as part of \code{\link{print.xbal}}. This plot does allow the
#' analyst an easy way to identify variables that might be the primary culprits
#' of overall imbalances and/or a way to assess whether certain important
#' covariates might be imbalanced even if the omnibus test reports that
#' the stratification overall produces balance.
#'
#' By default all variables and all stratifications are plotted. The scope
#' of the plot can be reduced by using the \code{\link{subset.xbal}} function to
#' make a smaller \code{xbal} object with only the desired variables or
#' stratifications.
#'
#' @param x An object returned by \code{\link{xBalance}}
#' @param xlab The label for the x-axis of the plot
#' @param statistic The statistic to plot. The default choice of standardized
#' difference is a good choice as it will have roughly the same scale for all
#' plotted variables.
#' @param thecols NOT YET IMPLEMENTED
#' @param thesymbols NOT YET IMPLEMENTED
#' @param absolute Convert the results to the absolute value of the statistic.
#' @param strata.labels A named vector of the from \code{c(strata1 = "Strata Label 1", ...)} 
#' that maps the stratification schemes to textual labels.
#' @param variable.labels A named vector of the from \code{c(var1 = "Var Label1", ...)} 
#' that maps the variables to textual labels.
#' @param ... additional arugments to pass to \code{\link{balanceplot}}
#' @seealso \code{\link{xBalance}} \code{\link{subset.xbal}} \code{\link{balanceplot}}
#' @example inst/examples/plot.xbal.R
#' @import abind
plot.xbal<-function(x,adjustxaxis=.25,segments=TRUE,legend=TRUE,
                    mar=c(3,3,2,0)+0.1,mgp=c(1.5,.5,0),tck=-.01,
                    xlab = "Standardized Differences",
                    statistic = "std.diff",
                    thecols=rainbow(length(which.strata)),
                    thesymbols=c(19,22,23,24,25)[1:length(which.strata)],
                    absolute = FALSE,
                    strata.labels = NULL,
                    variable.labels = NULL,
                    ...){

  # the helper .plot.xbal turns xb and the many of the arguments into a
  # an ordered variables by stratas table. The colnames are the labeled stratifications
  # the rownames are the labeled variables
  if (dim(x$results)[2] > 1) {
    # this means that the user is passing an xBalance object with more than one statistic
    # so we need to trim it down
    
    # but first we need to make sure the statistic exists
    if (!(statistic %in% dimnames(x$results)[[2]])) {
      stop("Unknown statistic: ", statistic)
    }
    x <- subset(x, stats = statistic)
  }

  x <- adrop(x$results, drop = 2)  

  if (!is.null(variable.labels)) {
    if (is.null(names(variable.labels))) {
      stop("Variable labels must be a named vector of the form c('var1' = 'Var One', ...)")
    }
    rownames(x) <- variable.labels[rownames(x)]
  }

  if (!is.null(strata.labels)) {
    if (is.null(names(strata.labels))) {
      stop("Strata labels must be a named vector of the form c('var1' = 'Var One', ...)")
    }
    colnames(x) <- strata.labels[colnames(x)]
  }

  if (absolute) {
    x <- abs(x)
  }


  return(balanceplot(x, xlab = xlab, ...))

  ### NOT RUN: (but saving while we transition to the more general balanceplot function 

  # nvars <- dim(theresults)[1]
  # nstrata <- dim(theresults)[2]
  # varlabels <- rownames(theresults)

  # ypos <- seq(length.out = nvars)
  # xrange <- range(theresults, na.rm = TRUE)
  # xrange <- xrange + xrange * adjustxaxis

  # ##Setup the margin, adjust for the lengths of the labels.
  # par(mar=mar,tck=tck,mgp=mgp) ##set default margins
  # mymai<-par('mai')

  # # when using the SVG device, strwidth throws fits, so we hope that the default mai[2] is good enough
  # if (names(dev.cur()) != "svg") {
  #   mymai[2]<-max(c(strwidth(varlabels ,units="inches"),mymai[2]))
  # }

  # ##Setup the plotting region
  # par(mai=mymai)

  # if(length(thecols)!=length(which.strata)){
  #   if(length(thecols)==1){
  #     thecols<-rep(thecols,length(which.strata))
  #   }
  #   if(length(thecols)>1){
  #     cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
  #   }
  # }
  # names(thecols)<-which.strata

  # if(length(thesymbols)!=length(which.strata)){
  #   if(length(thesymbols)==1){
  #     thesymbols<-rep(thesymbols,length(which.strata))
  #   }
  #   if(length(thesymbols)>1){
  #     cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
  #   }
  # }
  # names(thesymbols)<-which.strata
  # 
  # 
  # plot(xrange,range(ypos),axes=FALSE,pch=19,col="blue",
  #      ylab="",xlab=thexlab,type="n",...)
  # for(i in which.strata){
  #   points(theresults[,i],ypos,col=thecols[i],pch=thesymbols[i])
  #   }
  # if(segments&length(which.strata)>1){ ##segments are mainly useful for drawing the eye to changes in balance along a single variable across more than 1 stratification
  #   for(j in ypos){
  #     segments(min(theresults[j,]),j,
  #              max(theresults[j,]),j,col=gray(.7),lwd=.5) 
  #   }
  # }
  # axis(1,at=pretty(seq(xrange[1],xrange[2],length=5)))
  # axis(2,labels=varlabels,at=ypos,las=2,tick=FALSE)
  # lines(c(0,0),range(ypos)+c(-.025*length(ypos),.025*length(ypos)),col="grey",lwd=1)
  # ##segments(0,min(ypos),0,max(ypos),col="grey",lwd=1)
  # if(legend){
  #   legend(x="topright",#xrange[1],ypos[ypos==max(ypos)],
  #          legend=thestratalabs,
  #          col=thecols,
  #          pch=thesymbols,
  #          bty="n")
  # }
}

#' Create a plot of the balance on variables across different stratifications.
#'
#' This plotting function summarizes variable by stratification matrices. For
#' each variable (a row in the \code{x} argument), the values are under each
#' stratification (the columns of \code{x}) plotted on the same line. 
#' 
#' It is
#' conventional to standardize the differences to common scale (e.g. z-scores),
#' but this is not required. Plotting will automatically order the data from
#' largest imbalance to smallest based on the first column of \code{x}.
#' 
#' @param x A matrix of variables (rows) by stratifications (columns).
#' @param ordered Should the variables be ordered (within groups if any) from most to least imbalance on the first statistic?
#' @param xlab The label of the x-axis of the plot.
#' @param ... Additional arguments to pass to \code{\link{plot.default}}.
#' @seealso \code{\link{plot.xbal}} \code{\link{xBalance}}
#' @example inst/examples/balanceplot.R
#' @export
balanceplot <- function(x, ordered = F, xlab = "Balance", ...) {
  original.par <- par()

  nvars <- dim(x)[1]
  nstrat <- dim(x)[2]

  xrange <- range(x, na.rm = TRUE)
  xrange <- xrange + xrange * 0.25
  
  if (ordered) {
    # order X by the groups, and within groups order by the first column
    localorder <- order( x[,1])
    x <- x[localorder, , drop = F]
  }

  ypos <- 1:nvars

  mai <- par('mai')
  mai[2] <- max(strwidth(rownames(x), units = "inches")) + mai[2]
  mar <- par('mar')
  mar[3] <- nstrat + 2
  original.par <- par(mar = mar, mai = mai)
 
  plot(xrange, 
       range(ypos) + c(0,1) ,
       axes = FALSE,
       pch = 19,
       col = "blue",
       ylab = "",
       xlab = xlab,
       type = "n",
       ...)

  for(i in 1:nstrat) {
    points(x[,i], ypos, pch = i, cex = 0.5) # col =thecols[i],pch=thesymbols[i])
  }

  axis(1, at = pretty(seq(xrange[1], xrange[2], length = 5)))
  axis(2, labels = rownames(x), at = ypos, las = 2, tick = FALSE)
  abline(v = 0, col = "#333333")


  legend(x = mean(xrange),
         y = max(ypos) + 2,
         legend = colnames(x),
         pch = 1:nstrat,
         bty = "n",
         xpd = T,
         xjust = 0.5,
         yjust = 0)

  par(original.par)
} 
