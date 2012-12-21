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

  nvars <- dim(theresults)[1]
  nstrata <- dim(theresults)[2]
  varlabels <- rownames(theresults)

  ypos <- seq(length.out = nvars)
  xrange <- range(theresults, na.rm = TRUE)
  xrange <- xrange + xrange * adjustxaxis

  ##Setup the margin, adjust for the lengths of the labels.
  par(mar=mar,tck=tck,mgp=mgp) ##set default margins
  mymai<-par('mai')
  mymai[2]<-max(c(strwidth(varlabels ,units="inches"),mymai[2]))
  ##Setup the plotting region
  par(mai=mymai)

  if(length(thecols)!=length(which.strata)){
    if(length(thecols)==1){
      thecols<-rep(thecols,length(which.strata))
    }
    if(length(thecols)>1){
      cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
    }
  }
  names(thecols)<-which.strata

  if(length(thesymbols)!=length(which.strata)){
    if(length(thesymbols)==1){
      thesymbols<-rep(thesymbols,length(which.strata))
    }
    if(length(thesymbols)>1){
      cat("I dont know which colors belong with which strata. Please either provide a vector of columns the same length as the number of stratifications or a single color to be used for all stratifications. \n"); stop()
    }
  }
  names(thesymbols)<-which.strata
  
  
  plot(xrange,range(ypos),axes=FALSE,pch=19,col="blue",
       ylab="",xlab=thexlab,type="n",...)
  for(i in which.strata){
    points(theresults[,i],ypos,col=thecols[i],pch=thesymbols[i])
    }
  if(segments&length(which.strata)>1){ ##segments are mainly useful for drawing the eye to changes in balance along a single variable across more than 1 stratification
    for(j in ypos){
      segments(min(theresults[j,]),j,
               max(theresults[j,]),j,col=gray(.7),lwd=.5) 
    }
  }
  axis(1,at=pretty(seq(xrange[1],xrange[2],length=5)))
  axis(2,labels=varlabels,at=ypos,las=2,tick=FALSE)
  lines(c(0,0),range(ypos)+c(-.025*length(ypos),.025*length(ypos)),col="grey",lwd=1)
  ##segments(0,min(ypos),0,max(ypos),col="grey",lwd=1)
  if(legend){
    legend(x="topright",#xrange[1],ypos[ypos==max(ypos)],
           legend=thestratalabs,
           col=thecols,
           pch=thesymbols,
           bty="n")
  }
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
#' To ease display, variables can be grouped and the results plotted together.
#' When groups are present, variables are sorted within groups from largest to
#' smallest imbalance.
#'
#' @param x A matrix of variables (rows) by stratifications (columns).
#' @param groups An optional factor of variable groups.
#' @param ordered Should the variables be ordered (within groups if any) from most to least imbalance on the first statistic?
#' @param xlab The label of the x-axis of the plot.
#' @param ... Additional arguments to pass to \code{\link{plot.default}}.
#' @seealso \code{\link{plot.xbal}} \code{\link{xBalance}}
#' @example inst/examples/balanceplot.R
#' @export
balanceplot <- function(x, groups = NULL, ordered = F, xlab = "Balance", ...) {

  original.par <- par()

  nvars <- dim(x)[1]
  nstrat <- dim(x)[2]


  xrange <- range(x, na.rm = TRUE)
  xrange <- xrange + xrange * 0.25
  
  
  if (is.null(groups)) {
    groups <- rep(1, nvars)
  }
  
  groups <- as.factor(groups)

  if (ordered) {
    # order X by the groups, and within groups order by the first column
    localorder <- order(groups, x[,1])
    x <- x[localorder, , drop = F]
    groups <- groups[localorder]
  }

  grpnames <- levels(groups)
  pergroup <- table(groups)[grpnames]

  ypos <- 1:nvars + unlist(mapply(rep, 0:(length(grpnames) - 1), pergroup))

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

  if (!is.null(groups)) {
    text(mean(xrange), cumsum(pergroup) + (1:length(grpnames)), labels = grpnames)
  }

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
