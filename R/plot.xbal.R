plot.xbal<-function(x,adjustxaxis=.25,segments=TRUE,legend=TRUE,
                    mar=c(3,3,2,0)+0.1,mgp=c(1.5,.5,0),tck=-.01,
                    which.strata=dimnames(x$results)[["strata"]],thestratalabs=which.strata,
                    which.stat="std.diff", ##dimnames(x$results)[["stat"]],
                    which.vars=dimnames(x$results)[["vars"]],thevarlabs=which.vars,
                    thexlab="Standardized Differences",
                    thecols=rainbow(length(which.strata)),
                    thesymbols=c(19,22,23,24,25)[1:length(which.strata)],
                    absolute = FALSE,
                    ordered = FALSE,
                    ...){

  # the helper .plot.xbal turns xb and the many of the arguments into a
  # an ordered variables by stratas table. The colnames are the labeled stratifications
  # the rownames are the labeled variables
  theresults <- .plot.xbal(x, 
                           which.strata,
                           thestratalabs,
                           which.stat,
                           which.vars,
                           thevarlabs,
                           absolute,
                           ordered)
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

.plot.xbal <- function(x, 
                       
                       which.strata,
                       thestratalabs,
                       which.stat,
                       which.vars,
                       thevarlabs,
                       absolute,
                       ordered) {

  if (!(which.stat %in% dimnames(x$results)[["stat"]])){
    stop(paste(which.stat,' not among results recorded in xbal object.'))}
  
  if (length(which.stat) > 1) {
    stop("Only one statistic allowed per-plot")
  }

  tmp <- !(which.vars %in% dimnames(x$results)[["vars"]])
  if (any(tmp)) {
    stop(paste("Unknown variable(s):", paste(which.vars[tmp], collapse = ",")))
  }

  theresults <- x$results[which.vars, which.stat, which.strata, drop=FALSE]
  # theresults are still an array, but there is guaranteed only one statistic
  # so we can cast to a data.frame to get a vars by strata table
  theresults <- as.data.frame(theresults)
  
  # we assume that the labels are in the same order as the which.xxx vars
  # and add potentially pretty labels to the rows and columns
  rownames(theresults) <- thevarlabs
  colnames(theresults) <- thestratalabs
  

  if (absolute) {
    theresults <- abs(theresults) 
  }

  if(ordered) {
    # the data are ordered using the  statistic in the first stratifying factor
    tmp <- order(theresults[,1])  
    theresults <- theresults[tmp,, drop = FALSE]
  } 

  return(theresults)
}

#' @export
balanceplot <- function(x, groups = NULL, ...) {
  nvars <- dim(x)[1]
  nstrat <- dim(x)[2]


  xrange <- range(x, na.rm = TRUE)
  xrange <- xrange + xrange * 0.25
  
  ypos <- 1:nvars
  
  if (!is.null(groups)) {
    groups <- as.factor(groups)

    # order X by the groups, and within groups order by the first column
    localorder <- order(groups, x[,1])
    x <- x[localorder, ]
    groups <- groups[localorder]

    grpnames <- levels(groups)
    pergroup <- table(groups)[grpnames]

    ypos <- ypos + unlist(mapply(rep, 0:(length(grpnames) - 1), pergroup))
  } else {
    # redorder X by the first variable, could be an option in the future
    x <- x[order(x[,1]), ]
  }

  mai <- par('mai')
  mai[2] <- max(c(strwidth(rownames(x), units = "inches") + mai[2], mai[2]))
  par(mai = mai)
  
  plot(xrange, 
       range(ypos) + c(0,1) ,
       axes = FALSE,
       pch = 19,
       col = "blue",
       ylab = "",
       xlab = "Balance",
       type = "n",
       ...)

  for(i in 1:nstrat) {
    points(x[,i], ypos, pch = i, cex = 0.5) # col =thecols[i],pch=thesymbols[i])
  }

  axis(1, at = pretty(seq(xrange[1], xrange[2], length = 5)))
  axis(2, labels = rownames(x), at = ypos, las = 2, tick = FALSE)

  if (!is.null(groups)) {
    text(mean(xrange), cumsum(pergroup) + (1:length(grpnames)), labels = grpnames)
  }

} 
