plot.xbal<-function(x,adjustxaxis=.25,segments=TRUE,legend=TRUE,
                    mar=c(3,3,2,0)+0.1,mgp=c(1.5,.5,0),tck=-.01,
                    which.strata=dimnames(x$results)[["strata"]],thestratalabs=which.strata,
                    which.stats="std.diff", ##dimnames(x$results)[["stat"]],
                    which.vars=dimnames(x$results)[["vars"]],thevarlabs=which.vars,
                    thexlab="Standardized Differences",
                    thecols=rainbow(length(which.strata)),
                    thesymbols=c(19,22,23,24,25)[1:length(which.strata)],...){
  if (!("std.diff"%in%dimnames(x$results)[["stat"]]))
    stop('"std.diff" not among results recorded in xbal object.') 
  theresults<-x$results[which.vars,which.stats,which.strata,drop=FALSE]
  ypos<-seq(length=length(which.vars))
  xrange<-range(theresults[,which.stats,],na.rm=TRUE)
  xrange<-xrange+xrange*adjustxaxis

  ##Setup the margin, adjust for the lengths of the labels.
  par(mar=mar,tck=tck,mgp=mgp) ##set default margins
  mymai<-par('mai')
  mymai[2]<-max(c(strwidth(thevarlabs,units="inches"),mymai[2]))
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
    points(theresults[which.vars,which.stats,i],ypos,col=thecols[i],pch=thesymbols[i])
    }
  if(segments&length(which.strata)>1){ ##segments are mainly useful for drawing the eye to changes in balance along a single variable across more than 1 stratification
    for(j in ypos){
      segments(min(theresults[j,which.stats,which.strata]),j,
               max(theresults[j,which.stats,which.strata]),j,col=gray(.7),lwd=.5) 
    }
  }
  axis(1,at=pretty(seq(xrange[1],xrange[2],length=5)))
  axis(2,labels=thevarlabs,at=ypos,las=2,tick=FALSE)
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
