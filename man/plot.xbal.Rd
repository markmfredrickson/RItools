\name{plot.xbal}
\alias{plot}
%\alias{xBalance}
%\alias{xbal}
\alias{plot.xbal}
\title{Plotting xBalance Objects}
\description{
 A plot method for xBalance objects. WORK IN PROGRESS. PATCHES AND CODE
 SUGGESTIONS APPRECIATED.

 We aim this plot to be a diagnostic tool rather than a publication
 quality presentation tool.
}
\usage{
\method{plot}{xbal}(x,adjustxaxis=.25,segments=TRUE,legend=TRUE,
                    mar=c(3,3,2,0)+0.1,mgp=c(1.5,.5,0),tck=-.01,
                    which.strata=dimnames(x$results)[["strata"]],thestratalabs=which.strata,
                    which.stats="std.diff", ##dimnames(x$results)[["stat"]],
                    which.vars=dimnames(x$results)[["vars"]],thevarlabs=which.vars,
                    thexlab="Standardized Differences",
                    thecols=rainbow(length(which.strata)),
                    thesymbols=c(19,22,23,24,25)[1:length(which.strata)],...)
}

\arguments{
  \item{x}{An object of class "xbal" --- the result of a call to
    xBalance()}
  \item{adjustxaxis}{amount by which the x-axis should be expanded.}
  \item{segments}{Should thin horizontal lines be plotted connecting the
    statistics for the different stratifications.}
  \item{legend}{Should a legend be plotted?}
  \item{mar}{Preliminary margin setting}
  \item{mgp}{Preliminary setting for axis labels}
  \item{tck}{Length of the tick marks}
  \item{which.strata}{The stratification candidates to include in the
    printout. Default is all.}
  \item{thestratalabs}{The text labels for the strata.}
  \item{which.stats}{a character vector of length 1. The test statistics to include. Default is the
  standardized difference.}
  \item{which.vars}{The variables for which test information should be
    displayed. Default is all.}
  \item{thevarlabs}{The text labels for the variables.}
  \item{thexlab}{The label for the x-axis (should tell viewers about the
    statistic chosen).}
  \item{thecols}{A vector of colors either (a) one per strata or (b) one
    for all strata.}
  \item{thesymbols}{A vector of plotting symbols either (a) one per strata or (b) one
    for all strata.}
 \item{...}{other arguments to the \code{\link{plot.default}} function setting up the plotting region.}
}

%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\value{   
  The plot allows a quick visual comparison of the effect of different
  stratification designs on the comparability of different
  variables. This is not a replacement for the omnibus statistical test
  reported as part of \code{\link{print.xbal}}. This plot does allow the
  analyst an easy way to identify variables that might be the primary culprits
  of overall imbalances and/or a way to assess whether certain important
  covariates might be imbalanced even if the omnibus test reports that
  the stratification overall produces balance.
  }
\references{
\code{\link{xBalance}}  
}
\examples{
data(nuclearplants)

xb0<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,data=nuclearplants)

plot(xb0)

xb1<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,
         strata=data.frame(unstrat=factor(character(32)),
           pt=factor(nuclearplants$pt)),
         data=nuclearplants,
         report=c("adj.means","adj.mean.diffs","std.diffs", "z.scores", "chisquare.test","p.values"))

plot(xb1)

}
\keyword{hplot, print}
