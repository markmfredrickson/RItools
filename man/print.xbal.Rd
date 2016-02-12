\name{print.xbal}
\alias{print}
%\alias{xBalance}
%\alias{xbal}
\alias{print.xbal}
\title{Printing xBalance Objects}
\description{
 A print method for xBalance objects.
}
\usage{
\method{print}{xbal}(x,
which.strata=dimnames(x$results)[["strata"]],
which.stats=dimnames(x$results)[["stat"]],
which.vars=dimnames(x$results)[["vars"]],
print.overall=TRUE,digits = NULL,printme=TRUE,
show.signif.stars=getOption("show.signif.stars"),
show.pvals=!show.signif.stars,
horizontal=TRUE,
report = NULL,...)
}

\arguments{
  \item{x}{An object of class "xbal" --- the result of a call to
    xBalance()}
  \item{which.strata}{The stratification candidates to include in the
    printout. Default is all.}
  \item{which.stats}{The test statistics to include. Default is all those
    requested from the call to xBalance.}
  \item{which.vars}{The variables for which test information should be
    displayed. Default is all.}
  \item{print.overall}{Should the omnibus test be reported? Default is
    TRUE.}
  \item{digits}{To how many digits should the results be displayed?
    Default is \code{max(2,getOptions("digits")-4)}.}
  \item{printme}{Print the table to the console? Default is TRUE}
  \item{show.signif.stars}{Use stars to indicate z-statistics larger than
    conventional thresholds. Default is TRUE.}
  \item{show.pvals}{Instead of stars, use p-values to summarize the
    information in the z-statistics. Default is FALSE}
  \item{horizontal}{Display the results for different candidate
    stratifications side-by-side (Default, TRUE), or as a list for each
    stratification (FALSE).}
  \item{report}{Override the reporting arguments originally passed to
  \code{xBalance}. See the help file for \code{\link{xBalance}} for options.}
  \item{...}{Other arguments. Not currently used.}
}

%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\value{   
  \item{vartable}{The formatted table of variable-by-variable statistics
  for each stratification.}
  \item{overalltable}{If the overall Chi-squared statistic is requested,
  a formatted version of that table is returned.}
  }
\references{
\code{\link{xBalance}}  
}
\examples{
data(nuclearplants)

xb0<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n,data=nuclearplants)

print(xb0)

xb1<-xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
         data = nuclearplants,
         report = c("all"))

str(xb1)

print(xb1)

print(xb1, show.pvals = TRUE)

print(xb1, horizontal = FALSE)

## The following doesn't work yet. 
\dontrun{print(xb1,which.vars=c("date","t1"),which.stats=c("adj.means","z.scores","p.values"))}

## The following example prints the adjusted means

print(xb1,
      which.vars = c("date", "t1"),
      which.stats = c("Control", "Treatment", "z", "p"))

## Now, not asking for the omnibus test
xb2 <- xBalance(pr~ date + t1 + t2 + cap + ne + ct + bw + cum.n + strata(pt),
         data = nuclearplants,
         report = c("all"))

print(xb2, which.strata = "pt")

}
\keyword{print}
