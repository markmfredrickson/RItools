% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.xbal.R
\name{balanceplot}
\alias{balanceplot}
\title{Create a plot of the balance on variables across different stratifications.}
\usage{
balanceplot(
  x,
  ordered = FALSE,
  segments = TRUE,
  colors = "black",
  shapes = c(15, 16, 17, 18, 0, 1, 10, 12, 13, 14),
  segments.args = list(col = "grey"),
  points.args = list(cex = 1),
  xlab = "Balance",
  xrange = NULL,
  groups = NULL,
  tiptext = NULL,
  include.legend = TRUE,
  legend.title = NULL,
  plotfun = .balanceplot,
  ...
)
}
\arguments{
\item{x}{A matrix of variables (rows) by strata (columns).}

\item{ordered}{Should the variables be ordered from
most to least imbalance on the first statistic?}

\item{segments}{Should lines be drawn between points for each
variable?}

\item{colors}{Either a vector or a matrix of shape indicators
suitable to use as a \code{col} argument to the
\code{\link{points}} function. If the argument is a vector, the
length should be the same as the number of columns in \code{x}. If
the argument is a matrix, it should have the same dims as \code{x}.}

\item{shapes}{Either a vector or a matrix of shape indicators
suitable to use as a \code{pch} argument to the
\code{\link{points}} function. If the argument is a vector, the
length should be the same as the number of columns in \code{x}. If
the argument is a matrix, it should have the same dims as
\code{x}.
<!-- The suggested vector has been selected to work with
RSVGTipsDevice tool tips.-->}

\item{segments.args}{A list of arguments to pass to the
\code{\link{segments}} function.}

\item{points.args}{A list of arguments to pass to the \code{\link{points}} function.}

\item{xlab}{The label of the x-axis of the plot.}

\item{xrange}{The range of x-axis. By default, it is 1.25 times the range of \code{x}.}

\item{groups}{A factor that indicates the group of each row in
\code{x}. Groups are printed under a common header.}

\item{tiptext}{ignored (legacy argument retained for internal reasons)
<!-- If you are using the \code{RSVGTipsDevice} library for
rendering, you can include an array of the dimensions of x
with another dimension of length 2. For example, if there are 4
observations and 2 strata, the array should be 4 by 2 by 2. The
\code{tiptext[i, j, 1]} entry will be the first line of the tool
tip for the data in \code{x[i, j]}. Likewise for the second row of
the tool tip. -->}

\item{include.legend}{Should a legend be included?}

\item{legend.title}{An optional title to attach to the legend.}

\item{plotfun}{Function to do the plotting; defaults to [RItools:::.balanceplot]}

\item{...}{Additional arguments to pass to \code{\link{plot.default}}.}
}
\value{
Returns NULL, displays plot
}
\description{
This plotting function summarizes variable by stratification matrices. For
each variable (a row in the \code{x} argument), the values are under each
stratification (the columns of \code{x}) plotted on the same line.
}
\details{
It is conventional to standardize the differences to common scale
(e.g.  z-scores), but this is not required. When \code{ordered} is
set to true, plotting will automatically order the data from
largest imbalance to smallest based on the first column of
\code{x}.

You can fine tune the colors and shapes with the like named
arguments. Any other arguments to the \code{\link{points}} function
can be passed in a list as \code{points.args}. Likewise, you can
fine tune the segments between points with \code{segments.args}.
}
\examples{
set.seed(20121204)

# generate some balance data
nvars <- 10
varnames <- paste("V", letters[1:nvars])

balance_data <- matrix(c(rnorm(n = nvars, mean = 1, sd = 0.5), 
                         rnorm(n = nvars, mean = 0, sd = 0.5)),
                       ncol = 2)

colnames(balance_data) <- c("Before Adjustment", "After Matching")

rownames(balance_data) <- varnames

balanceplot(balance_data,
                      colors = c("red", "green"),
                      xlab = "Balance Before/After Matching")

# base R graphics are allowed

abline(v = colMeans(balance_data), lty = 3, col = "grey")

}
\seealso{
\code{\link{plot.xbal}}, \code{\link{xBalance}},
\code{\link{segments}}, \code{\link{points}}
}
