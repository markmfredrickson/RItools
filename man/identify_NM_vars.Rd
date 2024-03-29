% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Design.R
\name{identify_NM_vars}
\alias{identify_NM_vars}
\title{Identify vars recording not-missing (NM) info}
\usage{
identify_NM_vars(vnames)
}
\arguments{
\item{vnames}{character, variable names}
}
\value{
character vector of names of NM vars, possibly of length 0
}
\description{
ID variables recording NM information, from
names and positions in the variable list.
Presumption is that the NM cols appear at the
end of the list of vars and are encased in
\sQuote{()}.  If something in the code changes
to make this assumption untrue, then this
helper is designed to err on the side of not
identifying other columns as NM cols.
}
\author{
Hansen
}
\keyword{internal}
