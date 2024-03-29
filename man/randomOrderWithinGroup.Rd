% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Common.R
\name{randomOrderWithinGroup}
\alias{randomOrderWithinGroup}
\title{Randomly place the sample within group}
\usage{
randomOrderWithinGroup(group = NULL, ids = NULL, seed = 666)
}
\arguments{
\item{group}{Group label}

\item{ids}{Sample names}

\item{seed}{Default 666}
}
\description{
Randomly place the sample within group
}
\examples{
group = rep(c("G1","G2"),c(10,10))
ids = 1:20
loonR::randomOrderWithinGroup(group, ids)
}
