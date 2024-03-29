% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Common.R
\name{meltDataFrameByGroup}
\alias{meltDataFrameByGroup}
\title{Melt data frame by group}
\usage{
meltDataFrameByGroup(
  d.frame = NULL,
  group = NULL,
  na.rm = TRUE,
  variable_name = "Gene",
  group2 = NULL
)
}
\arguments{
\item{d.frame}{Row is sample, column is feature.}

\item{na.rm}{Default TRUE. Should NA values be removed from the data set?}

\item{variable_name}{Default "Gene"}

\item{group2}{The second group}
}
\description{
Melt data frame by group
}
\examples{
data(LIRI)
d.frame = LIRI[,3:6]
group = LIRI$status
head( loonR::meltDataFrameByGroup(d.frame, group) )
}
