% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Common.R
\name{scaleDF}
\alias{scaleDF}
\title{Scale a data.frame by column or row}
\usage{
scaleDF(
  df,
  byRow = FALSE,
  byColumn = FALSE,
  center = TRUE,
  scale = TRUE,
  maxUnit = 4,
  returnRaw = F
)
}
\arguments{
\item{byRow}{Default FALSE}

\item{byColumn}{Default FALSE, by column}

\item{center}{Default TRUE. Mean = 0}

\item{scale}{Default TRUE. 0-1 scale}

\item{maxUnit}{Default 4}

\item{returnRaw}{if TRUE, return raw dataframe}
}
\description{
Scale a data.frame by column or row
}
\examples{
Default by column
}
