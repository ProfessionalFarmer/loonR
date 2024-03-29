% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qPCR.R
\name{checkHousekeepingCombination}
\alias{checkHousekeepingCombination}
\title{Iterate all the combination of housekeeping genes}
\usage{
checkHousekeepingCombination(
  raw.ct.df,
  housekeeping.ct.df,
  group,
  alternative = "less",
  p.cutoff = 0.05,
  difference.cutoff = 0,
  sep = ", ",
  logTen = FALSE,
  expTwo = FALSE,
  scale = TRUE
)
}
\arguments{
\item{raw.ct.df}{row is sample}

\item{housekeeping.ct.df}{row is sample}

\item{group}{Second shuld be TRUE}

\item{alternative}{Default "less". One of c("less", "greater", "both")}

\item{p.cutoff}{Default 0.05}

\item{difference.cutoff}{Default 0. Should be greater than 0. The direction of difference depends on alternative.}

\item{sep}{Default ", "}

\item{logTen}{log10(2^(-ΔCT))}

\item{expTwo}{2^(-ΔCT)}

\item{scale}{}
}
\description{
Iterate all the combination of housekeeping genes
}
