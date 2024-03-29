% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNA.R
\name{limma_differential}
\alias{limma_differential}
\title{Differential analysis by LIMMA}
\usage{
limma_differential(
  df,
  group,
  rawcount = FALSE,
  voom = FALSE,
  pre.filter = 0,
  prop.expressed.sample = 0.5,
  cal.AUC = TRUE
)
}
\arguments{
\item{df}{raw count or log2(TPM+1), default log2(TPM+1). Must be data.frame, not matrix.}

\item{group}{factor, first control then experiment}

\item{rawcount}{true or false}

\item{voom}{true or false. If library size changed too much.}

\item{pre.filter}{Cutoff for mean log2(TPM)}

\item{prop.expressed.sample}{Default 0.5. Proportion of samples have a count greater than pre.filter}

\item{cal.AUC}{If to calculate AUC}
}
\description{
Differential analysis by LIMMA
}
\examples{
loonR::limma_differential(tpm.table, group)
}
