% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biomarker.R
\name{univariate_or}
\alias{univariate_or}
\title{Variate logistic analysis
Row: sample, Column: gene expression
score E.g.: Gene or miRNA expression, or risk score}
\usage{
univariate_or(d.frame, label)
}
\arguments{
\item{d.frame}{Data.frame --- Row: sample, Column: gene expression}

\item{label}{True Sample label}
}
\value{
c(OR, 2.5% CI, 97.5% CI)
}
\description{
Variate logistic analysis
Row: sample, Column: gene expression
score E.g.: Gene or miRNA expression, or risk score
}
\examples{
univariate_or(risk.df, label)
}