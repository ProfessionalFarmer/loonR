% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNA.R
\name{check_diff_gene}
\alias{check_diff_gene}
\title{Check gene expression}
\usage{
check_diff_gene(
  gene,
  genes_expr,
  group_list,
  color = "aaas",
  stat = NULL,
  comparisions = NULL,
  outlier.shape = 19,
  add = "jitter"
)
}
\arguments{
\item{stat}{stat method}

\item{comparisions}{list(c("A","B"))}

\item{outlier.shape}{Default 19. To hide outlier, specify outlier.shape = NA. When jitter is added, then outliers will be automatically hidden.}

\item{add}{Default is jitter}
}
\description{
Check gene expression
}
