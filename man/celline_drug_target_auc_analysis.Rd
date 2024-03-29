% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Drug.R
\name{celline_drug_target_auc_analysis}
\alias{celline_drug_target_auc_analysis}
\title{Analyze if gene is correlated with drug sensitivity}
\usage{
celline_drug_target_auc_analysis(
  drug.sen.df = NULL,
  gene.df = NULL,
  compare.method = "wilcox.test"
)
}
\arguments{
\item{drug.sen.df}{Row is drug, Col is sample. Pls perform log transformation}

\item{gene.df}{Row is gene Col is sample}

\item{compare.method}{wilcox.test t.test}
}
\description{
Analyze if gene is correlated with drug sensitivity
}
