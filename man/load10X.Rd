% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRNA.R
\name{load10X}
\alias{load10X}
\title{Read data: raw to FindVariableFeatures}
\usage{
load10X(
  paths = NA,
  sample_names = NA,
  min.cells = 3,
  min.features = 200,
  gene.column = 1,
  remove_doublet = T,
  max_nCount_RNA = 10000,
  max_nFeature_RNA = 8000,
  max_percent.mt = 20,
  integrate_CCA = FALSE,
  top_variable_features = 2000,
  remove_cc = T
)
}
\arguments{
\item{min.cells}{3}

\item{min.features}{200}

\item{gene.column}{1}

\item{max_nFeature_RNA}{8000}

\item{max_percent.mt}{20}

\item{integrate_CCA}{Default FALSE}

\item{top_variable_features}{2000}
}
\description{
Read data: raw to FindVariableFeatures
}
