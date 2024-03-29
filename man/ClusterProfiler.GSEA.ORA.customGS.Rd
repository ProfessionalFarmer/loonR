% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSE.R
\name{ClusterProfiler.GSEA.ORA.customGS}
\alias{ClusterProfiler.GSEA.ORA.customGS}
\title{Perform analysis using user custom gene sets}
\usage{
ClusterProfiler.GSEA.ORA.customGS(
  g,
  CustomGS = NULL,
  gse = FALSE,
  ova = FALSE,
  minGSSize = 10,
  qvalue = 0.05
)
}
\arguments{
\item{g}{GSE: values with genes symbols, OVA: symbols}

\item{CustomGS}{ClusterProfiler::read.gmt() term to symbols}

\item{gse}{Default FALSE, gene set enrichment analysis}

\item{ova}{Default FALSE, overrepresentation test}

\item{minGSSize}{10}

\item{qvalue}{0.05}
}
\description{
Perform analysis using user custom gene sets
}
