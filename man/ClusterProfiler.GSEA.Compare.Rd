% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSE.R
\name{ClusterProfiler.GSEA.Compare}
\alias{ClusterProfiler.GSEA.Compare}
\title{Perform GSE analysis across multiple group by cluster profiler}
\usage{
ClusterProfiler.GSEA.Compare(
  gene,
  minGSSize = 10,
  qvalue = 0.05,
  exp.gene.type = "ENSEMBL",
  GO = FALSE,
  KEGG = FALSE,
  MSigDb = TRUE,
  Hallmark = TRUE
)
}
\arguments{
\item{gene}{A list including multiple group. Vector names are gene.key and vector value are log fold chage}

\item{minGSSize}{10}

\item{qvalue}{0.05}

\item{exp.gene.type}{key types of input phenotype. Default ENSEMBL, can by keytypes(org.Hs.eg.db)}

\item{GO}{TRUE}

\item{KEGG}{TRUE}

\item{MSigDb}{TRUE}

\item{Hallmark}{TRUE}
}
\description{
Perform GSE analysis across multiple group by cluster profiler
}
