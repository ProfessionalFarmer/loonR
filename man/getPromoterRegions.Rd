% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methylation.R
\name{getPromoterRegions}
\alias{getPromoterRegions}
\title{Get promoters by GenomicFeatures package. Region needed}
\usage{
getPromoterRegions(
  upstream = 2000,
  downstream = 500,
  ann = "Ensembl",
  ref.genome = "hg38",
  ens.release = 99,
  addChr = TRUE
)
}
\arguments{
\item{upstream}{Default 2000}

\item{downstream}{Default 500}

\item{ann}{Default Ensembl}

\item{ref.genome}{Default hg38 (for refGene)}

\item{ens.release}{Default 99 (for Ensembl)}
}
\description{
Get promoters by GenomicFeatures package. Region needed
}
