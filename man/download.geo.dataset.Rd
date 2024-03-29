% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEO.R
\name{download.geo.dataset}
\alias{download.geo.dataset}
\title{Download dataset from GEO by accession ID and platform}
\usage{
download.geo.dataset(
  geo.accession.id,
  platform = NULL,
  destdir = "~/GSE",
  platform.available = TRUE,
  normalizeBetweenArrays = FALSE,
  id.mapping = TRUE
)
}
\arguments{
\item{geo.accession.id}{GEO Accession ID}

\item{platform}{Platform}

\item{destdir}{Default tempdir()}

\item{platform.available}{If GLP platform not available}

\item{normalizeBetweenArrays}{Default FALSE. If perform normlize between arrays. Check boxplot first}

\item{id.mapping}{Default TRUE in order to select the maximum expression value if a gene has muiltiple probes}
}
\value{
list(expression, phenotype, probe.annotation)
}
\description{
Download dataset from GEO by accession ID and platform
}
