% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Common.R
\name{show_hcluster}
\alias{show_hcluster}
\title{Perform hclustering analysis}
\usage{
show_hcluster(
  df,
  group = NULL,
  dist.method = "euclidean",
  hclust.method = "ward.D2",
  color.pla = "npg",
  main = "",
  cutree = 0
)
}
\arguments{
\item{df}{row is gene, col is sample}

\item{dist.method}{Default euclidean. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.}

\item{hclust.method}{Default ward.D2. The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).}

\item{color.pla}{Default npg}

\item{main}{Title}

\item{cutree}{Number of clusters}
}
\description{
Perform hclustering analysis
}
\examples{
data(LIRI)
loonR::show_hcluster(t(LIRI[,3:5]), LIRI$status)
}
