% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Subtyping.R
\name{subtypeAssociationAnalysis}
\alias{subtypeAssociationAnalysis}
\title{Subtype association: Calculate Jaccard similarity coefficient and perform hypergeometrix test}
\usage{
subtypeAssociationAnalysis(
  df,
  concateStudy = F,
  adjusted.hypergeometrix.p = F,
  cut.edge.byPval = 0.05,
  print.message = T
)
}
\arguments{
\item{df}{row is sample, column is study subtype. Plsease note rownames is sample ID}

\item{concateStudy}{If to add study name (from colnames) into result}

\item{adjusted.hypergeometrix.p}{Default FALSE, if to adjust hypergeometrix P value by BH method.}

\item{cut.edge.byPval}{Default 0.05}

\item{print.message}{If print message}
}
\description{
Subtype association: Calculate Jaccard similarity coefficient and perform hypergeometrix test
}
\examples{
data(Subtype.matrix)
res <- loonR::subtypeAssociationAnalysis(Subtype.matrix)
res$hyperGeoTest.Analysis
}
