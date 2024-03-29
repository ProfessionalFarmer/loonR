% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSE.R
\name{simplifyEnrichment}
\alias{simplifyEnrichment}
\title{GO annotation with word cloud}
\usage{
simplifyEnrichment(go_ids, measure = "Rel", cluster_method = "binary_cut")
}
\arguments{
\item{go_ids}{A vector of GO IDs}

\item{measure}{Semantic measure for the GO similarity, see https://rdrr.io/pkg/GOSemSim/man/termSim.html}

\item{cluster_method}{https://jokergoo.github.io/simplifyEnrichment/reference/cluster_terms.html}
}
\description{
GO annotation with word cloud
}
