% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCA.R
\name{ntbft.boot}
\alias{ntbft.boot}
\title{Net benefit Bootstrap method to correct overfitting}
\usage{
ntbft.boot(
  outcome,
  pred,
  boots = 500,
  xstart = 0.01,
  xstop = 0.99,
  step = 0.01,
  type = "treated",
  model.name = "Model"
)
}
\arguments{
\item{boots}{500}

\item{xstart}{the starting point of the threshold probability, the default value is 0.01.}

\item{xstop}{the end point of the threshold probability, the default value is 0.99}

\item{step}{a numerical value specifying the incremental step of the threshold probability, the default value is 0.01}

\item{type}{controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”}

\item{model.name}{}

\item{response}{Should be 0 and 1}
}
\description{
Net benefit Bootstrap method to correct overfitting
}
\examples{
data(LIRI)

pred <- as.vector( unlist( LIRI[,c(3)] ) )
ntbft.boot(LIRI$status, pred)

}
