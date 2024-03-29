% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCA.R
\name{ntbft}
\alias{ntbft}
\title{Net benifit}
\usage{
ntbft(
  response,
  pred,
  xstart = 0.01,
  xstop = 0.99,
  step = 0.01,
  type = "treated",
  model.name = "Model",
  exterdt = NULL
)
}
\arguments{
\item{response}{Should be 0 and 1}

\item{xstart}{the starting point of the threshold probability, the default value is 0.01.}

\item{xstop}{the end point of the threshold probability, the default value is 0.99}

\item{step}{a numerical value specifying the incremental step of the threshold probability, the default value is 0.01}

\item{type}{controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”}

\item{exterdt}{specify whether the prediction should be made in external dataset; this is useful for cross validation of the model. By default, the value is NULL, indicating the prediction is made in the same dataset as the training dataset}
}
\description{
Net benifit
}
\examples{
data(LIRI)

pred <- as.vector( unlist( LIRI[,c(3)] ) )
ntbft(LIRI$status, pred)

}
