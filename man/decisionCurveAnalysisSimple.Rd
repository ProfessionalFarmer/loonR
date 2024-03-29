% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCA.R
\name{decisionCurveAnalysisSimple}
\alias{decisionCurveAnalysisSimple}
\title{Decision curve analysis using risk probability and label}
\usage{
decisionCurveAnalysisSimple(
  label,
  pred,
  xstart = 0.01,
  xstop = 0.99,
  step = 0.01,
  type = "treated",
  model.name = "Model",
  boots = 0
)
}
\arguments{
\item{label}{Should be 0 and 1}

\item{xstart}{the starting point of the threshold probability, the default value is 0.01.}

\item{xstop}{the end point of the threshold probability, the default value is 0.99}

\item{step}{a numerical value specifying the incremental step of the threshold probability, the default value is 0.01}

\item{type}{controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”}

\item{model.name}{Default Model}

\item{boots}{100 Bootstrape to avoid overfitting}
}
\description{
Decision curve analysis using risk probability and label
}
\examples{
data(LIRI)

pred <- as.vector( unlist( LIRI[,c(3)] ) )
decisionCurveAnalysisSimple(LIRI$status, pred)

}
