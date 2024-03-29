% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCA.R
\name{rms_validate}
\alias{rms_validate}
\title{Resampling Validation of a Fitted Model's Indexes of Fit}
\usage{
rms_validate(rms.model, B = 100)
}
\arguments{
\item{B}{Bootstrap times}
}
\description{
Resampling Validation of a Fitted Model's Indexes of Fit
}
\examples{
data(LIRI)
m=loonR::build.logistic.model(LIRI[,c(3,4)],LIRI$status, rms = T)
m=m$model
rms::validate(m)
}
