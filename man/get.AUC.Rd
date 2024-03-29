% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biomarker.R
\name{get.AUC}
\alias{get.AUC}
\title{Obtain AUC value}
\usage{
get.AUC(pred, label, raw = TRUE, direction = "auto")
}
\arguments{
\item{pred}{Predicted score or probability}

\item{label}{label/class/group}

\item{raw}{Default TRUE, return raw auc result with description}

\item{direction}{<: (controls < t <= cases) ######  >: (controls > t >= cases)}
}
\description{
Obtain AUC value
}
\examples{
data("LIRI")
loonR::get.AUC(LIRI$ANLN, LIRI$status)

}
