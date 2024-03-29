% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biomarker.R
\name{build.elastic.net.regression}
\alias{build.elastic.net.regression}
\title{Build elastic net regression}
\usage{
build.elastic.net.regression(
  df,
  group,
  seed = 666,
  scale = TRUE,
  nfolds = 10,
  summaryFunction = "twoClassSummary"
)
}
\arguments{
\item{seed}{Default 666}

\item{scale}{Default TRUE}

\item{nfolds}{Default 10. number of folds}

\item{summaryFunction}{Default twoClassSummary for computes sensitivity, specificity and the area under the ROC curve. Please refer https://www.rdocumentation.org/packages/caret/versions/6.0-90/topics/defaultSummary}
}
\description{
Build elastic net regression
}
\examples{
data("LIRI")
reg.res <- loonR::build.elastic.net.regression(LIRI[,3:5],LIRI$status)

}
