% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biomarker.R
\name{lasso_best_lamda}
\alias{lasso_best_lamda}
\title{Get best lamda by perform multiple round lasso-sv
https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
https://stackoverflow.com/questions/62170965/looping-cv-glmnet-and-get-the-best-coefficients}
\usage{
lasso_best_lamda(
  d.matrix,
  group,
  family = "binomial",
  type.measure = "auc",
  nfolds = 5,
  nreps = 100,
  scale = TRUE,
  seed = 66
)
}
\arguments{
\item{d.matrix}{Row is sample}

\item{family}{Default binomial}

\item{type.measure}{class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.}

\item{nfolds}{Default 5}

\item{nreps}{Default 1000}

\item{scale}{Default TRUE}
}
\description{
Get best lamda by perform multiple round lasso-sv
https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
https://stackoverflow.com/questions/62170965/looping-cv-glmnet-and-get-the-best-coefficients
}
