% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Biomarker.R
\name{nomogram.plot}
\alias{nomogram.plot}
\title{Nomogram plot by rms}
\usage{
nomogram.plot(fit = NULL, data = NULL, fun.list = NA, lp = F)
}
\arguments{
\item{fit}{rms model}

\item{data}{The data used to build the model}

\item{fun.list}{an optional function to transform the linear predictors, and to plot on another axis. If more than one transformation is plotted, put them in a list, e.g. list(function(x) x/2, function(x) 2*x). Any function values equal to NA will be ignored.}

\item{lp}{线性预测If fun.list is NA, lp will be TRUE. Set to FALSE to suppress creation of an axis for scoring X beta}
}
\description{
Nomogram plot by rms
}
\examples{
# Logistic model
data(LIRI)
res = loonR::build.logistic.model(LIRI[,3:5], LIRI$status, rms = T, scale = F)
f1 = list(Risk = loonR::logit2prob)
loonR::nomogram.plot(res$model, res$data, lp = T)
loonR::nomogram.plot(res$model, res$data, f1, lp =F)

# Survial model
res = build.psm.regression.model(LIRI[,3:5],LIRI$status, LIRI$time, scale = F)

surv <- Survival(res$model) # This would also work if f was from cph
surv_100 <- function(x) surv(100, lp = x)
surv_300 <- function(x) surv(300, lp = x)

med <- Quantile(res$model)
med_f <- function(x) med(lp=x)

f.list=list(
  `Median Survival Time`= med_f,
  `Probability of 100-day Survival`=surv_100,
  `Probability of 300-day Survival`=surv_300
)
nomogram.plot(res$model, res$data, fun.list = f.list, lp =F)
}
