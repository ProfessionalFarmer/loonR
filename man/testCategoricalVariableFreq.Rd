% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNV.R
\name{testCategoricalVariableFreq}
\alias{testCategoricalVariableFreq}
\title{Categorical variable test}
\usage{
testCategoricalVariableFreq(
  sample1,
  sample2,
  s1.name = "Sample1",
  s2.name = "Sample2",
  f = "chisq.test"
)
}
\arguments{
\item{f}{Default chisq.test, can be fisher.test}
}
\description{
Categorical variable test
}
\examples{
g1 = sample(c("Amp","N","Del"), 10, replace = T)
g2 = c( sample(c("Amp","N","Del"), 15, replace = T), sample(c("Del"), 5, replace = T) )
res = loonR::testCategoricalVariableFreq(g1, g2)
res$P
}
