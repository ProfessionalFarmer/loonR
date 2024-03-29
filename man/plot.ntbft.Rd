% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCA.R
\name{plot.ntbft}
\alias{plot.ntbft}
\title{Plotting the Net benefit function. Coule be multiple curves}
\usage{
\method{plot}{ntbft}(
  nb,
  nolines = 2:dim(nb)[2],
  nobands = NULL,
  ymin = -0.1,
  ymax = max(nb[, c(nolines, nobands)], na.rm = T),
  legpos = c(0.9, 0.8),
  palette = "aaas"
)
}
\arguments{
\item{nb}{Object from}

\item{nolines}{the number of the columns of nb which should be plotted using lines. The default is to plot all columns (except the first one containing the threshold).}

\item{nobands}{the number of the columns of nb which should be plotted using bands (useful to plot confidence intervals). The default is to plot no bands}

\item{legpos}{a vector of two coordinates indicating where the legend should be in the graph}

\item{palette}{}
}
\description{
Plotting the Net benefit function. Coule be multiple curves
}
\examples{
p = loonR::decisionCurveAnalysisSimple(label, risk)
p$Plot
v.netbenefit = p$NetBenefit

# similar to getr.netbenefit

nb <- data.frame(threshold = v.netbenefit$threshold,
                 All = v.netbenefit$All,
                 None = v.netbenefit$None,
                 `Transcriptomic panel` = v.netbenefit$Model,
                 `Risk-stratification model` = r.netbenefit$Model,
                 stringsAsFactors = F,
                 check.names = F)

p = loonR:::plot.ntbft(nb, 2:5)
p

}
