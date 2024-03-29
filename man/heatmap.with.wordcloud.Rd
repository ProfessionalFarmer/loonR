% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap.R
\name{heatmap.with.wordcloud}
\alias{heatmap.with.wordcloud}
\title{Plot heatmap while perform word cloud analysis}
\usage{
heatmap.with.wordcloud(df)
}
\arguments{
\item{df}{row is sample, column is gene. Will perform gene correlation analysis instead of sample-wise}
}
\description{
Plot heatmap while perform word cloud analysis
}
\examples{
# https://jokergoo.github.io/simplifyEnrichment/articles/word_cloud_anno.html

df = data.frame(`a b c dd` = c(1,2,1,1,1), `a d c bb a c` =  c(3,2,3,1,1), `bb aa cc dd` = c(4,2,5,7,6))
M = cor(df)
# show correlation
corrplot(M, method = 'number', order = 'hclust')




}
