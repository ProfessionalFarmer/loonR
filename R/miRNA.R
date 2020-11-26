#' Title
#'
#' @param count.df
#' @param group
#' @param method "TMM","TMMwsp","RLE","upperquartile","none"
#' @param log
#'
#' @return
#' @export
#'
#' @examples
normalize_miRNACount_toCPM <- function(count.df, group, method = "TMM", log = TRUE, prior.count = 1){

  library(edgeR)

  dge <- DGEList(count.df, group = group)
  dge <- calcNormFactors(dge, method = method)
  #dge <- estimateCommonDisp(dge, verbose=TRUE)
  cpm <- cpm(dge, log = log, prior.count = prior.count)
  cpm
}


