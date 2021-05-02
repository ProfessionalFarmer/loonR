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
  dge <- estimateCommonDisp(dge, verbose=TRUE)
  cpm <- cpm(dge, log = log, prior.count = prior.count)
  cpm
}



#' Load miRDeep2 expression matrix
#'
#' @param file.path mirdeep2 expression. Should be merged before
#'
#' @return
#' @export
#'
#' @examples
read_mirdeep2_result <- function(file.path){

mirdeep2.res <- read.table(file.path, sep = "\t", header = T, check.names = F)
count <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_read_count")]
cpm <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_seq\\(norm\\)")]

colnames(count) <- stringr::str_remove_all(colnames(count), "_read_count")
colnames(cpm) <- stringr::str_remove_all(colnames(cpm), "_seq\\(norm\\)")
miRNA <- mirdeep2.res$miRNA
precursor <- mirdeep2.res$precursor

res <- list(Count = count, CPM = cpm, miRNA = miRNA, precursor = precursor)
res
}
