#' Normalize miRNA counts to CPM
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
#' @param add.names2df Default TRUE. If add miRNA and precusor names in the first two columns
#'
#' @return
#' @export
#'
#' @examples
read_mirdeep2_result <- function(file.path, add.names2df=TRUE){

mirdeep2.res <- read.table(file.path, sep = "\t", header = T, check.names = F)
count <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_read_count")]
cpm <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_seq\\(norm\\)")]

colnames(count) <- stringr::str_remove_all(colnames(count), "_read_count")
colnames(cpm) <- stringr::str_remove_all(colnames(cpm), "_seq\\(norm\\)")
miRNA <- mirdeep2.res$miRNA
precursor <- mirdeep2.res$precursor

if(add.names2df){
  count <- tibble::add_column(count, miRNA = miRNA, precursor = precursor, .before = 1)
  cpm   <- tibble::add_column(cpm, miRNA = miRNA, precursor = precursor, .before = 1)
}

res <- list(Count = count, CPM = cpm, miRNA = miRNA, precursor = precursor)
res
}


#' miRBase convert
#'
#' @param miR.IDs A character vector representing the source miRNA names needed to be convert.
#' @param miRBase.IDs miRBase accessions. "MIMAT0000095"
#' @param version A character value representing the target miRBase version corresponding the Accessions. Default "v22"
#'
#' @return
#' @export
#'
#' @examples
miRBaseConverter <- function(miR.IDs  =NULL, miRBase.IDs = NULL, version = NULL){

  if(!require("miRBaseConverter")){
    devtools::install_github("taoshengxu/miRBaseConverter")
    library("miRBaseConverter")
  }

  if(is.null(miR.IDs) & is.null(miRBase.IDs) ){
    stop("Pls specify miR or miRBase IDs")
  }

  if(is.null(miR.IDs)){
    if(is.null(version)){version="v22"}

    res = miRNA_AccessionToName(miRBase.IDs, targetVersion = version )
    colnames(res) = c("miRBase", "miR")
  }

  if(is.null(miRBase.IDs)){
    if(is.null(version)){version=checkMiRNAVersion(miR.IDs, verbose = TRUE)}
    res = miRNA_NameToAccession(miR.IDs, version = version)
    colnames(res) = c("miR", "miRBase")
    res = res[,c(2,1)]
  }
  res

}

