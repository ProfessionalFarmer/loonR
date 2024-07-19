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



#' Load miRDeep2 expression from a directory
#'
#' @param dir.path one folder is one sample
#' @param add.names2df Default TRUE. If add miRNA and precusor names in the first two columns
#' @param message If print message
#'
#' @return
#' @export
#'
#' @examples
#'
#' # 文件夹下有样本文件夹，样本表达为类似miRNAs_expressed_all_samples_06_09_2024.csv
#' dir.path="~/work/ProstaticCancer/analysis/321samples"
#' read_mirdeep2_result(dir.path)
#'
read_mirdeep2_resultDir <- function(dir.path, add.names2df=TRUE, message = F){

  samples = list.dirs(dir.path,recursive = F, full.names = F)
  samples = setdiff(samples, "multiqc")
  library(tidyverse)


  samples.df.list = lapply(samples, function(s){

    if(message){cat("Reading ", s,"\n")}


    exp.f = list.files(file.path(dir.path,s), "miRNAs_expressed_all_samples", full.names = T)
    exp.f = read.table(exp.f, header = T, quote = "", sep = "\t", comment.char = "%")
    # [1] "X.miRNA"    "read_count" "precursor"  "total" "seq"        "seq.norm."

    colnames(exp.f) = c("miRNA",
                        paste0(s, "_", "counts" ),
                        "precursor",
                        paste0(s, "_", "total" ),
                        paste0(s, "_", "seq" ),
                        paste0(s, "_", "CPM" )
                        )
    exp.f = exp.f[,c(1,3,2,6)]

    exp.f = exp.f %>% arrange(desc(miRNA), desc(precursor))

    exp.f
  })

  names(samples.df.list) = samples

  ########### check names
  samples.df.list = lapply(samples, function(s){

    if(message){cat("Checking ", s)  }

    s = samples.df.list[[s]]
    # 所有样本的前两列都要保持一致
    error = sum(s$miRNA!= samples.df.list[[1]]$miRNA)
    error = error + sum(s$precursor!= samples.df.list[[1]]$precursor)
    if(error!=0){
      stop(s, " Stop: miRNA and precusor not match")
    }
    if(message){ cat("-----Pass ", "\n") }
    s
  })
  names(samples.df.list) = samples


  print("Merging........")
  mirdeep2.res = reduce(samples.df.list, cbind)
  mirdeep2.res = mirdeep2.res[,c(1,2,
                                 seq(3, ncol(mirdeep2.res),4),
                                 seq(4, ncol(mirdeep2.res),4)
                                 )]

  count <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_counts")]
  cpm <- mirdeep2.res[, stringr::str_detect(colnames(mirdeep2.res), "_CPM")]

  colnames(count) <- stringr::str_remove_all(colnames(count), "_counts")
  colnames(cpm) <- stringr::str_remove_all(colnames(cpm), "_CPM")
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

