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



#' Get TCGA miRNA expression data
#'
#' @param tcga.project See ?RTCGA.miRNASeq::miRNASeq. For example: COAD
#' @param rawCount If to download raw count data
#' @param CPM If to download CPM data
#' @param log2 If to perform log2 transformation
#'
#' @return
#' @export
#'
#' @examples
get.TCGA.miRNAExpression <- function(tcga.project=NULL, rawCount=FALSE, CPM=FALSE, log2=FALSE){

  if(is.null(tcga.project)){
    warning(?RTCGA.miRNASeq::miRNASeq)
    stop("Please specify a TCGA project. See ?RTCGA.miRNASeq::miRNASeq")
  }
  if(!rawCount & !CPM){
    stop("Please specify download rawCount or CPM. Should only set one as TRUE")
  }else if(rawCount & CPM){
    stop("Should only set one as TRUE. Raw Count or CPM?")
  }



  if (!require(RTCGA))BiocManager::install("RTCGA")
  if (!require(RTCGA.miRNASeq))BiocManager::install("RTCGA.miRNASeq")
  if (!require(RTCGA.clinical))BiocManager::install("RTCGA.clinical")

  # Expression
  ##-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###
  expression.df = get( paste0(tcga.project,".miRNASeq") )

  if(rawCount){
    col.ind <- grep("read_count", expression.df$miRNA_ID)
  }else if(CPM){
    col.ind <- grep("reads_per_million_miRNA_mapped", expression.df$miRNA_ID)
  }

  expression.df = expression.df[col.ind,-c(1,2)]

  patient_ids <- rownames(expression.df)

  # convert string to numeric
  expression.df <- apply(as.matrix.noquote(expression.df),2,as.numeric)

  patient_ids -> rownames(expression.df)

  if(log2){
    expression.df = t(log2(expression.df+1))
  }else{
    expression.df = t( expression.df )

  }

  # Clinical
  ##-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###
  # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  group_list <- ifelse(substr(patient_ids,14,15)=='01','Tumor','Normal')

  library(dplyr)
  if(any( ! substr(patient_ids,14,15) %in% c(01,11)  )){
    stop("It seems some samples are not primary tumor or normal samples")
  }


  meta <- get( paste0(tcga.project,".clinical") )
  meta$patient.bcr_patient_barcode <- stringr::str_to_upper( meta$patient.bcr_patient_barcode )

  meta <- meta[match(substr(patient_ids,1,12),
                     meta$patient.bcr_patient_barcode),]

  row.names(meta) <- patient_ids
  meta$Label = group_list

  res = list(meta=meta,
             expr=expression.df,
             label=group_list)

  res
}





