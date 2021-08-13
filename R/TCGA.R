#' Get sample barcode and TCGA project map information
#'
#'
#' @param project A vector.
#'  TCGA-ACC, TCGA-COAD.
#'  See getGDCprojects()$project_id
#' @return A data.frame
#'  X1 proj samples
#'	1	1	TCGA-ACC	TCGA-OR-A5JZ
#'	2	1	TCGA-ACC	TCGA-OR-A5JK
#' @export
#'
#' @examples
get.sample.project.infomap <- function(project=NA){

  library(TCGAbiolinks)

  # If not specify project, load default
  if (is.na(project)){
    project <- grep("TCGA", sort(getGDCprojects()$project_id),
                    value = TRUE)
  }


  df <- plyr::adply(project, .margins = 1, .fun = function(proj) {
    samples <- TCGAbiolinks:::getSubmitterID(proj)
    return(data.frame(proj, samples))
  })

  return(df)

}


#' Convert FPKM data fram to TPM
#'
#' @param fpkm.exp.df
#'
#' @return
#' @export
#'
#' @examples
fpkm2tpm <- function(fpkm.exp.df){

  fpkm2tpm <- apply(fpkm.exp.df, 2, function(fpkm){
      tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
      tpm[which(is.na(tpm))] <- 0
      tpm
  })

  return(fpkm2tpm)
}



#' Download RAN expression by TCGAbiolinks
#'
#' @param project E.g. TCGA-LIHC
#' @param dir Raw data directory
#' @param remove.Raw If to remove raw data
#'
#' @return list(clinical = clinical, expression = expression.df, group = group)
#' Expression was log2 transformed.
#' Group is a data.frame including label, short name and sample barcode
#' Clinical is a list
#' @export
#'
#' @examples
tcgabiolinks.get.RNA.expression.log2tpm <- function(project, remove.Raw = FALSE, dir="~/rspace/GDCdata"){

  fpkm2tpm <- function(fpkm){
    tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
    tpm[which(is.na(tpm))] <- 0
    return(tpm)
  }


  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Transcriptome Profiling",
                                  workflow.type = "HTSeq - FPKM", # HTSeq - FPKM-UQ or HTSeq - Counts
                                  data.type = "Gene Expression Quantification",
                                  experimental.strategy = "RNA-Seq",
                                  legacy = FALSE
  )
  TCGAbiolinks::GDCdownload(query, directory = dir)
  project.data <- TCGAbiolinks::GDCprepare(query = query, directory = dir)



  # Expression
  expression.df <- SummarizedExperiment::assay(project.data)
  gene.information <- rowRanges(data)
  rm(project.data)
  # convert to tpm
  expression.df <- apply(expression.df, 2, fpkm2tpm)

  normal.sample <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                                       typesample = "NT")
  tumor.sample <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                                      typesample = "TP")

  group <- data.frame(
    Label = rep( c("Normal","Tumor"), c(length(normal.sample), length(tumor.sample))  ),
    Short.Name = substr(c(normal.sample,tumor.sample),1,12),
    Barcode = c(normal.sample,tumor.sample),
    stringsAsFactors = FALSE
  )

  # log2 transformation
  expression.df <- log2(expression.df[,c(normal.sample,tumor.sample)]+1)



  ## clinical
  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Clinical",
                                  data.type = "Clinical Supplement",
                                  data.format = "BCR Biotab")
  TCGAbiolinks::GDCdownload(query, directory = dir)
  clinical <- TCGAbiolinks::GDCprepare(query, directory = dir)



  if(remove.Raw){
    file.remove( paste(dir,"/",project,sep="",collapse = "")  )
  }

  result <- list(clinical = clinical,
                 expression = expression.df,
                 group = group,
                 Gene.info = gene.information,
                 query = query
  )


}




#' Download junction data by TCGAbiolinks. Pls note, V2 tag and Illumina HiSeq platform were selected
#'
#' @param project E.g. TCGA-LIHC
#' @param dir Raw data directory
#' @param remove.Raw If to remove raw data
#' @param tags V2
#'
#' @return list(clinical = clinical, expression = expression.df, group = group)
#' Group is a data.frame including label, short name and sample barcode
#' Clinical is a list
#' @export
#'
#' @examples
tcgabiolinks.get.junction.coverage <- function(project, remove.Raw = FALSE, dir="~/rspace/GDCdata", tag = 'v2', platform = "Illumina HiSeq"){
  library(TCGAbiolinks)
  query <- TCGAbiolinks::GDCquery(project = project,
                                  legacy = TRUE,
                                  data.category = "Gene expression",
                                  data.type = "Exon junction quantification",
                                  experimental.strategy = "RNA-Seq",
                                  platform = platform

  )

  query[[1]][[1]] <- query[[1]][[1]][grep(tag,query[[1]][[1]]$tags), ]



  TCGAbiolinks::GDCdownload(query, directory = dir)
  project.data <- TCGAbiolinks::GDCprepare(query = query, directory = dir)



  # Expression
  expression.df <- SummarizedExperiment::assay(project.data)
  gene.information <- rowRanges(data)
  rm(project.data)

  normal.sample <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                                       typesample = "NT")
  tumor.sample <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                                      typesample = "TP")

  group <- data.frame(
    Label = rep( c("Normal","Tumor"), c(length(normal.sample), length(tumor.sample))  ),
    Short.Name = substr(c(normal.sample,tumor.sample),1,12),
    Barcode = c(normal.sample,tumor.sample),
    stringsAsFactors = FALSE
  )


  ## clinical
  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Clinical",
                                  data.type = "Clinical Supplement",
                                  data.format = "BCR Biotab")
  TCGAbiolinks::GDCdownload(query, directory = dir)
  clinical <- TCGAbiolinks::GDCprepare(query, directory = dir)


  if(remove.Raw){
    file.remove( paste(dir,"/",project,sep="",collapse = "")  )
  }

  result <- list(clinical = clinical,
                 expression = expression.df,
                 group = group,
                 Gene.info = gene.information,
                 query = query
  )


}




#' Get TCGA RNA expression (RPKM) data
#'
#' @param tcga.project See ?RTCGA.miRNASeq::miRNASeq. For example: COAD
#' @param log2 If to perform log2 transformation
#'
#' @return
#' @export
#'
#' @examples
#' This function is similar to get.TCGA.miRNAExpression in Biomarker.R file.
#' Data source is illumina hiseq Level 3 RSEM normalized expression data. Data from 2015-11-01 snapshot.
#' The RNAseq gene expression level 3 data contains Reads per Kilobase per Million mapped reads (RPKM)
get.TCGA.RNA.FPKM <- function(tcga.project=NULL, rawCount=FALSE, log2=FALSE){

  if(is.null(tcga.project)){
    warning(?RTCGA.miRNASeq::miRNASeq)
    stop("Please specify a TCGA project. See ?RTCGA.miRNASeq::miRNASeq")
  }




  if (!require(RTCGA))BiocManager::install("RTCGA")
  if (!require(RTCGA.rnaseq))BiocManager::install("RTCGA.rnaseq")
  if (!require(RTCGA.clinical))BiocManager::install("RTCGA.clinical")

  # Expression
  ##-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###
  expression.df = get( paste0(tcga.project,".rnaseq") )

  row.names(expression.df) = expression.df$bcr_patient_barcode
  expression.df = expression.df[,-c(1)]

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




