
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
  if (project==NA){
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
                 Gene.info = gene.information
  )


}





