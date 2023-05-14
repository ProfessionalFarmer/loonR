#' Get sample barcode and TCGA project map information
#'
#'
#' @param project A vector. TCGA-ACC, TCGA-COAD. See TCGAbiolinks::getGDCprojects()$project_id
#' @return A data.frame
#'  X1 proj samples
#'	1	1	TCGA-ACC	TCGA-OR-A5JZ
#'	2	1	TCGA-ACC	TCGA-OR-A5JK
#' @export
#'
#' @examples
get.TCGA.sample.project.infomap <- function(project=NA){

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

  library(TCGAbiolinks)
  library(SummarizedExperiment)

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
  gene.information <- SummarizedExperiment::rowRanges(project.data)
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

  library(dplyr)
  if(any( ! substr(patient_ids,14,15) %in% c("01","11")  )){
    warning("It seems some samples are not primary tumor or normal samples\nWe have removed\n")
    patient_ids = patient_ids[ substr(patient_ids,14,15) %in% c("01","11") ]
    expression.df = expression.df[patient_ids,]
  }


  if(log2){
    expression.df = t(log2(expression.df+1))
  }else{
    expression.df = t( expression.df )

  }

  # Clinical
  ##-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###-##-###-###-##-###
  # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  group_list <- ifelse(substr(patient_ids,14,15)=='01','Tumor','Normal')

  ind = order(group_list)
  group_list = group_list[ind]
  expression.df = expression.df[,ind]
  rm(ind)

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


#' Get TCGA project clinical information
#'
#' @param project ESCA. TCGAbiolinks::getGDCprojects()$project_id without "TCGA_"
#' @param subtype if download subtype information
#' @param microsatellite if download microsatellite data
#'
#' @return
#' @export
#'
#' @examples
get.TCGA.Project.Clinical <- function(project, subtype=FASLE, microsatellite=FALSE){

  library(TCGAbiolinks)

  res = list()

  # 1 Clinical indexed data
  clinical = TCGAbiolinks::GDCquery_clinic(project = paste0(c("TCGA-",project), sep="", collapse = "" ), type = "clinical")
  res$clinical = clinical

  # 2 BCR Biotab: tsv files parsed from XML files
  query <- TCGAbiolinks::GDCquery(project = paste0(c("TCGA-",project), sep="", collapse = "" ),
                    data.category = "Clinical",
                    data.type = "Clinical Supplement",
                    data.format = "BCR Biotab")
  GDCdownload(query)
  clinical.BCRtab.all <- GDCprepare(query)
  res$clinical.BCRtab.all = clinical.BCRtab.all



  # 3 MSI-Mono-Dinucleotide Assay is performed to test a panel of four mononucleotide repeat loci
  if(microsatellite){
    query <- GDCquery(project = paste0(c("TCGA-",project), sep="", collapse = "" ),
                      data.category = "Other",
                      legacy = TRUE,
                      access = "open",
                      data.type = "Auxiliary test"  )
    GDCdownload(query)
    msi_results <- GDCprepare_clinic(query, "msi")

    res$msi = msi_results
  }

  if(subtype){
    subtype = TCGAbiolinks::TCGAquery_subtype(tumor = project)
    res$subtype.info = subtype
  }

  res

}


#' Download 450K methylation data
#'
#' @param project e.g. TCGA-COAD. See getGDCprojects()$project_id
#'
#' @return
#' @export
#'
#' @examples
tcgabiolinks.get.450k.methylation <- function(project){

  library(TCGAbiolinks)
  library(SummarizedExperiment)

  query <- GDCquery(project = project,
                    data.category = "DNA Methylation",
                    platform = "Illumina Human Methylation 450"
  )
  GDCdownload(query)
  met.dat <- GDCprepare(
    query = query,
    summarizedExperiment = TRUE
  )

  clin = data.frame( colData(met.dat) )

  met.dat <- assay(met.dat)

  list(clinical=clin, data=met.dat)

}



TCGA.ids = c("TCGA-3X-AAV9-01", "TCGA-3X-AAVA-01", "TCGA-3X-AAVE-01", "TCGA-4G-AAZO-01", "TCGA-4G-AAZT-01", "TCGA-W5-AA2G-01", "TCGA-W5-AA2I-01", "TCGA-W5-AA2I-11", "TCGA-W5-AA2O-01", "TCGA-W5-AA2Q-01", "TCGA-W5-AA2Q-11", "TCGA-W5-AA2R-01", "TCGA-W5-AA2R-11")


#' Extract TCGA sample types
#'
#' @param TCGA.ids
#'
#' @return
#' @export
#'
#' @examples
#' ids = c("TCGA-FV-A495-01", "TCGA-G3-A3CH-11", "TCGA-G3-A3CH-01", "TCGA-DD-A3A6-11", "TCGA-DD-A3A6-01", "TCGA-DD-A11A-11", "TCGA-DD-A11A-01", "TCGA-CC-A3MB-01", "TCGA-BC-A3KF-01")
#' loonR::extractTCGASampleType(ids)
extractTCGASampleType <- function(TCGA.ids){

  library(dplyr)
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

  code.table = data.table::data.table(
        Code = c("01","02","03","04","05","06","07",
                 "08","09","10","11","12","13",
                 "14","15","16","20","40","50","60",
                 "61","99"),
        Definition = c("Primary Solid Tumor",
                       "Recurrent Solid Tumor",
                       "Primary Blood Derived Cancer - Peripheral Blood",
                       "Recurrent Blood Derived Cancer - Bone Marrow",
                       "Additional - New Primary","Metastatic",
                       "Additional Metastatic","Human Tumor Original Cells",
                       "Primary Blood Derived Cancer - Bone Marrow","Blood Derived Normal",
                       "Solid Tissue Normal","Buccal Cell Normal",
                       "EBV Immortalized Normal","Bone Marrow Normal",
                       "sample type 15","sample type 16",
                       "Control Analyte",
                       "Recurrent Blood Derived Cancer - Peripheral Blood","Cell Lines",
                       "Primary Xenograft Tissue",
                       "Cell Line Derived Xenograft Tissue","sample type 99"),
                 Short.Letter.Code = c("TP","TR","TB","TRBM","TAP","TM",
                       "TAM","THOC","TBM","NB","NT","NBC",
                       "NEBV","NBM","15SH","16SH","CELLC",
                       "TRB","CELL","XP","XCL","99SH")
               )

   code = substr(TCGA.ids, 14, 15)
   short.names = substr(TCGA.ids, 1, 12)

   Sample.table = data.frame(ID=TCGA.ids, Short=short.names, Code=code)

   g.res = loonR::vector.group( v = TCGA.ids, g = code )

   names(g.res) = code.table$Definition[match(names(g.res), code.table$Code)]

   if(sum("Solid Tissue Normal" %in%  names(g.res) )==1){
     paired.samples = intersect(
       substr(g.res$`Solid Tissue Normal`,1,12),
       substr(g.res$`Primary Solid Tumor`,1,12)
       )
     paired.normal.sample = g.res$`Solid Tissue Normal`[match(paired.samples, substr(g.res$`Solid Tissue Normal`,1,12))]
     paired.tumor.sample  = g.res$`Primary Solid Tumor`[match(paired.samples, substr(g.res$`Primary Solid Tumor`,1,12))]

     paried.sample.df = data.frame(Normal = paired.normal.sample,
                                Tumor = paired.tumor.sample,
                                Short = paired.samples)

   }else{
     paried.sample.df = "No paired samples"
   }

   res = list(Sample.table = Sample.table,
              Grouped = g.res,
              Code.description = code.table,
              Paired.samples = paried.sample.df)
   res

}

#' Get shorter TCGA ID
#'
#' @param ids
#'
#' @return
#' @export
#'
#' @examples
tcgaShorterBarcode <- function(ids){
  substr(ids, 1, 16)
}






