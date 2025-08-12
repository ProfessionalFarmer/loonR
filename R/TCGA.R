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




#' Download TCGA data through TCGABiolinks
#'
#' @param project.id
#' @param mRNA Default T
#' @param miRNA Default T
#' @param met Default T
#' @param protein Default T
#' @param mutation Default T
#' @param cnv Default T
#' @param clin Default T
#' @param subtype
#' @param dir Default current directory
#'
#' @return
#' @export
#'
#' @examples
download.TCGA.Biolinks = function(project.id = NULL,
                                  mRNA = T, miRNA = T,
                                  met =T, protein =T, mutation = T,
                                  cnv = T, clin = T, subtype = T, dir = "./"){

  if(is.null(project.id)){
    stop("Pls run TCGAbiolinks::getGDCprojects()$project_id")
  }
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(DT)

  #Simple Nucleotide Variation  √
  #Sequencing Reads
  #Biospecimen
  #Clinical
  #Copy Number Variation  profiling √
  #DNA Methylation √
  #Proteome Profiling √
  #Somatic Structural Variation  √
  #Structural Variation ×
  #miRNA Expression Quantification  √
  #GDCquery( project = c("TCGA-LGG") )

  data = list()

  if(mRNA){
    query <- GDCquery(
           project = project.id,
           data.category = "Transcriptome Profiling",
           data.type = "Gene Expression Quantification",
           workflow.type = "STAR - Counts"
    )
    cat("Downloading mRNA\n\n")
    GDCdownload(query, files.per.chunk = 50)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".mRNA.rdata"))
  }

  if(miRNA){
    query <- GDCquery(
      project = project.id,
      data.category = "Transcriptome Profiling",
      data.type = "miRNA Expression Quantification",
      workflow.type = "BCGSC miRNA Profiling"
    )
    cat("Downloading miRNA\n\n")
    GDCdownload(query, files.per.chunk = 300)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".miRNA.rdata"))
  }

  if(met){
    query <- GDCquery(
      project = project.id,
      data.category = "DNA Methylation",
      data.type = "Methylation Beta Value",
      platform = c("Illumina Human Methylation 450")
    )
    cat("Downloading 450k\n\n")
    GDCdownload(query, files.per.chunk = 10)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".met.450k.rdata"))
  }

  if(protein & project.id!= "TCGA-LAML" & project.id!= "TCGA-THCA"){
    query <- GDCquery(
        project = project.id,
        data.category = "Proteome Profiling",
        data.type = "Protein Expression Quantification"
    )
    cat("Downloading protein\n\n")
    GDCdownload(query, files.per.chunk = 300)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".protein.rppa.rdata"))
  }

  if(mutation){
    query <- GDCquery(
      project = project.id,
      data.category = "Simple Nucleotide Variation",
      access = "open",
      data.type = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    )
    cat("Downloading mutation\n\n")
    GDCdownload(query, files.per.chunk = 50)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".mutation.maf.rdata"))
  }

  if(cnv){
    query <- GDCquery(
      project = project.id,
      data.category = "Copy Number Variation",
      data.type = "Masked Copy Number Segment"
    )
    cat("Downloading CNV\n\n")
    GDCdownload(query, files.per.chunk = 50)
    GDCprepare(query, save = T, save.filename = paste0(dir, "/", project.id,".cnv.seg.rdata"))

  }

  if(clin){
    query <- GDCquery(
      project = project.id,
      data.category = "Clinical",
      data.type = "Clinical Supplement",
      data.format = "BCR Biotab"
    )
    cat("Downloading clinical information\n\n")
    GDCdownload(query)
    clinical.BCRtab.all <- GDCprepare(query)
    save(clinical.BCRtab.all, file = paste0(dir, "/", project.id, ".clin.brc.rdata") )
  }
  subtype.available = c("TCGA-ACC", "TCGA-BRCA", "TCGA-BLCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SKCM", "TCGA-SARC", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")

  if(subtype & project.id%in%subtype.available){

    cat("Downloading subtype information\n\n")
    subtype <- TCGAquery_subtype(tumor = stringr::str_to_lower(stringr::str_remove(project.id, "TCGA-")))
    save(subtype, file = paste0(dir, "/", project.id,".subtype.rdata") )
  }


}

get.pancancer.subtype = function(){
  TCGAbiolinks::PanCancerAtlas_subtypes()
}

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






