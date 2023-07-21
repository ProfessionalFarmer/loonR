#' Get promoters by GenomicFeatures package. Region needed
#'
#' @param upstream Default 2000
#' @param downstream  Default 500
#' @param ann Default Ensembl
#' @param ref.genome Default hg38 (for refGene)
#' @param ens.release Default 99 (for Ensembl)
#'
#' @return
#' @export
#'
#' @examples
getPromoterRegions <- function(upstream=2000, downstream=500, ann = "Ensembl", ref.genome = "hg38", ens.release = 99, addChr = TRUE ){

  if(!require("GenomicFeatures")){
    BiocManager::install("GenomicFeatures")
  }
  if(!require("RMariaDB")){
    BiocManager::install("RMariaDB")
  }


  if (ann == "Ensembl"){
    # https://www.ensembl.org/info/data/mysql.html
    ensembl.ann <- makeTxDbFromEnsembl(organism = "Homo sapiens",
                                       release = ens.release,
                                       server = "ensembldb.ensembl.org")
    genes <- genes(ensembl.ann)
    if (addChr){
      if (!requireNamespace("diffloop", quietly = TRUE)){ BiocManager::install("diffloop") }
      diffloop::addchr(genes)
    }


  }else if(ann == "refGene"){
    refgene <- makeTxDbFromUCSC(genome = ref.genome, tablename="refGene")
    genes <- genes(refgene)

  }

  promoter <- promoters(genes, upstream = upstream, downstream = downstream, use.names=TRUE)
  promoter

}







#' Load methylation arrary from directory
#'
#' @param dir
#' @param arraytype EPIC or 450K
#' @param combat If run combat to remove batch effect
#' @param batchname Default c("Slide")
#'
#' @return
#' @export
#'
#' @examples
#' Reference: https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
loadMethylationArraryData <- function(dir, arraytype = 'EPIC', combat=FALSE, batchname=c("Slide")){
  library(ChAMP)
  myLoad <- champ.load(dir, arraytype = arraytype)

  # champ.QC() function and QC.GUI() function would draw some plots for user to easily check their data’s quality.
  champ.QC()
  myNorm <-champ.norm(arraytype = arraytype, cores=50, method="BMIQ")
  champ.SVD()

  # 5.6 Batch Effect Correction
  if(combat){
    myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd, batchname= batchname )
    myCombat
  }else{
    myNorm
  }

}





#' ChAMP QC Pipeline for beta matrix
#'
#' @param Sample.beta.df Column is sample
#' @param Sample.Group
#' @param Slide
#' @param arraytype 450K or EPIC
#' @param CpG.GUI If use CpG.GUI
#' @param QC.GUI
#' @param combat If perform combat to remove batch effect
#' @param batchname Batch variable name. Default is c("Slide")
#'
#' @return
#' @export
#'
#' @examples
#'
#' library(ChAMP)
#' testDir=system.file("extdata",package="ChAMPdata")
#' myImport <- champ.import(testDir)
#'
#' beta.df = myImport$beta
#' group = myImport$pd$Sample_Group
#' myLoad = loonR::ChAMP_QC_Pipeline_Frome_Beta_Value(Sample.beta.df=beta.df, Sample.Group=group)
#'
#' ########### Step 1
#' res = loonR::ChAMP_QC_Pipeline_Frome_Beta_Value
#' myLoad = res$myLoad
#'
#' ########### Step 2
#' res = res$followAnalysisImputeQCNormCombat(arraytype="450K)
#' head(res$myNorm)
ChAMP_QC_Pipeline_Frome_Beta_Value <- function(Sample.beta.df=NULL, Sample.Group='', Slide = '', arraytype=c("450K","EPIC"), CpG.GUI=FALSE, QC.GUI=FALSE, combat=FALSE, batchname=c("Slide") ){

  library(ChAMP)
  # https://www.bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
  res = list()

  if(is.null(Sample.beta.df)){
    stop("Please input sample beta data.frame")
  }

  arraytype <- match.arg(arraytype)

  pd = data.frame(
    Sample_Name = colnames(Sample.beta.df),
    rownames=colnames(Sample.beta.df),
    Sample_Group = Sample.Group,
    Slide = Slide
  )

  res$rawBeta = Sample.beta.df
  res$pd = pd

  # 5.2 filter 过滤SNP和XY上的
  warning("Start filtering \n")
  myLoad = champ.filter(
    beta = Sample.beta.df,
    pd = pd,
    filterXY=TRUE,
    filterMultiHit=TRUE,
    arraytype=arraytype,
    filterSNPs=TRUE,
    filterNoCG=TRUE,
    fixOutlier=FALSE,
    autoimpute=TRUE
  )

  followAnalysisImputeQCNormCombat <-  function(arraytype = NA){

    res = list()

    if(is.na(arraytype)){
      stop("Pls set arraytype in this function\n")
    }

    #if NAs are still existing
    tmp.myLoad = champ.impute()
    if(ncol(tmp.myLoad$beta)!=ncol(myLoad$beta)){
      stop("Some samples are filtered out, pls make sure not samples were filtered")
    }else{
      myLoad = tmp.myLoad
      rm(tmp.myLoad)
    }

    warning("Start imputation \n")
    if(sum(is.na(myLoad$beta))!=0){
      warning("LoonR:: Because there still has NA values after filtering, we need to use champ.impute")
      row.to.remove = rowSums(is.na(myLoad$beta)) > (ncol(myLoad$beta)/2)
      warning("LoonR:: ", sum(row.to.remove), " Probes have 0.5 or above NA will be removed")

      myLoad <- champ.impute(
        beta=as.matrix(myLoad$beta[!row.to.remove, ]),
        pd=myLoad$pd,
        method="Combine",
        k=5,
        ProbeCutoff=0.2,
        SampleCutoff=0.2
      )
    }
    res$myLoad = myLoad

    # 5.3 show CpG distribution, 也可以显示DMP的分布
    # This is a useful function to demonstrate the distribution of your CpG list. If you get a DMP list, you may use this function to check your DMP’s distribution on chromosome.
    if(CpG.GUI){
      CpG.GUI(CpG=rownames(myLoad$beta),arraytype=arraytype)
    }

    # QC check
    # champ.QC() function and QC.GUI() function would draw some plots for user to easily check their data’s quality.
    warning("Start QC \n")
    champ.QC()
    if(QC.GUI){
      QC.GUI(beta=myLoad$beta,arraytype=arraytype)
    }

    # 5.4 Normalization
    warning("Start normalization \n")
    myNorm <- champ.norm(beta=as.matrix(myLoad$beta), arraytype=arraytype, cores=1, method="BMIQ")
    res$myNorm = myNorm

    # QC check again
    warning("Start QC \n")
    champ.QC(beta = myNorm)
    if(QC.GUI){
      QC.GUI(beta=myNorm,arraytype=arraytype)
    }

    # 5.5
    warning("Start SVD \n")
    champ.SVD(beta=data.frame(myNorm, check.names = F),pd=myLoad$pd)

    # 5.6 Batch Effect Correction
    if(combat){
      myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd, batchname= batchname )
      res$myCombat = myCombat
    }
    res$probe.features = probe.features
    res
  }
  res = list(myLoad = myLoad, followAnalysisImputeQCNormCombat = followAnalysisImputeQCNormCombat)
  res
}

