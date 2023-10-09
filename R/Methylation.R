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
#' @param cores
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
#'
#' ########### Step 1
#' res = loonR::ChAMP_QC_Pipeline_Frome_Beta_Value(Sample.beta.df=beta.df, Sample.Group=group)
#' myLoad = res$myLoad
#'
#' ########### Step 2
#' res = res$followAnalysisImputeQCNormCombat(arraytype="450K")
#' head(res$myNorm)
ChAMP_QC_Pipeline_Frome_Beta_Value <- function(Sample.beta.df=NULL, Sample.Group='', Slide = '', arraytype=c("450K","EPIC"), CpG.GUI=FALSE, QC.GUI=FALSE, combat=FALSE, batchname=c("Slide"), cores = 50 ){

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

  followAnalysisImputeQCNormCombat <-  function(arraytype = NA, CpG.GUI = FALSE, cores = 1){

    res = list()

    if(is.na(arraytype)){
      stop("Pls set arraytype in this function\n")
    }

    #if NAs are still existing
    tmp.myLoad = champ.impute(beta = as.matrix(myLoad$beta), pd = myLoad$pd, SampleCutoff = 0.2)
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
    champ.QC(beta = myLoad$beta, resultsDir = "1stQC")
    if(QC.GUI){
      QC.GUI(beta = myLoad$bet, arraytype = arraytype)
    }

    # 5.4 Normalization
    warning("Start normalization \n")
    myNorm <- champ.norm(beta = myLoad$beta, arraytype=arraytype, cores=cores, method="BMIQ")
    res$myNorm = myNorm

    # QC check again
    warning("Start QC \n")
    champ.QC(beta = myNorm, resultsDir = "2ndQC.Norm")
    if(QC.GUI){
      QC.GUI(beta = myNorm, arraytype = arraytype)
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


#' RNA and methylation analysis
#'
#' @param rna.df row name should be ensembl ID
#' @param met.df row name should be probeID
#' @param group normal 在前, cancer 在后
#' @param genome hg19 or hg38
#' @param met.platform 450K or EPIC
#' @param beta.dif Default 0.2. Beta difference
#' @param dir.group1 It can be "hypo" which means the probes are hypomethylated in group1; "hyper" which means the probes are hypermethylated in group1
#' @param specify.probes Specify methylation probe instead of ELMER
#' @param TF.motif Default FALSE. If perform TF and motif analysi
#' @param run.diff If run methylation analysis
#'
#' @return
#' @export
#'
#' @examples
ELMER_RNA_Met_analysis <- function(rna.df = NULL, met.df = NULL, group = NULL, genome = "hg19", met.platform = "450K", beta.dif = 0.2, dir.group1 = "hyper",
                                   specify.probes = NULL, TF.motif = FALSE, run.diff = TRUE){

  if(!require(sesameData)){
    pacman::p_load_gh("zwdzwd/sesameData")
  }
  if(!require(ELMER.data)){
    pacman::p_load_gh("tiagochst/ELMER.data")
  }
  if(!require(ELMER)){
    pacman::p_load_gh("tiagochst/ELMER")
  }

  pacman::p_load(ELMER.data, ELMER)

  if( is.null(rna.df) | is.null(met.df) | is.null(group) ){
    stop("Exp or met data.frame, or group should not be NA")
  }

  if( ! loonR::AllEqual(colnames(rna.df), colnames(met.df) ) ){
    stop("RNA and met data.frame should be the same")
  }

  if(ncol(rna.df) != ncol(met.df) ){
    stop("Exp and met should have the same sample")
  }

  if(ncol(rna.df) != length(group) ){
    stop("Exp, met, group should have the same number of samples")
  }


  # https://www.bioconductor.org/packages/release/bioc/vignettes/ELMER/inst/doc/index.html
  ####################### get distal probes that are 2kb away from TSS on chromosome 1
  distal.probes <- get.feature.probe(
    genome = genome,
    met.platform = met.platform,
    rm.chr = paste0("chr",c("X","Y"))
  )

  sample.info = data.frame(row.names = colnames(rna.df),
                           primary = colnames(rna.df),
                           sample = colnames(rna.df),
                           group = group)


  assay <- c(
    rep("DNA methylation", ncol(met.df)),
    rep("Gene expression", ncol(rna.df))
  )
  primary <- c(colnames(met.df),colnames(rna.df))
  colname <- c(colnames(met.df),colnames(rna.df))
  sampleMap <- data.frame(assay,primary,colname)
  rm(assay, primary, colname)

  ####################### Creation of a MAE object
  mae <- createMAE(
    exp = rna.df, # ensembl ID
    met = met.df,
    colData = sample.info,
    sampleMap = sampleMap,
    save = TRUE,
    linearize.exp = TRUE,
    save.filename = "ELMERresult/mae.rda",
    filter.probes = distal.probes,
    met.platform = met.platform,
    genome = genome,
    TCGA = FALSE
  )
  res = list(mae = mae)


  ############# Identifying differentially methylated probes
  if(run.diff){
  library(SummarizedExperiment)
  library(DT)
  # A group from group.col. ELMER will run group1 vs group2.
  # That means, if direction is hyper, get probes hypermethylated in group 1 compared to group 2.
  sig.diff <- get.diff.meth(
    data = mae,
    group.col = "group",
    group1 = unique(group)[1], # group2 is normal ELMER run group1 vs group2. Hypo means hypomethylated probes in group 1
    group2 = unique(group)[2], # group2 is cancer.
    minSubgroupFrac = 1, # if supervised mode set to 1
    sig.dif = beta.dif,
    diff.dir = dir.group1, # Search for hypermethylated probes in group 1 (hypo in group2)
    cores = 50,
    dir.out ="ELMERresult",
    pvalue = 0.05
  )

  ############# result
  res$sig.diff = sig.diff
  }



  ############ Identifying putative probe-gene pairs
  if(is.null(specify.probes)){
    nearGenes <- GetNearGenes(
      data = mae,
      probes = sig.diff$probe,
      numFlankingGenes = 20
    ) # 10 upstream and 10 dowstream genes
  }else{
    nearGenes <- GetNearGenes(
      data = mae,
      probes = specify.probes,
      numFlankingGenes = 20
    ) # 10 u
    res$specified.probes = specify.probes
  }



  # https://www.bioconductor.org/packages/release/bioc/vignettes/ELMER/inst/doc/analysis_get_pair.html
  probe.gene.pair <- get.pair(
    data = mae,
    group.col = "group",
    group1 = unique(group)[1], # normal
    group2 = unique(group)[2], # cancer
    nearGenes = nearGenes,
    # If diff.dir is "hypo, U will be the group 1 and M the group2.
    # If diff.dir is "hyper" M group will be the group1 and U the group2.
    mode = "supervised",
    # It can be "hypo" which means the probes are hypomethylated in group1;
    # "hyper" which means the probes are hypermethylated in group1
    diff.dir = dir.group1,
    permu.dir = "ELMERresult/permu",
    permu.size = 100000, # Please set to 100000 to get significant results
    raw.pvalue = 0.05,
    Pe = 0.001, # Please set to 0.001 to get significant results
    # when preAssociationProbeFiltering: U > 0.3, M < 0.3
    filter.probes = TRUE, # See preAssociationProbeFiltering function
    filter.percentage = 0.05,
    filter.portion = 0.3,
    dir.out = "ELMERresult",
    cores = 1,
    label = dir.group1
  )

  res$probe.gene.pair = probe.gene.pair

  if(TF.motif){
  # 3.4 - Motif enrichment analysis on the selected probes
  enriched.motif <- get.enriched.motif(
    data = mae,
    probes = probe.gene.pair$Probe,
    dir.out = "ELMERresult",
    label = dir.group1,
    min.incidence = 10,
    lower.OR = 1.1
  )
  res$enriched.motif = enriched.motif


  # 3.5 - Identifying regulatory TFs
  TF <- get.TFs(
    data = mae,
    group.col = "group",
    group1 = unique(group)[1],
    group2 = unique(group)[2],
    mode = "supervised",
    enriched.motif = enriched.motif,
    dir.out = "ELMERresult",
    cores = 30,
    label = ELMERresult,
    save.plots = TRUE
  )
  res$TF = TF
  }


  res$PlotTutorial = "https://www.bioconductor.org/packages/release/bioc/vignettes/ELMER/inst/doc/plots_scatter.html"

  res
}



