


#' Get promoters by GenomicFeatures package
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

  library("GenomicFeatures")

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

#' Title
#'
#' @param dir
#' @param arraytype EPIC or 450K
#'
#' @return
#' @export
#'
#' @examples
#' Reference: https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
#'
loadMethylationArraryData <- function(dir, arraytype = 'EPIC'){
  library(ChAMP)
  myLoad <- champ.load(dir, arraytype = arraytype)

  # champ.QC() function and QC.GUI() function would draw some plots for user to easily check their data’s quality.
  champ.QC()
  myNorm <-champ.norm(arraytype = arraytype, cores=50)
  champ.SVD()


}




