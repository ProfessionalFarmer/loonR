

#' Perform over-representation test by cluster profiler
#'
#' @param gene A vector. Must gene symbol
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.OverRepresentationTest <-  function(gene, minGSSize = 10, qvalue = 0.05,
                                                    GO=TRUE, KEGG =TRUE, MSigDb=TRUE, Hallmark = TRUE){

  library(clusterProfiler)
  library(org.Hs.eg.db)
  result = list()

  if(GO){
    ego <- enrichGO(gene          = loonR::id_mapping(gene, key="SYMBOL", column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique(),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = qvalue,
                    minGSSize     = minGSSize,
                    readable      = TRUE)
    result$BP = ego
  }

  if(KEGG){
    kk <- enrichKEGG(gene         = loonR::id_mapping(gene, key="SYMBOL", column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique(),
                     organism     = 'hsa',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = qvalue,
                     minGSSize     = minGSSize)
    result$KEGG = kk
  }

  if(MSigDb){

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.symbols.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.symbols.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)
    kk <- enricher(gene,
             TERM2GENE = c2c5,
             minGSSize = minGSSize,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = qvalue)
    result$C2C5 <- kk
  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    kk <- enricher(gene,
                   TERM2GENE = H,
                   minGSSize = minGSSize,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = qvalue)
    result$Hallmark <- kk
  }
  # library(enrichplot)
  # dotplot(gene.gsea, showCategory=30)
  # emapplot(gene.gsea)
  result
}

#' Perform GSEA analysis by cluster profiler
#'
#' @param phenotype log2FC
#' @param geneNames gene symbols. Must be symbol
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE

#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.GSEA <- function(phenotype, geneNames, minGSSize = 10, qvalue = 0.05,
                                 GO=TRUE, KEGG =TRUE, MSigDb=TRUE, Hallmark = TRUE) {

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)

  result = list()

  geneList <- phenotype
  names(geneList) <- geneNames
  geneList4MsigDB <- geneList

  id.map <- loonR::id_mapping(geneNames, key="SYMBOL", column="ENTREZID")

  geneList <- geneList[names(geneList) %in% id.map$SYMBOL]
  names(geneList) <- id.map$ENTREZID[match(names(geneList), id.map$SYMBOL) ]

  geneList <- sort(geneList, decreasing = TRUE)
  geneList4MsigDB <- sort(geneList4MsigDB, decreasing = TRUE)


  if(GO){
    ego3 <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  nPerm        = 1000,
                  minGSSize    = minGSSize,
                  pAdjustMethod = "BH",
                  pvalueCutoff = qvalue,
                  verbose      = FALSE)

    result$GO = ego3
  }

  if(KEGG){

    kk2 <- gseKEGG(geneList     = geneList,
                   organism     = 'hsa',
                   nPerm        = 1000,
                   minGSSize    = minGSSize,
                   pAdjustMethod = "BH",
                   pvalueCutoff = qvalue,
                   verbose      = FALSE)

    result$KEGG = kk2

  }

  if(MSigDb){

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.entrez.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.entrez.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)
    kk <- GSEA(geneList     = geneList4MsigDB,
               TERM2GENE = c2c5,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               pvalueCutoff = qvalue,
               verbose      = FALSE)
    result$C2C5 <- kk

  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.entrez.gmt")

    kk <- GSEA(geneList     = geneList4MsigDB,
               TERM2GENE = H,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               pvalueCutoff = qvalue,
               verbose      = FALSE)
    result$Hallmark <- kk
  }

  result

}






#' Perform over-representation test across multiple group by cluster profiler
#'
#' @param gene A list including multiple group. Must gene symbol
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.OverRepresentationTest.Compare <- function(gene, minGSSize = 10, qvalue = 0.05,
                                                           GO=FALSE, KEGG =FALSE, MSigDb=TRUE, Hallmark = TRUE){

  gene.entrez.list <- lapply(gene, function(group.gene){
    loonR::id_mapping(group.gene, key="SYMBOL", column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique()
  })


  library(clusterProfiler)
  library(org.Hs.eg.db)
  result = list()


  if(GO){
    ego <- compareCluster(gene.entrez.list,
                          fun="enrichGO",
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = qvalue,
                          minGSSize     = minGSSize,
                          readable      = TRUE )
    result$BP = ego
  }

  if(KEGG){
    kk <- compareCluster(gene.entrez.list,
                         fun="enrichKEGG",
                         organism     = 'hsa',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = qvalue,
                         minGSSize     = minGSSize)
    result$KEGG = kk
  }

  if(MSigDb){

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.symbols.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.symbols.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)
    kk <- compareCluster(gene,
                          fun="enricher",
                          TERM2GENE = c2c5,
                          minGSSize = minGSSize,
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = qvalue)
    result$C2C5 <- kk
  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    kk <- compareCluster(gene,
                         fun="enricher",
                         TERM2GENE = H,
                         minGSSize = minGSSize,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = qvalue)

    result$Hallmark <- kk
  }
  # library(enrichplot)
  # dotplot(gene.gsea, showCategory=30)
  # emapplot(gene.gsea)
  result

  # another
  #f = function(de){ enricher(de, TERM2GENE=m_t2g) }
  #ck <- compareCluster(Entrez ~ MESCC, data=diff.analysis.clean.res, fun=f)
}



#' Perform GSE analysis across multiple group by cluster profiler
#'
#' @param gene A list including multiple group. Must vector names are gene symbol and vector value are log fold chage
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.GSEA.Compare <- function(gene, minGSSize = 10, qvalue = 0.05,
                                        GO=FALSE, KEGG =FALSE, MSigDb=TRUE, Hallmark = TRUE){
  library(dplyr)
  gene.entrez.list <- lapply(gene, function(group.gene){
    # covert symbol to entrezid
    names(group.gene) %<>% loonR::id_mapping(
                      key="SYMBOL",
                      column="ENTREZID") %>% pull(ENTREZID)
    # sort
    group.gene <- sort(group.gene, decreasing = TRUE)
    group.gene
  })


  library(clusterProfiler)
  library(org.Hs.eg.db)

  result = list()


  if(GO){
    f <- function(entrez.sort.by.fc){

      gseGO(geneList     = entrez.sort.by.fc,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "BP",
                    keyType      = "ENTREZID",
                    minGSSize    = minGSSize,
                    maxGSSize    = 500,
                    pvalueCutoff = qvalue,
                    pAdjustMethod = "BH",
                    verbose      = FALSE)

    }

    ego <- compareCluster(gene.entrez.list, fun=f)
    result$BP = ego
  }

  if(KEGG){
    f <- function(entrez.sort.by.fc){

      gseKEGG(geneList     = entrez.sort.by.fc,
              organism     = 'hsa',
              minGSSize    = minGSSize,
              pAdjustMethod = "BH",
              pvalueCutoff = qvalue,
              verbose      = FALSE)
    }

    kk <- compareCluster(gene.entrez.list, fun=f)
    result$KEGG = kk
  }

  if(MSigDb){

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.entrez.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.entrez.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)

    f <- function(entrez.sort.by.fc){
      GSEA(geneList     = entrez.sort.by.fc,
                 TERM2GENE = c2c5,
                 minGSSize = minGSSize,
                 pAdjustMethod = "BH",
                 pvalueCutoff = qvalue,
                 verbose      = FALSE)
    }


    kk <- compareCluster(gene.entrez.list, fun=f)
    result$C2C5 <- kk
  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.entrez.gmt")

    f <- function(entrez.sort.by.fc){
      GSEA(geneList     = entrez.sort.by.fc,
           TERM2GENE = H,
           minGSSize = minGSSize,
           pAdjustMethod = "BH",
           pvalueCutoff = qvalue,
           verbose      = FALSE)
    }
    kk <- compareCluster(gene.entrez.list, fun=f)

    result$Hallmark <- kk
  }

  result

  # another
  #f = function(de){ enricher(de, TERM2GENE=m_t2g) }
  #ck <- compareCluster(Entrez ~ MESCC, data=diff.analysis.clean.res, fun=f)
}




#' Perform analysis using user custom gene sets
#'
#' @param g GSE: values with genes symbols, OVA: symbols
#' @param CustomGS ClusterProfiler::read.gmt() term to symbols
#' @param gse Default FALSE, gene set enrichment analysis
#' @param ova Default FALSE, overrepresentation test
#' @param minGSSize 10
#' @param qvalue 0.05
#'
#' @return
#' @export
#'
#' @examples
clusterProfilter.GSEA.ORA.customGS <- function(g, CustomGS = NULL, gse=FALSE, ova=FALSE, minGSSize = 10, qvalue = 0.05){

  if(is.null(customGS)){
    stop("pls set gset")
  }
  if( (!gse) && (!ova)){
    stop("pls specify gene set enrichment analysis or overrepresentation test")
  }

  if(gse){
    g = sort(g, decreasing = TRUE)

    res <- GSEA(geneList     = g,
         TERM2GENE = customGS,
         minGSSize = minGSSize,
         pAdjustMethod = "BH",
         pvalueCutoff = qvalue,
         verbose      = FALSE)

  }else if(ova){
    res <- enricher(g,
                    TERM2GENE = customGS,
                    minGSSize = minGSSize,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = qvalue)

  }

  res

}




#' Perform GSE analysis across multiple group by cluster profiler and user customed gene sets
#'
#' @param gene
#' @param customGS ClusterProfiler::read.gmt()  term to symbols
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param gse Default FALSE, gene set enrichment analysis
#' @param ova Default FALSE, overrepresentation test
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.GSEA.ORA.customGS.Compare <- function(gene, customGS=NULL, minGSSize = 10, qvalue = 0.05, gse=FALSE, ova=FALSE){

  if(is.null(customGS)){
    stop("pls set gset")
  }

  library(dplyr)
  if( (!gse) && (!ova)){
    stop("pls specify gene set enrichment analysis or overrepresentation test")
  }

  library(clusterProfiler)
  library(org.Hs.eg.db)

  if(gse){
    g.list <- lapply(gene, function(group.gene){
      # # covert symbol to entrezid
      # names(group.gene) %<>% loonR::id_mapping(
      #   key="SYMBOL",
      #   column="ENTREZID") %>% pull(ENTREZID)
      # sort
      group.gene <- sort(group.gene, decreasing = TRUE)
      group.gene
    })

    f <- function(symbol.sorted.byFC.vector){
      GSEA(geneList     = symbol.sorted.byFC.vector,
           TERM2GENE = customGS,
           minGSSize = minGSSize,
           pAdjustMethod = "BH",
           pvalueCutoff = qvalue,
           verbose      = FALSE)
    }

    res <- compareCluster(g.list, fun=f)

  }else if(ova){

    res <-compareCluster(g.list,
                   fun="enricher",
                   TERM2GENE = customGS,
                   minGSSize = minGSSize,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = qvalue)

  }

  res
}

