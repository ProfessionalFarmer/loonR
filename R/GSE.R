

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
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

    result$GO = ego3
  }

  if(KEGG){

    kk2 <- gseKEGG(geneList     = geneList,
                   organism     = 'hsa',
                   nPerm        = 1000,
                   minGSSize    = minGSSize,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)

    result$GO = kk2

  }

  if(MSigDb){

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.symbols.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.symbols.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)
    kk <- GSEA(geneList     = geneList4MsigDB,
               TERM2GENE = c2c5,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               verbose      = FALSE)
    result$C2C5 <- kk

  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    kk <- GSEA(geneList     = geneList4MsigDB,
               TERM2GENE = H,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
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
}


