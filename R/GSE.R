#' Perform over-representation test by cluster profiler
#'
#' @param gene A vector. Must gene symbol
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE
#' @param exp.gene.type RNA expression ID ENSEMBL. keytypes(org.Hs.eg.db)
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.OverRepresentationTest <-  function(gene, minGSSize = 10, qvalue = 0.05,
                                                    GO=TRUE, KEGG =TRUE, MSigDb=TRUE, Hallmark = TRUE, exp.gene.type="ENSEMBL"){

  library(clusterProfiler)
  library(org.Hs.eg.db)
  result = list()

  if(GO){
    ego <- enrichGO(gene          = loonR::id_mapping(gene, key=exp.gene.type, column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique(),
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
    kk <- enrichKEGG(gene         = loonR::id_mapping(gene, key=exp.gene.type, column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique(),
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
    kk <- enricher(loonR::id_mapping(gene, key=exp.gene.type, column="SYMBOL") %>% filter(!is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique(),
             TERM2GENE = c2c5,
             minGSSize = minGSSize,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = qvalue)
    result$C2C5 <- kk
  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    kk <- enricher(loonR::id_mapping(gene, key=exp.gene.type, column="SYMBOL") %>% filter(!is.na(SYMBOL)) %>% pull(SYMBOL) %>% unique(),
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
#' @param exp.gene.type RNA expression ID ENSEMBL. keytypes(org.Hs.eg.db)
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.GSEA <- function(phenotype, geneNames, minGSSize = 10, qvalue = 0.05,
                                 GO=TRUE, KEGG =TRUE, MSigDb=TRUE, Hallmark = TRUE, exp.gene.type = "SYMBOL") {

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)

  result = list()

  geneList <- phenotype
  names(geneList) <- geneNames
  geneList4MsigDB <- geneList

  id.map <- loonR::id_mapping(geneNames, key=exp.gene.type, column=c("ENTREZID", "SYMBOL") )

  geneList <- geneList[names(geneList) %in% as.character(unlist(id.map[,exp.gene.type])) ]
  names(geneList) <- id.map$ENTREZID[match(names(geneList), as.character(unlist(id.map[,exp.gene.type]))) ]
  geneList <- sort(geneList, decreasing = TRUE)


  # geneList4MsigDB <- geneList4MsigDB[names(geneList4MsigDB) %in% as.character(unlist(id.map[,exp.gene.type])) ]
  # names(geneList4MsigDB) <- id.map$SYMBOL[match(names(geneList4MsigDB), as.character(unlist(id.map[,exp.gene.type]))) ]
  # geneList4MsigDB <- sort(geneList4MsigDB, decreasing = TRUE)


  if(GO){
    ego3 <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  nPerm        = 2000,
                  minGSSize    = minGSSize,
                  pAdjustMethod = "BH",
                  pvalueCutoff = qvalue,
                  verbose      = FALSE)

    result$GO = ego3
  }

  if(KEGG){

    kk2 <- gseKEGG(geneList     = geneList,
                   organism     = 'hsa',
                   nPerm        = 2000,
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
    kk <- GSEA(geneList     = geneList,
               TERM2GENE = c2c5,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               nPerm        = 2000,
               pvalueCutoff = qvalue,
               verbose      = FALSE)
    result$C2C5 <- kk

  }

  if(Hallmark){
    H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.entrez.gmt")

    kk <- GSEA(geneList     = geneList,
               TERM2GENE = H,
               minGSSize = minGSSize,
               pAdjustMethod = "BH",
               pvalueCutoff = qvalue,
               nPerm        = 2000,
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
#' @param exp.gene.type RNA expression ID ENSEMBL. keytypes(org.Hs.eg.db)
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.OverRepresentationTest.Compare <- function(gene, minGSSize = 10, qvalue = 0.05,
                                                           GO=FALSE, KEGG =FALSE, MSigDb=TRUE, Hallmark = TRUE, exp.gene.type = "ENSEMBL"){

  gene.entrez.list <- lapply(gene, function(group.gene){
    loonR::id_mapping(group.gene, key=exp.gene.type, column="ENTREZID") %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique()
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

    c5 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c5.all.v7.3.entrez.gmt")
    c2 <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/c2.all.v7.3.entrez.gmt")
    #H <- read.gmt("/data/home2/Zhongxu/R/x86_64-pc-linux-gnu-library/4.0/clusterProfiler/extdata/h.all.v7.4.symbols.gmt")

    c2c5 <- rbind(c2,c5)
    kk <- compareCluster(gene.entrez.list,
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

    kk <- compareCluster(gene.entrez.list,
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
#' @param gene A list including multiple group. Vector names are gene.key and vector value are log fold chage
#' @param minGSSize 10
#' @param qvalue 0.05
#' @param GO TRUE
#' @param KEGG TRUE
#' @param MSigDb TRUE
#' @param Hallmark TRUE
#' @param exp.gene.type key types of input phenotype. Default ENSEMBL, can by keytypes(org.Hs.eg.db)
#'
#' @return
#' @export
#'
#' @examples
ClusterProfiler.GSEA.Compare <- function(gene, minGSSize = 10, qvalue = 0.05, exp.gene.type="ENSEMBL",
                                        GO=FALSE, KEGG =FALSE, MSigDb=TRUE, Hallmark = TRUE){
  library(dplyr)
  gene.entrez.list <- lapply(gene, function(group.gene){
    # covert symbol to entrezid
    names(group.gene) %<>% loonR::id_mapping(
                      key=exp.gene.type,
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
                    nPerm        = 2000,
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
              nPerm        = 2000,
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
                 nPerm        = 2000,
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
           nPerm        = 2000,
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
ClusterProfiler.GSEA.ORA.customGS <- function(g, CustomGS = NULL, gse=FALSE, ova=FALSE, minGSSize = 10, qvalue = 0.05){

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
         nPerm        = 2000,
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
           nPerm        = 2000,
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


#' Perform GSEA across multiple groups (1 vs other)
#'
#' @param rna.df.log
#' @param group
#' @param prefix Default "Group"
#' @param customGS User customed gene set. qusage::read.gmt. Should be converted to ENTREZID
#' @param exp.gene.type RNA expression ID ENSEMBL. keytypes(org.Hs.eg.db)
#' @param cutoff.log10 Default 3. Minimux or maximum log10 value. Useful when meet inf or draw heatmap
#' @param cal.auc If calcuate AUC value
#' @param permutation Number of permutation
#'
#' @return
#' @export
#'
#' @examples
#' This function will perform GSEA analysis. Default geneset: GOBP, C2, C5, KEGG, HALLMARK
compare_GSE.HTSAnalyzer <- function(rna.df.log, group, prefix="Group", customGS=NULL, exp.gene.type="ENSEMBL", cutoff.log10 = 3, cal.auc = FALSE, permutation=1000){
  library(HTSanalyzeR2)
  library(org.Hs.eg.db)

  if(is.null(customGS)){
     ## generate gene set collection
     GO_BP <- GOGeneSets(species="Hs", ontologies=c("BP"))
     PW_KEGG <- KeggGeneSets(species="Hs")
     MSig_C2 <- MSigDBGeneSets(collection = "C2", species = "Hs")
     MSig_C5 <- MSigDBGeneSets(collection = "C5", species = "Hs")
     MSig_H <- MSigDBGeneSets(collection = "H", species = "Hs")
     ## combine all needed gene set collections into a named list for further analysis
     ListGSC <- list(GO_BP=GO_BP, PW_KEGG=PW_KEGG, MSig_C2=MSig_C2, MSig_C5=MSig_C5, MSig_H=MSig_H)

     # paramter
     nPermutations = permutation
     minGeneSetSize = 10
   }else{
     # library(qusage)
     # CustomGS <- read.gmt("/data/home2/Zhongxu/work/baidu2/raw/geneset.gmt")
     # # convert gene symbol to id
     # CustomGS <- lapply(CustomGS,function(x) {
     #      tmp <- mapIds(org.Hs.eg.db, keys = x,
     #               keytype = "SYMBOL", column = "ENTREZID")
     #      return(tmp[!is.na(tmp)])
     #   }
     # )
     ListGSC <- list(customGS=customGS)

     # parameter
     nPermutations = permutation
     minGeneSetSize = 5
   }


  function.analysis.res <- lapply(unique(group), function(x){

   print(paste("Now, ", prefix, x))

   ###### lapply start
    # differential analysis
   true.ind = which(group==x)
   false.ind = which(group!=x)
   limma.df = rna.df.log[ , c(false.ind, true.ind)]
   limma.diff <- loonR::limma_differential(limma.df, rep(c(FALSE,TRUE),
                                                         c(length(false.ind), length(true.ind))),
                                           cal.AUC = cal.auc)

   ## prepare input for analysis
   phenotype <- limma.diff$logFC
   # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html   -log10(res$logFC) * sign(res$logFC)
   #phenotype <- -log10(limma.diff$adj.P.Val) * sign(limma.diff$logFC)
   names(phenotype) <- row.names(limma.diff)


   phenotype <- sort(phenotype, decreasing = TRUE)

   ## initiate a *GSCA* object
   gsca <- GSCA(listOfGeneSetCollections=ListGSC, geneList=phenotype)

   ## preprocess
   gsca1 <- preprocess(gsca, species="Hs", initialIDs = exp.gene.type, # keytypes(org.Hs.eg.db)  # ENSEMBL   SYMBOL
                       keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                       orderAbsValue = FALSE)


   ## analysis
   if (requireNamespace("doParallel", quietly=TRUE)) {
     doParallel::registerDoParallel(cores=100)
   }

   set.seed(111)
   ## support parallel calculation using multiple cores
   gsca2 <- analyze(gsca1,
                    para=list(pValueCutoff = 1, pAdjustMethod = "BH",
                              nPermutations = nPermutations, minGeneSetSize = minGeneSetSize,
                              exponent = 1),
                    doGSOA = FALSE)  # Here we do gsoa

   ## append gene sets terms
   if(is.null(customGS)){
     gsca3 <- appendGSTerms(gsca2, msigdbGSCs=c("MSig_C2","MSig_H","MSig_C5"), goGSCs=c("GO_BP"), keggGSCs=c("PW_KEGG"))
   }else{
     gsca3 <- appendGSTerms(gsca2)
   }
   gsca3
   #### lapply end
   }
  )
  names(function.analysis.res) <- paste0(prefix,unique(group))
  result = list(rawResult = function.analysis.res)



  # combine gse result
  gse.result = lapply(names(function.analysis.res), function(tmp.group.name){

        tmp.gse.result = function.analysis.res[[tmp.group.name]]@result$GSEA.results

        # append group name and gene set name
        tmp.gse.result.after.append.information = lapply(names(tmp.gse.result), function(tmp.gs.name) {
            y = tmp.gse.result[[tmp.gs.name]]
            y$GroupName <- replicate(nrow(y), tmp.group.name) # x is group name
            y$GSType = tmp.gs.name
            y$GS.ID= row.names(y)
            # when customGS used, Gene.Set.Term is null, so set it.
            if(!is.null(customGS)){y$Gene.Set.Term = row.names(y)}
            return(y)
          })
       names(tmp.gse.result.after.append.information) <- names(tmp.gse.result)

       # combined together
       tmp.res <- do.call("rbind", tmp.gse.result.after.append.information)
       rm(tmp.gse.result.after.append.information)

       tmp.res
  })
  # GSEA result by group
  names(gse.result) = names(function.analysis.res)
  result$gseResultByGroup = gse.result


  # merge all groups' result into a signle data.frame
  gse.res.single.table <- do.call(rbind, gse.result)
  result$gseResultSingleTable = gse.res.single.table


  # prepare heatmap
  gse.res.single.table$Adjusted.Pvalue <- -log10(gse.res.single.table$Adjusted.Pvalue)
  gse.res.single.table[ which(is.infinite(gse.res.single.table$Adjusted.Pvalue) ),  ]$Adjusted.Pvalue <- cutoff.log10  # Minimum 5
  gse.res.single.table$Adjusted.Pvalue[ gse.res.single.table$Adjusted.Pvalue >= cutoff.log10 ] = cutoff.log10

  # add direction information
  gse.res.single.table$Adjusted.Pvalue <- gse.res.single.table$Adjusted.Pvalue * sign(gse.res.single.table$Observed.score)

  # subset data frame
  gse.res.single.table <- gse.res.single.table[,c("Gene.Set.Term", "GSType", "Adjusted.Pvalue", "GroupName")]

  # cast data frame
  gse.res.single.table <- reshape::cast(gse.res.single.table, Gene.Set.Term + GSType ~ GroupName, value = "Adjusted.Pvalue", fun.aggregate = mean)
  row.names(gse.res.single.table) <- gse.res.single.table$Gene.Set.Term


  result$heatmap.df.row.type = gse.res.single.table$GSType
  names(result$heatmap.df.row.type) = gse.res.single.table$Gene.Set.Term
  result$heatmap.df.includeGSType =gse.res.single.table

  gse.res.single.table <- as.data.frame(gse.res.single.table[,-c(1,2)])
  result$heatmap.df = gse.res.single.table


  result

}



#' Read in gene set information from .gmt files
#'
#' @param gmt.path Please make sure the second column is description. The name of a file to read data values from. Should be a tab-separated text file, with one row per gene set. Column 1 has gene set names (identifiers), column 2 has gene set descriptions, remaining columns are gene ids for genes in that geneset
#'
#' @return A list, where each index represents a separate gene set.
#' @export
#'
#' @examples
load.gmt <- function(gmt.path){
  # https://www.rdocumentation.org/packages/qusage/versions/2.4.0/topics/read.gmt
  qusage::read.gmt(gmt.path)
  # https://www.rdocumentation.org/packages/GSA/versions/1.03.1/topics/GSA.read.gmt
  # GSA::GSA.read.gmt(gmt.path)
}



#' ssGSEA
#'
#' @param expr column is sample
#' @param geneset gmt path or a vector
#'
#' @return
#' @export
#'
#' @examples
#'
ssGSEA <- function(expr, geneset){

  if(length(geneset) == 1 & is.character(geneset)){
    gene.set <- loonR::load.gmt(geneset)
  }else if(is.vector(geneset)){
    gene.set = list(GeneSet = geneset)
  }else if(is.list(geneset)){
    geneset = geneset
  }else{
    stop("Gene set not right")
  }

  library(GSVA)
  library(GSEABase)

  res.ssgsea <- gsva(
    as.matrix(expr),
    gene.set,
    method = "ssgsea",
    # By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
    # When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".
    kcdf = "Gaussian",
    min.sz = 10)

  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  #定义ssGSEA_Score矫正函数
  res.ssgsea = normalize(res.ssgsea)#对ssGSEA_Score进行矫正
  res.ssgsea

}


#' GO annotation with word cloud
#'
#' @param go_ids A vector of GO IDs
#' @param measure Semantic measure for the GO similarity, see https://rdrr.io/pkg/GOSemSim/man/termSim.html
#'
#' @param cluster_method https://jokergoo.github.io/simplifyEnrichment/reference/cluster_terms.html
#'
#' @return
#' @export
#'
#' @examples
simplifyEnrichment <- function(go_ids, measure = "Rel", cluster_method = "binary_cut"){
  # 可以是其他的id而不仅仅是GO id，参考https://jokergoo.github.io/simplifyEnrichment/articles/simplifyEnrichment.html
  # Semantic measures can be used for the similarity of GO terms.
  # However, there are still a lot of ontologies (e.g. MsigDB gene sets)
  # that are only represented as a list of genes where the similarity between gene sets
  # are mainly measured by gene overlap. simplifyEnrichment provides the term_similarity()
  # and other related functions (term_similarity_from_enrichResult(), term_similarity_from_KEGG(),
  # term_similarity_from_Reactome(), term_similarity_from_MSigDB()
  # and term_similarity_from_gmt()) which calculate the similarity of terms by the gene overlapping,
  # with methods of Jaccard coefficient, Dice coefficient, overlap coefficient and kappa coefficient.

  if(!require(simplifyEnrichment)){
    devtools::install_github("jokergoo/simplifyEnrichment")
    require(simplifyEnrichment)
  }
  mat = GO_similarity(go_ids, measure = measure)

  df = simplifyGO(mat, cluster_method = cluster_method)

  compare_clustering_methods(mat)
  compare_clustering_methods(mat, plot_type = "heatmap")


}





