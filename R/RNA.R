#' Differential analysis by LIMMA
#'
#' @param df raw count or log2(TPM+1), default log2(TPM+1)
#' @param group factor, first control then experiment
#' @param rawcount true or false
#' @param voom true or false. If library size changed too much.
#' @param pre.filter Cutoff for mean log2(TPM)
#'
#' @return
#' @export
#'
#' @examples loonR::limma_differential(tpm.table, group)
limma_differential <- function(df, group, rawcount = FALSE, voom = FALSE, pre.filter = 0){

  library(limma)
  library(edgeR)

  if (rawcount){
    dge <- DGEList(counts = df)
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=1) # log transformation
    df <- logCPM
  }

  # https://lashlock.github.io/compbio/R_presentation.html
  if(pre.filter!=0){
    # Pre-filtering the dataset
    # http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization
    keep <- rowSums(df) > pre.filter
    df <- df[keep,]
  }

  group <- factor( group, levels = unique(as.character(group)), labels = c("Control","Experiment") )

  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(df)


  if(voom){
    df <- as.data.frame( voom(df, design, plot=FALSE)  )
  }

  AUC <- apply(df, 1,function(x){
    suppressMessages(roc <- pROC::roc(group, x)  )
    roc$auc
  })

  fit <- limma::lmFit(df, design)  # limma 因为TPM不需要normalize，所以不用voom函数。v应该是log2转换之后的
  contrast.matrix <- limma::makeContrasts(Experiment-Control,levels = design)  # Low high的顺序决定谁比谁
  fit <- limma::contrasts.fit(fit, contrast.matrix)
  fit <- limma::eBayes(fit, trend=TRUE)

  tempOutput = limma::topTable(fit,  n=Inf) # coef
  DEG_voom = na.omit(tempOutput)
  # 关联基因
  DEG_voom$REF = row.names(DEG_voom)

  DEG_voom$AUC = AUC[row.names(DEG_voom)]

  DEG_voom

}


#' Differential analysis by DESEQ2
#'
#' @param rawcount Only support raw count. Please note
#' @param group factor, first control then experiment
#' @param return.normalized.df
#' @param pre.filter if filter low expression gene
#'
#' @return
#' @export
#'
#' @examples
DESeq2_differential <- function(rawcount, group, pre.filter = 0, return.normalized.df = FALSE){

   group = factor( group, levels = unique(as.character(group)), labels = c("Control","Experiment") )

   # https://lashlock.github.io/compbio/R_presentation.html
   if(pre.filter!=0){
     # Pre-filtering the dataset
     # http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization
     keep <- rowSums(rawcount) > pre.filter
     rawcount <- rawcount[keep,]
   }

   library(DESeq2)
   # some values in assay are not integers
   countData =  as.data.frame(sapply(rawcount, as.integer),row.names = row.names(rawcount) )
   # Must add first column
   countData = data.frame(Gene=row.names(countData),countData)


   # Construct DESEQDataSet Object
   dds <- DESeqDataSetFromMatrix(
     countData=countData,
     colData=data.frame(row.names = colnames(rawcount),
                        Sample = colnames(rawcount),
                        Group = group),
     design=~Group, tidy = TRUE)

   # estimation of size factors（estimateSizeFactors) --> estimation of dispersion（estimateDispersons) --> Negative Binomial GLM fitting and Wald statistics（nbinomWaldTest）
   dds <- DESeq(dds)

   # normalized count
   normalized.count <- counts(dds, normalized=TRUE)

   if(return.normalized.df){
     return(normalized.count)
   }

   # result
   res <- results(dds)
   # sort by p-value
   res <- as.data.frame(res[order(res$padj),])

   AUC <- apply(normalized.count, 1,function(x){
     suppressMessages(roc <- pROC::roc(group, x)  )
     roc$auc
   })

   res$AUC = AUC[row.names(res)]



   res
}



#' Title Remove redundant gene expression data, and select the maximum one
#'
#' @param expression.df Please note the first column must be gene names
#'
#' @return A clean expression data.frame
#' @export
#'
#' @examples loonR::unique_gene_expression(normlized.exp.df)
unique_gene_expression <- function(expression.df){

  expression.df <- aggregate(expression.df[,-c(1)],
                             by = list(gene = expression.df[,c(1)]),
                             FUN = max,
                             na.rm = TRUE)

  row.names(expression.df) <- expression.df[,c(1)]
  expression.df <- expression.df[, -c(1)]
  expression.df

}



#' Volcano Plot
#'
#' @param x Log2 Fold Change
#' @param y Adjusted P value
#' @param xlab Default "Log2"
#' @param ylab Default "-log10(Adjusted P)"
#' @param lg2fc Default cufoff 1
#' @param p Default cutoff 0.05
#' @param restrict.vector  A TURE/FALSE factor. Only show TRUE point in the vector.
#'
#' @return
#' @export
#'
#' @examples
volcano_plot <- function(x,y, xlab="Log2 Fold Change", ylab="-log10(Adjusted P)",
                         lg2fc = 1, p = 0.05, restrict.vector=NA, label = NA){
  # add text
  # https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

  df = data.frame(log2=x, P=y, label = label, stringsAsFactors = FALSE)
  df$`Differential Expression` <- "No"
  df$`Differential Expression`[ df$log2 > lg2fc  & df$P < p] <- "Up"
  df$`Differential Expression`[ df$log2 < -lg2fc & df$P < p] <- "Down"

  if(!is.na(restrict.vector)){
    df$`Differential Expression`[!restrict.vector] <- 'No'
  }

  t = unlist( table(df$`Differential Expression`) )
  df$`Differential Expression`[ df$`Differential Expression` == "Up" ] = paste("Up ","(",t["Up"],")",sep="")
  df$`Differential Expression`[ df$`Differential Expression` == "Down" ] = paste("Down ","(",t["Down"],")",sep="")
  df$`Differential Expression`[ df$`Differential Expression` == "No" ] = paste("No ","(",t["No"],")",sep="")

  df$`Differential Expression` <- factor(df$`Differential Expression`)


  df$P <- -log10(df$P)

  ggscatter(df,
            x="log2", y="P",
            xlab = xlab, ylab = ylab,
            color = "Differential Expression", palette = c("blue", "gray",  "red"),
            legend = "right", label = "label", repel = TRUE
  ) +
  geom_vline(xintercept=c(-lg2fc, lg2fc), col="gray") +
  # rremove("legend") +
  geom_hline(yintercept=-log10(p), col="gray")

}





