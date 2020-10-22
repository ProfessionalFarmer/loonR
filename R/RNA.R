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

  tempOutput = limma::topTable(fit,  n=Inf, adjust.method="BH",) # coef
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
  res <- results(dds, pAdjustMethod = "BH")
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
#' @param label Point names which you want to show in plot. If you don't want to show, set NA
#'
#' @return
#' @export
#'
#' @examples loonR::volcano_plot( tissue.exp.df.res$logFC, tissue.exp.df.res$adj.P.Val, lg2fc = 0.5, p = 0.05, label = label, restrict.vector = (tissue.exp.df.res$AUC > 0.7 & tissue.exp.df.res$AveExpr > 10)  )
volcano_plot <- function(x,y, xlab="Log2 Fold Change", ylab="-log10(Adjusted P)",
                         lg2fc = 1, p = 0.05, restrict.vector=NA, label = NA){
  # add text
  # https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

  df = data.frame(log2=x, P=y, label = label, stringsAsFactors = FALSE)
  df$Significant <- "No"
  df$Significant[ df$log2 > lg2fc  & df$P < p] <- "Up"
  df$Significant[ df$log2 < -lg2fc & df$P < p] <- "Down"

  if(!is.na(restrict.vector)){
    df$Significant[!restrict.vector] <- 'No'
  }

  palette = c()
  if(sum(df$Significant=="Down")!=0){
    palette = c("blue")
  }
  if(sum(df$Significant=="No")!=0){
    palette = c(palette, "gray")
  }
  if(sum(df$Significant=="Up")!=0){
    palette = c(palette, "red")
  }


  t = unlist( table(df$Significant) )
  df$Significant[ df$Significant == "Up" ] = paste("Up ","(",t["Up"],")",sep="")
  df$Significant[ df$Significant == "Down" ] = paste("Down ","(",t["Down"],")",sep="")
  df$Significant[ df$Significant == "No" ] = paste("No ","(",t["No"],")",sep="")

  df$Significant <- factor(df$Significant)


  df$P <- -log10(df$P)

  ggscatter(df,
            x="log2", y="P",
            xlab = xlab, ylab = ylab,
            color = "Significant", palette = palette,
            legend = "right", label = "label", font.label = c(12, "plain", "black"), repel = TRUE
  ) +
    geom_vline(xintercept=c(-lg2fc, lg2fc), col="gray") +
    # rremove("legend") +
    geom_hline(yintercept=-log10(p), col="gray")

}



#' M-A plot
#'
#' @param M Fold change or log ratio
#' @param A Average exprssion
#' @param p Adjusted P
#' @param m.cutoff Default 0
#' @param a.cutoff Default 0
#' @param p.cutoff Default 0.05
#' @param mlab Default 'Log2 Fold Change'
#' @param alab Default 'Average Expression'
#' @param restrict.vector  A TURE/FALSE factor. Only show TRUE point in the vector.
#' @param label Point names which you want to show in plot. If you don't want to show, set NA
#'
#' @return
#' @export
#'
#' @examples MA_plot(tissue.exp.df.res$logFC, tissue.exp.df.res$AveExpr, tissue.exp.df.res$adj.P.Val)
MA_plot <- function(M, A, p, m.cutoff=0, a.cutoff=0, p.cutoff=0.05,
                    mlab="Log2 Fold Change", alab="Average Expression", restrict.vector=NA, label = NA){

  df = data.frame(M=M, A=A, P = p, label = label, stringsAsFactors = FALSE)

  df$Significant <- "No"
  df$Significant[ df$M > m.cutoff  & df$P < p.cutoff & df$A > a.cutoff ] <- "Up"
  df$Significant[ df$M < m.cutoff  & df$P < p.cutoff & df$A > a.cutoff ] <- "Down"

  if(!is.na(restrict.vector)){
    df$Significant[!restrict.vector] <- 'No'
  }

  t = unlist( table(df$Significant) )
  df$Significant[ df$Significant == "Up" ] = paste("Up ","(",t["Up"],")",sep="")
  df$Significant[ df$Significant == "Down" ] = paste("Down ","(",t["Down"],")",sep="")
  df$Significant[ df$Significant == "No" ] = paste("No ","(",t["No"],")",sep="")

  df$Significant <- factor(df$Significant)

  df$P <- -log10(df$P)

  ggscatter(df,
            x="A", y="M",
            xlab = alab, ylab = mlab,
            color = "Significant", palette = c("blue", "gray",  "red"),
            legend = "right", label = "label", font.label = c(12, "plain", "black"), repel = TRUE
  )



}








