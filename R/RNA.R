#' Differential analysis by LIMMA
#'
#' @param df raw count or log2(TPM+1), default log2(TPM+1). Must be data.frame, not matrix.
#' @param group factor, first control then experiment
#' @param rawcount true or false
#' @param voom true or false. If library size changed too much.
#' @param pre.filter Cutoff for mean log2(TPM)
#' @param cal.AUC If to calculate AUC
#'
#' @return
#' @export
#'
#' @examples loonR::limma_differential(tpm.table, group)
limma_differential <- function(df, group, rawcount = FALSE, voom = FALSE, pre.filter = 0, cal.AUC = TRUE){

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

  if(cal.AUC){
      AUC <- apply(df, 1,function(x){
        suppressMessages(roc <- pROC::roc(group, x)  )
        roc$auc
      })
  }

  fit <- limma::lmFit(df, design)  # limma 因为TPM不需要normalize，所以不用voom函数。v应该是log2转换之后的
  contrast.matrix <- limma::makeContrasts(Experiment-Control,levels = design)  # Low high的顺序决定谁比谁
  fit <- limma::contrasts.fit(fit, contrast.matrix)
  fit <- limma::eBayes(fit, trend=TRUE)

  tempOutput = limma::topTable(fit,  n=Inf, adjust.method="BH",) # coef
  DEG_voom = na.omit(tempOutput)
  # 关联基因
  DEG_voom$REF = row.names(DEG_voom)

  if(cal.AUC){ DEG_voom$AUC = AUC[row.names(DEG_voom)] }

  DEG_voom

}



#' Calculate p value by t.test manually
#'
#' @param df
#' @param group
#' @param cal.AUC
#' @param exclude.zore
#' @param alternative c("two.sided", "less", "greater")
#'
#' @return
#' @export
#'
#' @examples
#'
#'
ttest_differential <- function(df, group, cal.AUC = TRUE, exclude.zore = FALSE, exclude.na = TRUE, alternative = "two.sided", cores = 40){

  cat("Pls note: Second unique variable is defined as experiment group\n")

  library(foreach)
  library(parallel)
  library(doParallel)

  g1 = unique(group)[1]
  g2 = unique(group)[2]

  registerDoParallel(cores=cores)
  res <- foreach::foreach(ind = rownames(df), .combine = rbind) %dopar% {

    exp.val = df[ind, ]

    f.group = group
    f.exclude = rep(FALSE, length(exp.val))

    # remove zore
    if (exclude.zore){
      f.exclude = exp.val ==0
    }
    # remove na
    if (exclude.na){
      f.exclude = f.exclude | is.na(exp.val)
    }

    # remove zore
    f.group = f.group[!f.exclude]
    exp.val = exp.val[!f.exclude]

    f.g1.no = sum(f.group==g1)
    f.g2.no = sum(f.group==g2)

    if(loonR::AllEqual(exp.val) | loonR::AllEqual(exp.val[f.group==g2]) | loonR::AllEqual(exp.val[f.group==g1]) ){
      p.val=NA
      t.statistic = NA
    } else{
      t.res <- t.test(exp.val[f.group==g2], exp.val[f.group==g1], alternative = alternative)
      p.val = t.res$p.value
      t.statistic = round(t.res$statistic, 3)
    }



    mean.diff = round( mean(exp.val[f.group==g2]) - mean(exp.val[f.group==g1]) , 3)
    mean = mean(exp.val)


    if(cal.AUC){
      if(is.na(p.val)){
        auc = NA
      }else{
        auc = round(loonR::get.AUC(exp.val, f.group),3)
      }
      f.res = c(ind, p.val, t.statistic, f.g1.no, f.g2.no, mean.diff, mean, auc)
      names(f.res) = c("Name", "P", "t statistic", g1, g2, "Difference", "Average", "AUC")
    }else{
      f.res = c(ind, p.val, t.statistic, f.g1.no, f.g2.no, mean.diff, mean)
      names(f.res) = c("Name", "P", "t statistic", g1, g2, "Difference", "Average")
    }

    f.res

  }
  res = data.frame(res, stringsAsFactors = F, check.names = F)
  res$`BH-Adjusted P` <- p.adjust(res$P, method = "BH")
  res$Difference <- as.numeric(res$Difference)

  row.names(res) <- res$Name
  res

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
DESeq2_differential <- function(rawcount, group, pre.filter = 0, return.normalized.df = FALSE, cal.AUC = TRUE){

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

  if(cal.AUC){
    AUC <- apply(normalized.count, 1,function(x){
      suppressMessages(roc <- pROC::roc(group, x)  )
      roc$auc
    })

    res$AUC = AUC[row.names(res)]
  }
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
unique_gene_expression <- function(expression.df, f = "max"){

  expression.df <- aggregate(expression.df[,-c(1)],
                             by = list(gene = expression.df[,c(1)]),
                             FUN = eval(f),
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
#' @param title
#'
#' @return
#' @export
#'
#' @examples loonR::volcano_plot( tissue.exp.df.res$logFC, tissue.exp.df.res$adj.P.Val, lg2fc = 0.5, p = 0.05, label = label, restrict.vector = (tissue.exp.df.res$AUC > 0.7 & tissue.exp.df.res$AveExpr > 10)  )
volcano_plot <- function(x, y, xlab="Log2 Fold Change", ylab="-log10(Adjusted P)",
                         lg2fc = 1, p = 0.05, restrict.vector=NA, label = NA, title = ''){
  # add text
  # https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

  library(ggpubr)

  df = data.frame(log2=as.numeric(x), P=as.numeric(y), label = label, stringsAsFactors = FALSE)
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
            xlab = xlab, ylab = ylab, title = title,
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
  df$Significant[ df$M < -m.cutoff  & df$P < p.cutoff & df$A > a.cutoff ] <- "Down"

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


#' Read Salmon output
#'
#' @param dirPath Simple set the directory which contains Salmon output folder
#'
#' @return A tximport oject
#' @export
#'
#' @examples
#' Directory tree: dir/sample_quant/quant.sf
#'
load.salmon.matrix <- function(dirPath){

  library(tximport)
  sample.salmon.pathes <- list.files(path = dirpath, full.names = TRUE)
  sample.salmon.pathes <- list.files(path = dirpath, full.names = TRUE, pattern = "isoforms.results")

  sample.names <- basename(sample.salmon.pathes)
  sample.names <- unlist(lapply( strsplit(base.name,"_"), function(x) {x[1]} ))

  cat("No. of samples:",length(sample.names),"\n")

  names(sample.salmon.pathes) <- sample.names

  tpm <- tximport(sample.salmon.pathes, type = "salmon", txIn = TRUE, txOut = TRUE)

  tpm

}



#' Title
#'
#' @param dirpath Simple set the directory which contains RSEM output folder
#' @param isoform specify data type. isoforma specific or not.
#' @param subdirs If existing in subdirectory
#'
#' @return
#' @export
#'
#' @examples
#' Directory tree subdirs = TRUE: dir/sample/prefix.genes.results or
#' Directory tree subdirs = TRUE: dir/sample/prefix.isoforms.results
#' Directory tree subdirs = FALSE: dir/sample.prefix.genes.results or
#' Directory tree subdirs = FALSE: dir/sample.prefix.isoforms.results
#'
load.rsem.matrix <- function(dirpath, isoform = FALSE, subdirs = TRUE){

  library(tximport)

  if (subdirs){
    subdirs <- list.files(path = dirpath, full.names = TRUE)
    if (isoform){
      sample.rsem.pathes <- list.files(path = subdirs, full.names = TRUE, pattern = "isoforms.results")
    } else{
      sample.rsem.pathes <- list.files(path = subdirs, full.names = TRUE, pattern = "genes.results")
    }

  }else{
    if (isoform){
      sample.rsem.pathes <- list.files(path = dirpath, full.names = TRUE, pattern = "isoforms.results")
    } else{
      sample.rsem.pathes <- list.files(path = dirpath, full.names = TRUE, pattern = "genes.results")
    }
  }
  base.name <- basename(sample.rsem.pathes)
  sample.names <- unlist(lapply( strsplit(base.name,"\\."), function(x) {x[1]} ))


  cat("No. of samples:",length(sample.names),"\n")

  names(sample.rsem.pathes) <- sample.names

  expr <- tximport(sample.rsem.pathes, type = "rsem", txIn = isoform, txOut = isoform)

  expr



}




#' Draw expression by dotplot
#' Need to be update
#'
#' @param df Row is gene, col is sample
#' @param group
#' @param nrow # of row to plot expression
#' @param stat.method wilcox.test or t.test
#' @param show.stat.p
#' @param rotate.x Default 0. numeric value specifying the rotation angle. 90 for vertical x-axis text.
#'
#' @return
#' @export
#'
#' @examples loonR::draw.expression.dotplot(candidate.combined.df, group )
draw.expression.dotplot <- function(df, group, nrow = 4, stat.method = "wilcox.test", show.stat.p = TRUE, rotate.x = 0){
  library(ggpubr)
  # expression in all samples
  candidate.plots.all.samples <- lapply(row.names(df), function(name){

    exp <- as.numeric( t(df[name,])  )
    tmp.df <- data.frame(Group = group, Expression = exp, stringsAsFactors = FALSE)

    p <- ggdotplot(tmp.df, y="Expression", x= "Group", add = "boxplot", title = name,
                   color = "Group", palette = loonR::get.palette.color("aaas", length(unique(group)), 0.7), xlab = "", show.legend = FALSE, legend = '',
                   short.panel.labs = FALSE, ylim = c(floor(min(tmp.df$Expression)), ceiling(max(tmp.df$Expression))+0.5) ) +
        rotate_x_text(angle = rotate.x)

    if(show.stat.p){
      p = p + stat_compare_means(
            aes(label = paste0("p = ", ..p.format..), method = eval(stat.method), comparisons = list(unique(group) ),
            label.y= (max(Expression)+1) ) )
    }

    p

  })
  cowplot::plot_grid(plotlist=candidate.plots.all.samples, nrow = nrow)

}



#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
log2dfQantileNormalization <- function(df){
  limma::normalizeBetweenArrays(df, method = "quantile")
}

