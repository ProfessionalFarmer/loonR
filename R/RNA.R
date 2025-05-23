#' Differential analysis by LIMMA
#'
#' @param df raw count or log2(TPM+1), default log2(TPM+1). Must be data.frame, not matrix.
#' @param group factor, first control then experiment
#' @param rawcount true or false
#' @param voom true or false. If library size changed too much.
#' @param pre.filter Cutoff for mean log2(TPM)
#' @param cal.AUC If to calculate AUC
#' @param prop.expressed.sample Default 0.5. Proportion of samples have a count greater than pre.filter
#'
#' @return
#' @export
#'
#' @examples loonR::limma_differential(tpm.table, group)
limma_differential <- function(df, group, rawcount = FALSE, voom = FALSE, pre.filter = 0, prop.expressed.sample = 0.5, cal.AUC = TRUE){

  library(limma)
  library(edgeR)

  # https://lashlock.github.io/compbio/R_presentation.html
  # Pre-filtering the dataset
  # http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization
  keep <- ( rowSums(df > pre.filter) / ncol(df) ) >= prop.expressed.sample
  df <- df[keep,]

  if (rawcount){
    dge <- DGEList(counts = df)
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=1) # log transformation
    df <- logCPM
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
        ifelse(roc$auc > 0.5, roc$auc, 1-roc$auc)
      })
  }

  fit <- limma::lmFit(df, design)  # limma 因为TPM不需要normalize，所以不用voom函数。v应该是log2转换之后的
  contrast.matrix <- limma::makeContrasts(Experiment-Control, levels = design)  # Low high的顺序决定谁比谁
  fit <- limma::contrasts.fit(fit, contrast.matrix)
  fit <- limma::eBayes(fit, trend=TRUE)

  tempOutput = limma::topTable(fit,  n=Inf, adjust.method="BH", coef = "Experiment - Control") # coef
  DEG_voom = na.omit(tempOutput)
  # 关联基因
  DEG_voom$REF = row.names(DEG_voom)

  if(cal.AUC){ DEG_voom$AUC = AUC[row.names(DEG_voom)] }

  DEG_voom

}




#' Calculate p value by MannWhitneyU (WilcoxonRankSum) Test manually
#'
#' @param df Row is gene, column is sample
#' @param group
#' @param cal.AUC If to calculate AUC
#' @param exclude.zore If to exclude zore
#' @param alternative c("two.sided", "less", "greater")
#' @param exclude.na If to exclude na
#' @param cores number of cores to use
#' @param paired Default FALSE
#'
#' @return
#' @export
#'
#' @examples
#' Nonparametric Tests of Group Differences
MannWhitneyU_WilcoxonRankSumTest_differential <- function(df, group, cal.AUC = TRUE, exclude.zore = FALSE, exclude.na = TRUE, alternative = "two.sided", cores = 40, paired=FALSE){

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
      t.res <- wilcox.test(exp.val[f.group==g2], exp.val[f.group==g1], alternative = alternative, paired=paired)
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
      names(f.res) = c("Name", "P", "W statistic", g1, g2, "Difference", "Average", "AUC")
    }else{
      f.res = c(ind, p.val, t.statistic, f.g1.no, f.g2.no, mean.diff, mean)
      names(f.res) = c("Name", "P", "W statistic", g1, g2, "Difference", "Average")
    }

    f.res

  }
  res = data.frame(res, stringsAsFactors = F, check.names = F)
  res$P <- as.numeric(res$P)
  res$`BH-Adjusted P` <- p.adjust(res$P, method = "BH")
  res$Difference <- as.numeric(res$Difference)

  row.names(res) <- res$Name
  res


}



#' Kruskal Wallis Test One Way Anova by Ranks or F test
#'
#' @param df Row is gene, column is sample
#' @param group
#' @param exclude.zore Default FALSE
#' @param exclude.na Default TRUE
#' @param cores Default 40
#' @param parameter.test Default FALSE, so perform Nonparametric Tests: KruskalWallisTest
#'
#' @return
#' @export
#'
#' @examples
#' https://jbhender.github.io/Stats506/F18/GP/Group3.html
#' Nonparametric Tests: KruskalWallisTest
#' parametric Tests: aov  ANOVA F test
oneway.anova <- function(df, group, exclude.zore = FALSE, exclude.na = TRUE, cores = 40, parameter.test=FALSE){

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

      test.df <- data.frame(Var=exp.val, Group=factor(f.group))

      if(parameter.test){
        # create parametric one-way ANOVA model
        t.res = aov(Var~Group, data = test.df)
        t.statistic = round(unlist(summary(t.res))[["F value1"]],3)
        p.val = unlist(summary(t.res))[["Pr(>F)1"]]

        f.res = c(ind, p.val, t.statistic )
        names(f.res) = c("Name", "P", "F value")

        #  the Tukey post hoc test
        detailed.res <- TukeyHSD(t.res)$Group
        diff <- unlist(detailed.res[,1])
        diff.pval <- unlist(detailed.res[,4])

        names(diff) <- paste("Diff", row.names(detailed.res))
        names(diff.pval) <- paste("P", row.names(detailed.res))

        f.res <- c(f.res, diff, diff.pval)

      }else{
        t.res <- kruskal.test(test.df$Var~test.df$Group)
        p.val = t.res$p.value
        t.statistic = round(t.res$statistic, 3)

        f.res = c(ind, p.val, t.statistic)
        names(f.res) = c("Name", "P", "Kruskal-Wallis chi-squared")


        if(!require(pgirmess)){install.packages("pgirmess")}
        #https://stats.stackexchange.com/questions/17342/is-there-a-nonparametric-equivalent-of-tukey-hsd
        # Multiple comparison test between treatments or treatments versus control after Kruskal-Wallis test.
        detailed.res <- pgirmess::kruskalmc(test.df$Var, test.df$Group)$dif.com
        diff <- unlist(detailed.res[,1])
        diff.pval <- unlist(detailed.res[,3])

        names(diff) <- paste("Diff", row.names(detailed.res))
        names(diff.pval) <- paste("P", row.names(detailed.res))

        f.res <- c(f.res, diff, diff.pval)

      }

      f.res

    }
    res = data.frame(res, stringsAsFactors = F, check.names = F)
    res$P <- as.numeric(res$P)
    res$`BH-Adjusted P` <- p.adjust(res$P, method = "BH")

    row.names(res) <- res$Name
    res

}



#' Calculate p value by t.test manually
#'
#' @param df Row is gene, column is sample
#' @param group
#' @param cal.AUC
#' @param exclude.zore
#' @param alternative c("two.sided", "less", "greater")
#' @param exclude.na
#' @param cores
#' @param paired Default FALSE
#'
#' @return
#' @export
#'
#' @examples
#' Parametric Tests of Group Differences
ttest_differential <- function(df, group, cal.AUC = TRUE, exclude.zore = FALSE, exclude.na = TRUE, alternative = "two.sided", cores = 40, paired=FALSE){

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
      t.res <- t.test(exp.val[f.group==g2], exp.val[f.group==g1], alternative = alternative, paired=paired)
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
        if(auc<0.5){auc = 1 - auc}
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
  res$P <- as.numeric(res$P)
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
#' @param prop.expressed.sample Default 0.5. Proportion of samples have a count greater than pre.filter
#' @param cal.AUC
#'
#' @return
#' @export
#'
#' @examples
DESeq2_differential <- function(rawcount, group, prop.expressed.sample = 0.5, pre.filter = 0, return.normalized.df = FALSE, cal.AUC = TRUE){

  group = factor( group, levels = unique(as.character(group)), labels = c("Control","Experiment") )
  rnames <- row.names(rawcount)
  rawcount <- loonR::convertDfToNumeric(rawcount)
  row.names(rawcount) <- rnames
  rm(rnames)

  # https://lashlock.github.io/compbio/R_presentation.html

  # Pre-filtering the dataset
  # http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization
  keep <- ( rowSums(rawcount > pre.filter)/ncol(rawcount) ) >= prop.expressed.sample
  rawcount <- rawcount[keep,]


  library(DESeq2)
  # some values in assay are not integers
  countData =  as.data.frame(sapply(rawcount, as.integer),row.names = row.names(rawcount) )

  ## Must add first column if tidy=TRUE
  countData = data.frame(Gene=row.names(countData),countData)

  names(countData) = stringr::str_replace(names(countData), "\\.","-")


  # Construct DESEQDataSet Object
  dds <- DESeqDataSetFromMatrix(
    countData=countData,
    colData=data.frame(
      row.names = colnames(rawcount),
      Sample = colnames(rawcount),
      Group = group),
      design = ~ Group,
      tidy = TRUE
    )

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
  res <- as.data.frame(res[order(res$pvalue),])

  if(cal.AUC){
    AUC <- apply(normalized.count, 1,function(x){
      suppressMessages(roc <- pROC::roc(group, x)  )
      ifelse(roc$auc > 0.5, roc$auc, 1-roc$auc)
    })

    res$AUC = AUC[row.names(res)]
  }
  res$REF <- row.names(res)

  res
}





#' Title Remove redundant gene expression data, and select the maximum one
#'
#' @param expression.df Please note the first column must be gene names if gene = Null
#' @param f Default "max"
#' @param gene If the first columnn is not gene name, input gene name here
#' @param method Default 1. 1 for aggragate, 2 for dplyr
#'
#' @return A clean expression data.frame
#' @export
#'
#' @examples
#' d <- read.table(text=
#' 'Name     Month  Rate1     Rate2
#' Aira       1      12        23
#' Aira       2      18        73
#' Aira       3      19        45
#' Ben        1      53        19
#' Ben        2      22        87
#' Ben        3      19        45
#' Cat        1      22        87
#' Cat        2      67        43
#' Cat        3      45        32', header=TRUE)
#' loonR::unique_gene_expression(d)
#'
unique_gene_expression <- function(expression.df, f = "max", gene = NULL, method=1){

  if(!is.null(gene)){
    expression.df = data.frame(
      Gene = gene,
      expression.df
    )
  }else{
    expression.df = data.frame(expression.df, check.names = F, stringsAsFactors = F)
    colnames(expression.df)[1] = c("Gene")
  }

  expression.df = switch(method,
         "1" = {
           # Method 1 Slow, So I used dplyr method
           # https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame
           expression.df <- aggregate(expression.df[,-c(1)],
                                      by = list(gene = expression.df[,c(1)]),
                                      FUN = eval(f),
                                      na.rm = TRUE)
           expression.df
         },
         "2" = {
           # Method 2
           # https://stackoverflow.com/questions/21644848/summarizing-multiple-columns-with-dplyr
           library(dplyr)
           expression.df <- expression.df  %>%
                            group_by(Gene) %>%
                            summarise(across(everything(), list(getFunction(f)), na.rm = TRUE ))
           expression.df = data.frame(expression.df, check.names = F, stringsAsFactors = F)
           colnames(expression.df) = stringr::str_remove_all(colnames(expression.df),"_1$")

           expression.df
         },
         stop(paste0("No handler for method ", id))
  )


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
#' @param col.pal Default: c("blue", "gray", "red"). A vector with 3 elements. For significantly down, not significant, significantly up
#' @param point.size Default 2
#' @param alpha Color alpha. Default 1
#'
#' @return
#' @export
#'
#' @examples loonR::volcano_plot( tissue.exp.df.res$logFC, tissue.exp.df.res$adj.P.Val, lg2fc = 0.5, p = 0.05, label = label, restrict.vector = (tissue.exp.df.res$AUC > 0.7 & tissue.exp.df.res$AveExpr > 10)  )
volcano_plot <- function(x, y, xlab="Log2 Fold Change", ylab="-log10(Adjusted P)",
                         lg2fc = 1, p = 0.05, restrict.vector=NA, label = NA,
                         title = '', col.pal=c("blue", "gray", "red"),
                         point.size = 2, alpha = 1 ){
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
    palette = col.pal[1]
  }
  if(sum(df$Significant=="No")!=0){
    palette = c(palette, col.pal[2])
  }
  if(sum(df$Significant=="Up")!=0){
    palette = c(palette, col.pal[3])
  }
  palette = loonR::get.palette.color(palette, alpha = alpha)

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
            legend = "right", label = "label",
            font.label = c(12, "plain", "black"),
            repel = TRUE, size = point.size
  ) +
    geom_vline(xintercept=c(-lg2fc, lg2fc), col="gray") +
    # rremove("legend") +
    geom_hline(yintercept=-log10(p), col="gray") +
    cowplot::theme_cowplot(font_family = "Arial")

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
#' @param col.pal Default: c("blue", "gray", "red"). A vector with 3 elements. For significantly down, not significant, significantly up
#'
#' @return
#' @export
#'
#' @examples MA_plot(tissue.exp.df.res$logFC, tissue.exp.df.res$AveExpr, tissue.exp.df.res$adj.P.Val)
MA_plot <- function(M, A, p, m.cutoff=0, a.cutoff=0, p.cutoff=0.05,
                    mlab="Log2 Fold Change", alab="Average Expression",
                    restrict.vector=NA, label = NA, col.pal=c("blue", "gray", "red") ){

  library(ggpubr)

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
            color = "Significant", palette = col.pal,
            legend = "right", label = "label", font.label = c(12, "plain", "black"), repel = TRUE
  ) + cowplot::theme_cowplot(font_family = "Arial")



}


#' Read Salmon output
#'
#' @param dirPath Simple set the directory which contains Salmon output folder
#' @param isoform specify data type. isoforma specific or not.
#' @param countsFromAbundance countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"),
#'
#' @return A tximport oject
#' @export
#'
#' @examples
#' Directory tree: dir/sample_quant/quant.sf
load.salmon.matrix <- function(dirPath, isoform = TRUE, countsFromAbundance = "no", dir.pattern = ""){

  warning("Pls note isoform paramter: ", isoform)


  library(tximport)
  sample.salmon.pathes <- list.files(path = dirPath, full.names = TRUE, pattern = dir.pattern )
  sample.names <- basename(sample.salmon.pathes)
  sample.names <- unlist(lapply( strsplit(sample.names, dir.pattern), function(x) {x[1]} ))

  if(isoform){
    sample.salmon.pathes <- list.files(path = sample.salmon.pathes, full.names = TRUE, pattern = "quant.sf")
  }else{
    sample.salmon.pathes <- list.files(path = sample.salmon.pathes, full.names = TRUE, pattern = "quant.genes.sf")
  }

  cat("No. of samples:",length(sample.names),"\n")

  names(sample.salmon.pathes) <- sample.names

  if(isoform){
    tpm <- tximport(sample.salmon.pathes, type = "salmon", txIn = isoform, txOut = isoform, countsFromAbundance = countsFromAbundance)
  }else{
    tpm <- tximport(sample.salmon.pathes, type = "salmon", txIn = isoform, txOut = isoform, geneIdCol=1, countsFromAbundance = countsFromAbundance)
  }


  tpm

}



#' Load quantification result from RSEM by tximport
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
load.rsem.matrix <- function(dirpath, isoform = FALSE, subdirs = TRUE){

  warning("Pls note isoform paramter: ", isoform)

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

  # add isoform percent 20210502
  if(isoform){
    library(dplyr)

    iso.list <- lapply(names(sample.rsem.pathes), function(sample){
        sample.iso.df <- read.table(sample.rsem.pathes[sample], header = T, sep = "\t")
        iso.pct <- sample.iso.df %>% dplyr::select(transcript_id, gene_id, IsoPct)
        colnames(iso.pct) <- c("transcript", "gene", sample)
        iso.pct
    })
    names(iso.list) <- names(sample.rsem.pathes)
    iso.pct.res <- iso.list %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("transcript", "gene") ), .)
    row.names(iso.pct.res) <- iso.pct.res$transcript

    # corrsponded with tximport
    iso.pct.res <- iso.pct.res[ row.names(expr$abundance), ]

    expr$iso.gene.map = iso.pct.res[,c(1,2)]

    expr$IsoPct <- iso.pct.res[,-c(1,2)]
  }

  # stor FPKM
  library(dplyr)

  if(isoform){
    FPKM.list <- lapply(names(sample.rsem.pathes), function(sample){
      sample.iso.df <- read.table(sample.rsem.pathes[sample], header = T, sep = "\t")
      iso.pct <- sample.iso.df %>% dplyr::select(transcript_id, FPKM)
      colnames(iso.pct) <- c("transcript", sample)
      iso.pct
    })
    names(FPKM.list) <- names(sample.rsem.pathes)
    FPKM.res <- FPKM.list %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("transcript") ), .)
    row.names(FPKM.res) <- FPKM.res$transcript
    # corrsponded with tximport
    FPKM.res <- FPKM.res[ row.names(expr$abundance), ]
    expr$FPKM <- FPKM.res[,-c(1)]
  }

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
#' @param ylab
#' @param ylim Default NA. User can specify.
#' @param alternative c("two.sided", "less", "greater")
#'
#' @return
#' @export
#'
#' @examples loonR::draw.expression.dotplot(candidate.combined.df, group )
draw.expression.dotplot <- function(df, group, nrow = 4, stat.method = "wilcox.test", show.stat.p = TRUE, rotate.x = 0, ylab="Expression", ylim = NULL, alternative = "two.sided" ){
  library(ggpubr)
  # expression in all samples
  candidate.plots.all.samples <- lapply(row.names(df), function(name){

    exp <- as.numeric( t(df[name,])  )
    tmp.df <- data.frame(Group = group, Expression = exp, stringsAsFactors = FALSE)

    # when is null, use max
    if(is.null(ylim)){
      ylim = c(floor(min(tmp.df$Expression)), ceiling(max(tmp.df$Expression))+2.5)
    }else{
      ylim = c(ylim[1],ylim[2]+2.5)
    }

    p <- ggdotplot(tmp.df, y="Expression", x= "Group", add = "boxplot", title = name,
                   color = "Group", palette = loonR::get.palette.color("aaas", length(unique(group)), 0.7), xlab = "", ylab = ylab, show.legend = FALSE, legend = '',
                   short.panel.labs = FALSE, ylim = ylim ) +
        rotate_x_text(angle = rotate.x)



    if(show.stat.p){
      p = p + stat_compare_means(
        method = stat.method, paired = FALSE,
        method.args = list(alternative = alternative),
        label.y= (max(ylim)-2) )
    }
    print(list(unique(group) ))
    p

  })
  cowplot::plot_grid(plotlist=candidate.plots.all.samples, nrow = nrow)

}


#' Quantile normalization for log2 data frame
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
log2dfQuantileNormalization <- function(df){
  res <- limma::normalizeBetweenArrays(df, method = "quantile")
  res <- as.data.frame(res)
  res

}


#' Filter genes by kinks of criteria
#'
#' @param df
#' @param expression.cutoff
#' @param proportion.expressed.samples
#' @param mean.expression
#'
#' @return
#' @export
#'
#' @examples
filterGenes <- function(df, expression.cutoff = 0, proportion.expressed.samples=0.7, mean.expression=0){

  keep1 = rowMeans(df) > mean.expression
  keep2 = (rowSums(df > expression.cutoff) / ncol(df) ) > proportion.expressed.samples
  df[keep1&keep2,]

}


#' If perform log2 transformation
#'
#' @param df
#'
#' @return Boolean variable: TRUE or FALSE
#' @export
#'
#' @examples
iflog2 <- function(df){
  df <- loonR::convertDfToNumeric(df)
  # if log2 transform.    VALUE	Normalized log2 data represent microRNA expression levels
  qx <- as.numeric(quantile(df, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    print("Should perform log2 transformation")
  }else{
    print("Note: should not perform log2 transformation")
  }
  LogC
}





#' Perform differential analysis (1 vs others)
#'
#' @param rna.df.log
#' @param group
#' @param prefix Default "Group"
#' @param cal.auc If calculate AUC value
#'
#' @return
#' @export
#'
#' @examples
compare_differential.analysis <- function(rna.df.log, group, prefix="Group", cal.auc=FALSE){

  function.analysis.res <- lapply(unique(group), function(x){

    print(paste("Now, ", prefix, x))

    ###### lapply start
    # differential analysis
    true.ind = which(group==x)
    false.ind = which(group!=x)
    limma.df = rna.df.log[ , c(false.ind, true.ind)]
    limma.diff <- loonR::limma_differential(limma.df, rep(c(FALSE,TRUE),c(length(false.ind), length(true.ind))), cal.AUC = cal.auc )

    limma.diff
  }
  )
  names(function.analysis.res) <- paste0(prefix,unique(group))
  result = list(diffResult = function.analysis.res)

  phenotype.res <- lapply(function.analysis.res, function(limma.diff){

    ## prepare input for analysis
    phenotype <- limma.diff$logFC
    # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html   -log10(res$logFC) * sign(res$logFC)
    #phenotype <- -log10(limma.diff$adj.P.Val) * sign(limma.diff$logFC)
    names(phenotype) <- row.names(limma.diff)
    phenotype <- sort(phenotype, decreasing = TRUE)
    phenotype
  })

  result$phenotype = phenotype.res
  result

}



#' Perform differential analysis (1 vs others) by t test
#'
#' @param rna.df.log
#' @param group
#' @param prefix Default "Group"
#' @param cal.auc If calculate AUC value
#'
#' @return
#' @export
#'
#' @examples
compare_differential.ttest.analysis <- function(rna.df.log, group, prefix="Group", cal.auc=FALSE){

  function.analysis.res <- lapply(unique(group), function(x){

    print(paste("Now, ", prefix, x))

    ###### lapply start
    # differential analysis
    true.ind = which(group==x)
    false.ind = which(group!=x)
    limma.df = rna.df.log[ , c(false.ind, true.ind)]
    limma.diff <- loonR::ttest_differential(limma.df, rep(c(FALSE,TRUE),c(length(false.ind), length(true.ind))), cal.AUC = cal.auc )

    limma.diff
  }
  )
  names(function.analysis.res) <- paste0(prefix,unique(group))
  result = list(diffResult = function.analysis.res)


  result

}



#' QuantileNormalization
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
log2dfQuantileNormalization <- function(df){
  limma::normalizeBetweenArrays(df, method = "quantile")
}


#' Download dataset from ArrayExpress
#'
#' @param accession_id
#' @param Processed.data processed.1.zip
#' @param relationship sdrf.txt
#' @param Array.design adf.txt
#'
#' @return
#' @export
#'
#' @examples
#' # https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6134/
#' accession_id = "E-MTAB-6134"
#' Processed.data = "E-MTAB-6134.processed.1.zip"
#' relationship = "E-MTAB-6134.sdrf.txt"
#' Array.design = "A-GEOD-13667.adf.txt"
#'
#' loonR::downloadArrayExpressDataset(accession_id)
#'
loadArrayExpressDataset <- function(accession_id = NULL, Processed.data=NULL, relationship = NULL, Array.design=NULL){

  if(missing(accession_id)|missing(Processed.data)|missing(relationship)|missing(Array.design)){
    stop("Pls set all the paramemters")
  }

  ############## read expression
  cat("Load processed data\n")
  tmp_dir = tempdir()
  dir.create(tmp_dir)


  oldw <- getOption("warn")
  options(warn = -1)



  if(endsWith(Processed.data,"zip")){
    # just a list of files inside master.zip
    tsvfile <- as.character(unzip(Processed.data, list = TRUE)$Name)

    unzip(Processed.data, exdir=tmp_dir)

    tsvfile = file.path(tmp_dir, tsvfile)

    # load the first file "file1.csv"
    processed.df = data.table::fread(tsvfile, data.table = F)

  }else{
    processed.df = data.table::fread(Processed.data, data.table = F)
  }
  options(warn = oldw)

  rownames(processed.df) = processed.df$V1
  processed.df = processed.df[,-c(1)]

  cat("Load Sample and data relationship data\n")
  ############## relationship
  if(endsWith(relationship,"zip")){
    # just a list of files inside master.zip
    tsvfile <- as.character(unzip(relationship, list = TRUE)$Name)

    # load the first file "file1.csv"
    unzip(relationship, exdir=tmp_dir)

    tsvfile = file.path(tmp_dir, tsvfile)

    relationship.df = data.table::fread(tsvfile, sep="\t")

  }else{
    relationship.df = data.table::fread(relationship, sep="\t")
  }
  rownames(relationship.df) = relationship.df$`Source Name`


  cat("Load Array design data\n")
  ################ read array design
  if(endsWith(Array.design,"zip")){
    # just a list of files inside master.zip
    tsvfile <- as.character(unzip(Array.design, list = TRUE)$Name)
    # load the first file "file1.csv"
    unzip(Array.design, exdir=tmp_dir)

    tsvfile = file.path(tmp_dir, tsvfile)

    probe.annotation.df = data.table::fread(tmp_tsv, sep="\t", fill= TRUE, blank.lines.skip = TRUE)
  }else{

    probe.annotation.df = data.table::fread(Array.design, sep="\t", fill= TRUE, blank.lines.skip = TRUE )
  }

  row_start_ind = which(probe.annotation.df == "Reporter Name", arr.ind = TRUE)[1,1]
  probe.annotation.df = probe.annotation.df[-c(1:(row_start_ind-1)),]

  colnames(probe.annotation.df) = as.character(unlist(unclass(probe.annotation.df[1,])))
  probe.annotation.df = probe.annotation.df[-c(1),]

  res = list(
    Accession = accession_id,
    Expression = processed.df,
    Phenotype = relationship.df,
    ProbeAnnotation = probe.annotation.df
  )
  res

}



#' Check gene expression
#'
#' @param gene
#' @param genes_expr
#' @param group_list
#' @param color
#' @param stat stat method
#' @param comparisions list(c("A","B"))
#' @param outlier.shape Default 19. To hide outlier, specify outlier.shape = NA. When jitter is added, then outliers will be automatically hidden.
#' @param add Default is jitter
#'
#' @return
#' @export
#'
#' @examples
check_diff_gene <- function (gene, genes_expr, group_list, color ="aaas", stat = NULL, comparisions = NULL, outlier.shape = 19, add = "jitter")
{

  library(ggpubr)
  if (length(gene) == 1) {
    if (!gene %in% rownames(genes_expr)) {
      stop(paste0(gene, " in not in your expression matrix"))
    }
    df = data.frame(value = as.numeric(genes_expr[gene,]),
                    group = group_list)
    p = ggpubr::ggboxplot(df, x = "group", y = "value", color = "group", group = "group",
                      palette = color, add = add,
                      shape = "group", outlier.shape = 19)

    if(!is.null(stat)){
      p = p + stat_compare_means(
        comparisons = comparisions
      )
    }
    p



  }
  else {
    cg = gene
    cg = cg[cg %in% rownames(genes_expr)]
    warning(paste0("Only ", length(cg), " in ", length(gene),
                   " genes are in your expression matrix\n"),
            "The following not found: ",
            paste0( setdiff(gene,cg),collapse = ", " ) ,
            "\n" )

    if (length(cg) < 1) {
      stop("None of the gene in your expression matrix")
    }

    df = tibble::as_tibble( t(genes_expr[cg,])  )
    df$group = group_list

    df.melt = reshape2::melt(df, id.vars= "group", na.rm = T, variable.name = "gene" )

    p = ggpubr::ggboxplot(df.melt, x = "group", y = "value", color = "group", group = "group",
                          palette = color, add = add,
                          shape = "group", outlier.shape = 19)

    p = facet(p, facet.by = "gene", nrow = 1)

    if(!is.null(stat)){
      p = p + stat_compare_means(
        comparisons = comparisions
      )
    }
    p


  }
}





#' Format ICGC dataframe
#'
#' @param filepath
#'
#' @return
#' @export
#'
#' @examples
readICGCExpFile = function(filepath){

  # https://zhuanlan.zhihu.com/p/544873240
  # https://zhuanlan.zhihu.com/p/545361642

  library(tidyverse)
  library(data.table)
  exp_seq <- read_delim(filepath, "\t", escape_double = FALSE, trim_ws = TRUE)

  clin = unique(exp_seq[,1:7])

  ########### raw count
  exp <- unique(exp_seq[,c("icgc_specimen_id","gene_id","raw_read_count")])

  colnames(exp) <- c("icgc_specimen_id","gene_id","raw_read_count")

  dat <- dcast(data.table(exp), gene_id~icgc_specimen_id,
               value.var="raw_read_count", fun.aggregate = max)


  ######## normalized
  exp <- unique(exp_seq[,c("icgc_specimen_id","gene_id","normalized_read_count")])

  colnames(exp) <- c("icgc_specimen_id","gene_id","normalized_read_count")

  normalized.dat <- dcast(data.table(exp), gene_id~icgc_specimen_id,
                          value.var="normalized_read_count", fun.aggregate = max)


  ########### convert to TPM
  ###### get gene length
  library(GenomicFeatures)
  txdb <- GenomicFeatures::makeTxDbFromEnsembl(organism = "Homo sapiens",
                                               server = "ensembldb.ensembl.org")

  exons_gene <- exonsBy(txdb, by = "gene")
  # https://www.biostars.org/p/83901/
  red.exonic <- reduce(exons_gene)
  exons_gene_lens <- vapply(width(red.exonic), sum, numeric(1))

  gene <- as.matrix(exons_gene_lens) %>% as.data.frame() %>% rownames_to_column()
  colnames(gene)<-c('gene_id','length')
  gene$gene_id <- str_split(gene$gene,"[.]",simplify = T)[,1] # 删除版本的“.”和后的数字
  gene$length <- as.numeric(gene$length)
  gene <- unique(gene) # 去重

  # tmp <- inner_join(dat,  gene, by = 'gene_id') # 匹配gene_count与gene_length


  ## colunmn to rownames
  counts = tibble::column_to_rownames(dat, var = "gene_id")
  normalized.dat = tibble::column_to_rownames(normalized.dat, var = "gene_id")

  ###### convert to tpm
  # https://support.bioconductor.org/p/91218/
  x <- counts/gene[match(rownames(counts),gene$gene_id),]$length
  x = na.omit(x)

  tpm = t(t(x)*1e6/colSums(x))



  res = list(clin = clin, rawCount = counts, normalized.data = normalized.dat, TPM = tpm)
  res


}




#' Convert count to TPM
#'
#' @param mat row is gene, col is sample
#' @param gene_length
#'
#' @return
#' @export
#'
#' @examples
count_to_tpm<- function(mat, gene_length){
  mat<- mat/gene_length
  total_per_sample<- colSums(mat)
  mat<- t(t(mat)/total_per_sample)
  return(mat*1000000)
}


#' Beautiful ridge plot
#'
#' @param term
#' @param value
#' @param logP
#' @param xlab
#' @param Legend
#' @param ylab
#'
#' @return
#' @export
#'
#' @examples
#' term = sample(LETTERS[1:10], 100, TRUE)
#' value = rnorm(100,0,1)
ridge_plot = function(term = NULL, value = NULL, logP = NULL, xlab = "Log2FoldChange", Legend = "-Log10(pvalue)", ylab = ""){

  if(is.null(term) | is.null(value)){
    stop("Pls provide term and value")
  }
  if(is.null(logP)){
    df = data.frame(Description = term, logfc = value, log10P = 1)
  }else{
    df = data.frame(Description = term, logfc = value, log10P = logP)
  }

  custom_colors <- colorRampPalette(c("#8075ad", "#f5edf0","#f5da73", "#ffb800"))(100)

  ggplot(df, aes(x = logfc, y = Description, fill = log10P)) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.5) +
    labs(x = xlab, y = ylab) +
    scale_fill_gradientn(colors = custom_colors,name = Legend) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )


}

