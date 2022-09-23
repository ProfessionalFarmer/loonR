#' Analyze if gene is correlated with drug sensitivity
#'
#' @param drug.sen.df Row is drug, Col is sample. Pls perform log transformation
#' @param gene.df Row is gene Col is sample
#' @param compare.method wilcox.test t.test
#'
#' @return
#' @export
#'
#' @examples
celline_drug_target_auc_analysis <- function(drug.sen.df=NULL, gene.df = NULL, compare.method = "wilcox.test"){


  warning("Before input the drug sensitivity value, pls filter NA. In the following analysis we will not consider NA samples")
  if(is.null(drug.sen.df)|is.null(gene.df)){
    stop("Pls set drug sensititivity data and gene expression data. E.g. from https://discover.nci.nih.gov/rsconnect/cellminercdb/")
  }

  library(dplyr)
  library(foreach)

  library(utils)
  pb <- utils::txtProgressBar(style = 3)

  cat("Start AUC analysis\n")
  gene.names = rownames(gene.df)

  gene.auc.analysis.res = foreach(drugName = rownames(drug.sen.df), .combine = rbind) %do% {

    drug.sensitivity = unclass(unlist(drug.sen.df[drugName,]))
    drug.sensitivity = drug.sensitivity[!is.na(drug.sensitivity)]

    cellline.available.names = names(drug.sensitivity)

    # if lower than median, it's sensitive
    # based on median, split into T/F group
    drug.sensitivity = drug.sensitivity < median(drug.sensitivity)

    drug.res = foreach(gene = rownames(gene.df), .combine = c ) %do% {
      gene.expr.values = unclass(unlist(gene.df[gene, cellline.available.names]))
      loonR::get.AUC(gene.expr.values, drug.sensitivity)
    }

    setTxtProgressBar(pb, match(drugName,rownames(drug.sen.df))/nrow(drug.sen.df))

    drug.res = c(drugName, round(drug.res,3) )
    names(drug.res) = c("Drug",gene.names)
    drug.res

  }
  gene.auc.analysis.res = data.frame(gene.auc.analysis.res, check.names = F)
  rownames(gene.auc.analysis.res) = gene.auc.analysis.res$Drug
  gene.auc.analysis.res = gene.auc.analysis.res[,-c(1)]

  ##################################### a

  cat("Start significance analysis\n")


  gene.sig.analysis.res = foreach(drugName = rownames(drug.sen.df), .combine = rbind) %do% {

    drug.sensitivity = unclass(unlist(drug.sensitivity.df[drugName,]))
    drug.sensitivity = drug.sensitivity[!is.na(drug.sensitivity)]

    cellline.available.names = names(drug.sensitivity)

    drug.sensitivity = drug.sensitivity < median(drug.sensitivity)

    tmp.gene.exp = gene.df[,cellline.available.names]

    if(startsWith(compare.method,"wilco")){
      diff.res = loonR::MannWhitneyU_WilcoxonRankSumTest_differential(tmp.gene.exp, drug.sensitivity, cores = 1, cal.AUC = F)
    }else if(startsWith(compare.method,"t.te")){
      diff.res = loonR::ttest_differential(tmp.gene.exp, drug.sensitivity, cores = 1, cal.AUC = F)
    }

    setTxtProgressBar(pb, match(drugName,rownames(drug.sen.df))/nrow(drug.sen.df))

    diff.res = diff.res[gene.names,]

    drug.res = c(drugName, diff.res$P)
    names(drug.res) = c("Drug",gene.names)
    drug.res

  }
  close(pb)

  gene.sig.analysis.res = data.frame(gene.sig.analysis.res, check.names = F)
  rownames(gene.sig.analysis.res) = gene.sig.analysis.res$Drug
  gene.sig.analysis.res = gene.sig.analysis.res[rownames(gene.auc.analysis.res),-c(1)]


  # numeric
  gene.sig.analysis.res = loonR::convertDfToNumeric(gene.sig.analysis.res)
  gene.auc.analysis.res = loonR::convertDfToNumeric(gene.auc.analysis.res)


  # merge two table
  a = gene.auc.analysis.res
  a$Name = rownames(gene.auc.analysis.res)
  a= reshape2::melt(a, value.name = "AUC")


  b = gene.sig.analysis.res
  b$Name = rownames(gene.sig.analysis.res)
  b = reshape2::melt(b, value.name = "P")

  merged.table = dplyr::full_join(a,b, by=c("Name"="Name","variable"="variable"))
  rm(a, b)
  merged.table$BH.adj.P = p.adjust(merged.table$P, method = "BH")


  colnames(merged.table)[1:2] = c("Drug","Gene")
  split.byDrug = plyr::dlply(merged.table, "Drug", identity)


  data(CTRP.drug.annotation)
  merged.table = dplyr::left_join(merged.table, drug.annotation, by=c("Drug"="NAME"))

  res = list(
    AUC = gene.auc.analysis.res,
    Sig = gene.sig.analysis.res,
    Merged = merged.table,
    split.byDrug = split.byDrug
  )

 res
}


#' Analyze if a durg sensitivity is correlated with group/subtype
#'
#' @param drug.sen.df
#' @param group
#' @param compare.method Default wilcox.test
#' @param group.name Default Group
#'
#' @return
#' @export
#'
#' @examples
celline_drug_group_auc_analysis <- function(drug.sen.df=NULL, group = NULL, group.name = "Group", compare.method = "wilcox.test"){

  # 与前面的相反，前面是根据drug sensitivity 分成两组，用基因表达值来判断
  # 这俩是用drug sensitivity来判断group/subtype
  warning("Before input the drug sensitivity value, pls filter NA. In the following analysis we will not consider NA samples")
  if(is.null(drug.sen.df)|is.null(group)){
    stop("Pls set drug sensititivity data and group information (vector). E.g. from https://discover.nci.nih.gov/rsconnect/cellminercdb/")
  }

  library(dplyr)
  library(foreach)

  library(utils)
  pb <- utils::txtProgressBar(style = 3)

  cat("Start AUC analysis\n")

  drug.auc.analysis.res = foreach(drugName = rownames(drug.sen.df), .combine = rbind) %do% {

    drug.sensitivity = unclass(unlist(drug.sen.df[drugName,]))

    availabel.sample.ind = !is.na(drug.sensitivity)

    drug.sensitivity = drug.sensitivity[availabel.sample.ind]
    durg.group = group[availabel.sample.ind]

    # > for (controls > t >= cases)
    auc = loonR::get.AUC(drug.sensitivity, durg.group, raw = FALSE, direction = ">")

    setTxtProgressBar(pb, match(drugName,rownames(drug.sen.df))/nrow(drug.sen.df))

    drug.res = c(drugName, round(auc,3) )
    names(drug.res) = c("Drug", group.name)
    drug.res

  }
  drug.auc.analysis.res = data.frame(drug.auc.analysis.res, check.names = F)
  rownames(drug.auc.analysis.res) = drug.auc.analysis.res$Drug

  close(pb)

  ##################################### a

  cat("Start significance analysis\n")

  if(startsWith(compare.method,"wilco")){
    drug.sig.analysis.res = loonR::MannWhitneyU_WilcoxonRankSumTest_differential(drug.sen.df, group, cores = 1, cal.AUC = F)
  }else if(startsWith(compare.method,"t.te")){
    drug.sig.analysis.res = loonR::ttest_differential(drug.sen.df, group, cores = 1, cal.AUC = F)
  }

  # numeric
  drug.sig.analysis.res = loonR::convertDfToNumeric(drug.sig.analysis.res)
  warning("No worry. NA introduced because we try convert data.frame to numeric. We have solved.")
  drug.sig.analysis.res$Name = rownames(drug.sig.analysis.res)
  colnames(drug.sig.analysis.res)[1] = "Drug"


  drug.sig.analysis.res$P = formatC(drug.sig.analysis.res$P, digits = 2, format = 'e')
  drug.sig.analysis.res$`BH-Adjusted P` = formatC(drug.sig.analysis.res$`BH-Adjusted P`, digits = 2, format = 'e')
  drug.sig.analysis.res$Difference =  round(drug.sig.analysis.res$Difference, 3)
  drug.sig.analysis.res$Average = round(drug.sig.analysis.res$Average, 3)


  merged.table = dplyr::full_join(drug.sig.analysis.res, drug.auc.analysis.res, by=c("Drug"="Drug") )

  data(CTRP.drug.annotation)
  merged.table = dplyr::left_join(merged.table, drug.annotation.20220409, by=c("Drug"="NAME"))


  print("Pls note, differential analysis didn't consider sensitivity direction\n")
  print("Pls also note, we converted AUC to 1-AUC since lower IC50 is sensitive\n")

  res = list(
    AUC = drug.auc.analysis.res,
    Sig = drug.sig.analysis.res,
    Merged = merged.table
  )

  res
}






#' Calculate drug sensitivity based on oncoPredict
#'
#' @param train.df path or readRDS("DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
#' @param train.res readRDS("DataFiles/Training Data/CTRP2_Res.rds")
#' @param predict.df
#' @param minNumSamples
#'
#' @return
#' @export
#'
#' @examples
oncoPredict <- function(train.df=NULL, train.res = NULL, predict.df = NULL, minNumSamples = 10){

  # https://www.cancer.org/cancer/esophagus-cancer/treating/chemotherapy.html

  if(is.vector(train.df)){
    train.df = readRDS(train.df)
  }


  if(is.null(train.df)| is.null(train.res)){
    stop("Pls down load training data from https://osf.io/c6tfx/ Ref: https://github.com/maese005/oncoPredict")
  }
  warning("Pls check log2 transformation")

  library(oncoPredict)

  calcPhenotype(trainingExprData = train.df,
                trainingPtype = train.res,
                testExprData = predict.df,
                batchCorrect = 'eb',  #   "eb" for ComBat
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = minNumSamples,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData' )


  drug.predict.res = read.csv("./calcPhenotype_Output/DrugPredictions.csv")
  rownames(drug.predict.res) = drug.predict.res$X
  drug.predict.res = drug.predict.res[,-c(1)]
  drug.predict.res
}


#' Drug and gene correlation analysis
#'
#' @param drug.sen.df A data.frame. Row is drug, col is cell line
#' @param gene.expr A vector with name. Name is cell line
#'
#' @return
#' @export
#'
#' @examples
cellline_drug_gene_correlation <- function(drug.sen.df=NULL, gene.expr=NULL){

  # 看IC50或者AUC与基因表达之间的关系，负相关的为潜在的药物
  warning("Before input the drug sensitivity value, pls filter NA. In the following analysis we will not consider NA samples")
  if(is.null(drug.sen.df)|is.null(gene.expr)){
    stop("Pls set drug sensititivity data.frame or gene.expression vector, E.g. from https://discover.nci.nih.gov/rsconnect/cellminercdb/")
  }


  ### check names
  if(!identical(names(gene.expr), colnames(drug.sen.df))){
    stop("Gene expression vector should be named the same as drug sensitivity data.frame")
  }


  library(foreach)

  library(utils)
  pb <- utils::txtProgressBar(style = 3)

  cat("Start correlation analysis\n")

  drug.auc.analysis.res = foreach(drugName = rownames(drug.sen.df), .combine = rbind) %do% {

    drug.sensitivity = unlist(drug.sen.df[drugName,])

    na.ratio = sum(is.na(drug.sensitivity)) / length(drug.sensitivity)
    na.ind = is.na(unlist(drug.sen.df[drugName,]))

    # correlation analysis
    cor.res = cor.test(drug.sensitivity, gene.expr, use = "complete.obs", method = "spearman" )

    setTxtProgressBar(pb, match(drugName,rownames(drug.sen.df))/nrow(drug.sen.df))


    # AUC and anova P
    tmp.df <- data.frame(
       Gene = gene.expr[!na.ind],
       # # 注意，这里是小于，值越小,AUC越高
       Drug = drug.sensitivity[!na.ind] < median(drug.sensitivity[!na.ind])
    )



    auc = loonR::get.AUC(tmp.df$Gene, tmp.df$Drug, raw = F)
    anova.p = anova(lm(Drug~Gene, tmp.df))$"Pr(>F)"[1]

    drug.res = c(drugName,
                 round(cor.res$estimate,3), cor.res$p.value,
                 na.ratio,
                 auc, anova.p )
    names(drug.res) = c("Drug",
                        "Estimate", "P", "NA ratio",
                        "AUC", "ANOVA.P")
    drug.res

  }
  drug.auc.analysis.res = data.frame(drug.auc.analysis.res, check.names = F)

  drug.auc.analysis.res$P = as.numeric(drug.auc.analysis.res$P)
  drug.auc.analysis.res$Estimate = as.numeric(drug.auc.analysis.res$Estimate)
  drug.auc.analysis.res$`NA ratio` = as.numeric(drug.auc.analysis.res$`NA ratio`)
  drug.auc.analysis.res$AUC = as.numeric(drug.auc.analysis.res$AUC)
  drug.auc.analysis.res$ANOVA.P = as.numeric(drug.auc.analysis.res$ANOVA.P)

  rownames(drug.auc.analysis.res) = drug.auc.analysis.res$Drug

  close(pb)

  drug.auc.analysis.res$adjustedP = p.adjust(drug.auc.analysis.res$P)
  drug.auc.analysis.res = drug.auc.analysis.res[order(drug.auc.analysis.res$Estimate),]

  #  annotation 20220409
  # drug.annotation = readr::read_tsv("/path/to/data_CTRP-Broad-MIT_act.txt")
  # drug.annotation = drug.annotation[,c("NAME","MOA","CLINICAL.STATUS")]
  # save(drug.annotation, file = "~/Rpackage/loonR/data/CTRP.drug.annotation.rda")
  data(CTRP.drug.annotation)

  drug.auc.analysis.res = dplyr::left_join(drug.auc.analysis.res, drug.annotation, by=c("Drug"="NAME"))


  cor.plot = function(drug){
    library(ggpubr)

    na.ind = is.na(unlist(drug.sen.df[drug,]))
    y = unlist(drug.sen.df[drug,])[!na.ind]
    x = gene.expr[!na.ind]

    df = data.frame(
      Drug = y,
      Gene = x
    )
    ggscatter(df, x = "Gene", y = "Drug", cor.method = "spearman",
              cor.coef = T, add = "reg.line", title = drug ) +
        xlab("Gene expression")

  }

  list(res = drug.auc.analysis.res,
       drug.plot = cor.plot)

}


