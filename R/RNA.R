



#' Differential analysis by LIMMA
#'
#' @param df raw count or log2(TPM+1), default log2(TPM+1) 
#' @param group 
#' @param rawcount true or false
#' @param voom true or false. If library size changed too much.
#'
#' @return
#' @export
#'
#' @examples
limma_differential <- function(df, group, rawcount = FALSE, voom = FALSE){
  
  library(limma)
  library(edgeR)
  
  if (rawcount){
    dge <- DGEList(counts = df)
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=1)
    df <- logCPM
  }
  

  group <- factor( group, levels = unique(as.character(group)), labels = c("Control","Experiment") )
  
  design <- model.matrix(~0+group)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(df)
  
  
  if(voom){
    v <- voom(df, design, plot=FALSE)
  }
  
  fit <- limma::lmFit(df, design)
  # limma 因为TPM不需要normalize，所以不用voom函数。v应该是log2转换之后的
  # v = voom(isogmat, design = design, normalize.method = 'none')
  contrast.matrix <- limma::makeContrasts(Experiment-Control,levels = design)  # Low high的顺序决定谁比谁
  fit <- limma::contrasts.fit(fit, contrast.matrix)
  fit <- limma::eBayes(fit, trend=TRUE)
  
  tempOutput = limma::topTable(fit,  n=Inf) # coef
  DEG_voom = na.omit(tempOutput)
  # 关联基因
  DEG_voom$REF =row.names(DEG_voom)
  
  DEG_voom
}

