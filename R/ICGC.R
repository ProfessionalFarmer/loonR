#' clean ICGC expression matrix file
#'
#' @param file.path file path
#'
#' @return
#' @export
#'
#' @examples
loadICGC.exp = function(file.path){

  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  expr <- fread(file.path,data.table = F)%>%
    dplyr::select(icgc_specimen_id,gene_id,normalized_read_count)%>%
    group_by(icgc_specimen_id,gene_id)%>%
    summarise_all(max)%>%
    pivot_wider(names_from = 'icgc_specimen_id',values_from = "normalized_read_count")%>%
    summarise_all(function(x){ifelse(is.na(x),0,x)})
  table(is.na(expr))
  #######将ENSG名转换为gene名
  library(org.Hs.eg.db)
  gene_id <- expr$gene_id
  geneIDselect <-select(org.Hs.eg.db,
                        keys=gene_id,
                        columns="SYMBOL",
                        keytype="ENSEMBL" )
  colnames(geneIDselect)[1] <- c("gene_id")
  expr2 <- merge(geneIDselect,expr,by='gene_id')[,-1]
  ######去除重复的行，重命名矩阵，然后就完成了
  expr3 <- expr2[complete.cases(expr2[,1]),]
  expr3 <- distinct(expr3,SYMBOL,.keep_all = T)
  rownames(expr3) <- expr3[,1]
  expr3 <- expr3[,-1]
  expr3

}

#' clean ICGC methylation data file
#'
#' @param file.path
#'
#' @return
#' @export
#'
#' @examples
loadICGC.methy = function(file.path){

  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  expr <- fread(file.path,data.table = F)%>%
    dplyr::select(icgc_specimen_id,probe_id,methylation_value)%>%
    group_by(icgc_specimen_id,probe_id)%>%
    pivot_wider(names_from = 'icgc_specimen_id',values_from = "methylation_value")%>%
    as.data.frame()

  rownames(expr) <- expr$probe_id
  expr <- expr[,-1]
  expr

}






