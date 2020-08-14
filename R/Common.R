

#' Export figure to PPT
#'
#' @param obj ggplot object
#' @param file default is "~/test.pptx"
#' @param append
#'
#' @return NA
#' @export
#'
#' @examples export2ppt(ggplot2.obj)
#'
export2ppt <- function(obj,file="~/test.pptx", append=TRUE){
  library(export)
  graph2ppt(obj, file=file, append=append)
}



#' Title
#'
#' @param df
#' @param group
#' @param palette Default npg. The color palette to be used for coloring or filling by groups.
#' Allowed values include "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...;
#' or custom color palette e.g. c("blue", "red");
#' and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#'
#' @return ggplot2 object
#' @export
#'
#' @examples PlotPCA(expression, group)
PlotPCA <- function(df, group, palette = 'npg'){


  df_pca <- prcomp(df) #计算主成分
  df_pcs <-data.frame(df_pca$x,
                      class = factor(group)
  ) #定义分组

  #定义百分比
  percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

  library(ggplot2)
  library(ggpubr)

  p <- ggscatter(df_pcs, x="PC1", y="PC2", color="class", palette = palette) +
          xlab(percentage[1]) +
          ylab(percentage[2])
  p


}


