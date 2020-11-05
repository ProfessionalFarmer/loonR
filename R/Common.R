

#' Export figure to PPT
#' https://github.com/tomwenseleers/export
#'
#' @param obj ggplot object
#' @param file default is "~/test.pptx"
#' @param append Default TRUE
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




#' PCA plot
#' Please note: scale will be performed autmaticly
#'
#' @param df Row: sample, Column: gene expression
#' @param group
#' @param palette Default npg. The color palette to be used for coloring or filling by groups.
#' Allowed values include "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...;
#' or custom color palette e.g. c("blue", "red");
#' and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param ellipse Whether to add elipse.
#' @param legend.title  legend.title
#' @param main.title main title
#' @param alpha
#' @param return.percentage If TRUE, reture PCA percentage instead of PCA lot
#'
#' @return
#' @export
#'
#' @examples plotPCA(df, group, "aaas")
plotPCA <- function(df, group, palette = 'npg', ellipse = FALSE, legend.title = "Class", main.title = "", alpha=1, return.percentage = FALSE, pre.filter = 0.01){

  df <- df[, colMeans(df) > pre.filter]

  # Compute PCA
  df_pca <- prcomp(df, scale = TRUE) #计算主成分,强制scale

  # Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
  # factoextra::fviz_eig(res.pca)
  df_pcs <-data.frame(df_pca$x,
                      Class = factor(group) #定义分组
  )

  #解释方差比例
  pcvar <- apply(df_pca$x,2,var)
  pcvar <- pcvar/sum(pcvar)

  #利用标准差的结果计算,与上面结果一致
  #pcvar <- df_pca$sdev^2/sum(df_pca$sdev^2)

  pcvar <- round(pcvar*100,1)
  percentage <-paste(colnames(df_pcs)," (", paste(as.character(pcvar), "%", ")", sep=""),sep="")


  library(ggplot2)
  library(ggpubr)

  p <- ggscatter(df_pcs, x="PC1", y="PC2", color="Class", palette = get.palette.color(palette, n=length( levels(factor(group)) ), alpha=alpha), ellipse = ellipse) +
          xlab(percentage[1]) +
          ylab(percentage[2])

  p <- ggpar(p, legend = "right", legend.title = legend.title, main = main.title)

  p <- p + cowplot::theme_cowplot(font_family = "Arial")

  if(return.percentage){
    df_pcs
  }else{
    p
  }


}

#' Count each event type and draw pie chart
#'
#' @param data A data.frame or list object.
#' @param color ggsci color palette
#' @param colid If provide data.frame, column id need to be set
#' @param alpha Alpha value in plot
#' @param title Pie title
#' @param border Border color
#'
#' @return
#' @export
#'
#' @examples
#' plotPie(ioe.events.df$Type, title = "# of events")
#' or plotPie(ioe.events.df, col = 2, title = "# of events")
#'
plotPie <- function(data, color = "jco", colid = 2, alpha =1 , title = "", border="white"){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))
  if(n.color >= 10 ){
    stop("Please check, too many colors (More than 9)")
  }

  library(ggpubr)
  myPalette <- loonR::get.palette.color(palette = color, alpha = alpha, n = n.color)

  Prop <- unclass(table(data))

  lbls <- names(  unclass(table(data))  )
  lbls.bak <-lbls
  pct <- round(Prop/sum(Prop)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # add % to labels
  lbls <- paste(lbls,paste(" (",Prop,")",sep=""),sep="") # add value

  if (title == ""){
    title <- paste(" Total ", sum(Prop), sep = "" )
  }else{
    title <- paste(title, " (Total ", sum(Prop),")",sep = "" )
  }

  # draw
  # You can change the border of each area with the classical parameters:
  tmp.pie.df <- data.frame(Type=lbls,
                           Prop = as.numeric(Prop),
                           Label = lbls.bak,
                           stringsAsFactors = F)

  p <- ggpie(tmp.pie.df, "Prop", label = "Label", fill = "Type",
        color = border, palette = myPalette, title = title, legend = "right" , legend.title = "",
        font.family = "Arial")

  p
}

