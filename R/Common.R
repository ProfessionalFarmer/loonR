

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
#' @param pre.filter
#' @param lable Must be a vector or NULL
#'
#' @return
#' @export
#'
#' @examples plotPCA(df, group, "aaas")
plotPCA <- function(df, group, palette = 'npg', ellipse = FALSE, legend.title = "Class",
                    main.title = "", alpha=1, return.percentage = FALSE,
                    pre.filter = 0.01, label = NULL, plot3D = FALSE){

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

  if(plot3D){
    library(plotly)
    p = plot_ly(df_pcs,
            x = ~PC1, y = ~PC2, z = ~PC3, color = ~Class,  # c('#BF382A', '#0C4B8E')
            colors = loonR::get.palette.color(palette, n=length( levels(factor(group)) ), alpha=alpha) )

  }else{
      library(ggplot2)
      library(ggpubr)

      p <- ggscatter(df_pcs, x="PC1", y="PC2", color="Class",
                     palette = loonR::get.palette.color(palette, n=length( levels(factor(group)) ), alpha=alpha),
                     ellipse = ellipse,
                     label = label) +
              xlab(percentage[1]) +
              ylab(percentage[2])

      p <- ggpar(p, legend = "right", legend.title = legend.title, main = main.title)

      p <- p + cowplot::theme_cowplot(font_family = "Arial")
  }
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
#' @param label Whether to show labels
#'
#' @return
#' @export
#'
#' @examples
#' plotPie(ioe.events.df$Type, title = "# of events")
#' or plotPie(ioe.events.df, col = 2, title = "# of events")
#'
plotPie <- function(data, color = "jco", colid = 2, alpha =1 , title = "", border="white" , label = FALSE){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))
  if(n.color >= 9 & color != "Most" ){
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
  if(!label){
    lbls.bak = ""
  }
  tmp.pie.df <- data.frame(Type=lbls,
                           Prop = as.numeric(Prop),
                           Label = lbls.bak,
                           stringsAsFactors = F)

  p <- ggpie(tmp.pie.df, "Prop", fill = "Type", label = "Label",
        color = border, palette = myPalette, title = title, legend = "right" , legend.title = "",
        font.family = "Arial")

  p
}




#' Perform hclustering analysis
#'
#' @param df row is gene, col is sample
#' @param group
#' @param dist.method Default euclidean
#' @param hclust.method Default ward.D2
#' @param color.pla Default npg
#' @param main Title
#'
#' @return
#' @export
#'
#' @examples loonR::show_hcluster(data.frame, group)
show_hcluster <- function(df, group, dist.method = "euclidean", hclust.method = "ward.D2", color.pla = "npg", main = ""){

  library(factoextra)

  sample.dist <- dist(t(df), method = dist.method )
  sample.dist_hc <- hclust(d = sample.dist, method =hclust.method )
  p <- fviz_dend(sample.dist_hc, cex = 0.6,
            label_cols = factor(group[sample.dist_hc$order],
                                labels = loonR::get.palette.color(color.pla, length(unique(group)),0.7)
                                ),
            main = main
  )
  p

}


#' Plot bar with mean and standard error
#'
#' @param values vector
#' @param group vector
#' @param title
#' @param xlab
#' @param ylab
#' @param color
#' @param comparisons list( c("N", "T") )
#' @param method wilcox.test or t.test
#' @param label.y
#'
#' @return
#' @export
#'
#' @examples
plotBarWithErr <- function(values, group, title = "", xlab = "X label", ylab = "Mean", color = "aaas", comparisons = '', method = "wilcox.test", label.y = NULL){
  library(ggpubr)
  tmp.df <- data.frame(Value = values, Group = as.factor(group), stringsAsFactors = FALSE, check.names = F)
  colnames(tmp.df) <- c(ylab, xlab)

  p <- ggbarplot(tmp.df, x = xlab, y = ylab, add = "mean_se", fill = xlab, palette = color, title = title )
  if(comparisons != ''){# label = "p.signif"
    p <- p + stat_compare_means( method = method, comparisons =comparisons, label.y = label.y )
  }
  p
}


#' Boxplot with jitter
#'
#' @param values  vector
#' @param group  vector
#' @param title
#' @param xlab
#' @param ylab
#' @param color
#' @param comparisons list( c("N", "T") )
#' @param method wilcox.test or t.test
#' @param label.y
#'
#' @return
#' @export
#'
#' @examples
plotJitterBoxplot <- function(values, group, title = "", xlab = "X label", ylab = "Value", color = "aaas", comparisons = '', method = "wilcox.test", label.y = NULL, add = "jitter"){
  library(ggpubr)
  tmp.df <- data.frame(Value = values, Group = as.factor(group), stringsAsFactors = FALSE, check.names = F)
  colnames(tmp.df) <- c(ylab, xlab)

  p <- ggboxplot(tmp.df, x=xlab, y=ylab,
                 outlier.shape = NA, title = title,
                 add = add,  color = xlab,
                 palette = color )

  if(comparisons != ''){# label = "p.signif"
    p <- p + stat_compare_means(method = method, comparisons =comparisons, label.y = label.y )
  }
  p

}





#' Title
#'
#' @param df row is feature, column is sample
#' @param group
#' @param color
#' @param class
#' @param label Default FALSE
#'
#' @return
#' @export
#'
#' @examples
plotSilhouette <- function(df, group, color = "aaas", class = "Class", label=FALSE, alpha = 0.8){

  if (is.null(group)) {
    message("No 'group' value defined")
    stop()
  }


  library(factoextra)
  library(cluster)
  set.seed(123)

  sil <- silhouette( as.numeric(as.character(factor(group, levels = unique(group), labels = 1:length(unique(group)))   )  ),
                     dist(t(df), method = "euclidean")  )

  fviz_silhouette(sil, label = label) +
    labs(fill = class,  labels = unique(group) ) +
    scale_fill_manual( values=loonR::get.palette.color(color, length(unique(group)), alpha = alpha )  ) +
    scale_color_manual(values=loonR::get.palette.color(color, length(unique(group)), alpha = alpha )  ) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

}



#' Title
#'
#' @param row.group
#' @param col.group
#' @param row.prefix
#' @param col.prefix
#' @param lower.tail Default FALSE
#'
#' @return
#' @export
#'
#' @examples
hyperGeoTest <- function(row.group, col.group, row.prefix = "", col.prefix = "", lower.tail = FALSE){

  # https://www.omicsclass.com/article/324
  # 1-phyper(抽取样本中属于“特定类别”的数量-1,总样本中“特定类别”的数量, 总样本数-总样本中“特定类别”的数量, 从总样本中随机抽取的数量,)

    # 超几何检验，与原来的分组比较
    geomatrix  = unclass(table(row.group, col.group))
    # perform geometrix,把p值放在相同矩阵的数据框中
    tmpgeo = matrix(nrow=length(row.names(geomatrix)),ncol=length(colnames(geomatrix)))
    colnames(tmpgeo) = paste(col.prefix, colnames(geomatrix),sep="" )
    rownames(tmpgeo) = paste(row.prefix, rownames(geomatrix),sep="" )
    for(i in 1:length(row.names(tmpgeo))  ){ # row
      for(j in 1:length(colnames(tmpgeo))){  # column
        # 白球的个数，白球的总个数，黑球的总个数，抽的球（不是黑球，是球）个数
        p = phyper(geomatrix[i,j]-1, sum(geomatrix[i,]), sum(geomatrix)-sum(geomatrix[i,]), sum(geomatrix[,j]), lower.tail = lower.tail   )
        tmpgeo[i,j] = p
      }
    }
   tmpgeo = as.data.frame(tmpgeo)

   tmpgeo.log = -log10(tmpgeo)
   pheatmap::pheatmap(tmpgeo.log,
                      cluster_rows = F,
                      cluster_cols = F,
                      color = c (rep("#FFFFFF", 26), colorRampPalette(c("#FFFFFF", "#0269A4" ))(70) ),
                      breaks=unique(c(seq(0,5, length=100-1  ))),
                      display_numbers = format(tmpgeo, trim = TRUE, digits = 3, scientific = 3),
                      main = ""
   )
   # reture formated value
   format(tmpgeo, trim = TRUE, digits = 3, scientific = 3)


}




#' Draw scatter plot
#'
#' @param xvalue Vector of values
#' @param yvalue The same length wtih x value
#' @param xlab
#' @param ylab
#' @param group
#' @param color Default jco
#' @param title
#' @param margin Default TRUE, show margin plot
#'
#' @return ggplot2 object
#' @export
#'
#' @examples loonR::drawScatter(sample.info$before_reads, sample.info$mirdeep2_mapped)
drawScatter <- function(xvalue, yvalue, xlab = "X", ylab = "Y", group = NA, color = "jco", title = "", margin = TRUE, xlim = NA){

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

  if (anyNA(group)){
    group = rep("Group", length(xvalue))

  }
  df = data.frame(x=xvalue, y=yvalue, Type = group)


  if (anyNA(xlim)){
    xlim = c(0,max(xvalue)*1.05)
  }

  library(cowplot)
  library(ggplot2)

  # Main plot
  # pmain <- ggplot(df, aes(x = x, y = y, color = Type)) +
  #   geom_point() + theme_bw() + labs(xlab=xlab, ylab=ylab, title = title) +
  #   ggpubr::color_palette(color)  +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))

  # Main plot
  pmain <- ggpubr::ggscatter(df, x="x", y="y", color = "Type", palette = color, xlim=xlim, xlab = xlab, ylab = ylab)


  if (!margin){
    return(pmain)
  }

  # Marginal densities along x axis
  xdens <- axis_canvas(pmain, axis = "x") +
    geom_density(data = df, aes(x = x, fill = Type),
                 alpha = 0.7, size = 0.2) +
    ggpubr::fill_palette(color)

  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes(x = y, fill = Type),
                 alpha = 0.7, size = 0.2) +
    coord_flip()+
    ggpubr::fill_palette(color)


  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

  ggdraw(p2)


}







