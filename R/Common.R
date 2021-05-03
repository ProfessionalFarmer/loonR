

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
#' @param label
#' @param plot3D
#' @param show.sample.name
#' @param point.size Default 2
#' @param pre.filter
#'
#' @return
#' @export
#'
#' @examples plotPCA(df, group, "aaas")
plotPCA <- function(df, group, palette = 'npg', ellipse = FALSE, legend.title = "Group",
                    main.title = "", alpha=1, return.percentage = FALSE,
                    pre.filter = 0.01, label = NULL, plot3D = FALSE,
                    show.sample.name = FALSE, point.size = 2){

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
      if(show.sample.name){
        label = row.names(df_pcs)
      }
      p <- ggscatter(df_pcs, x="PC1", y="PC2", color="Class",
                     palette = loonR::get.palette.color(palette, n=length( levels(factor(group)) ), alpha=alpha),
                     ellipse = ellipse, size = point.size,
                     label = label, repel = show.sample.name) +
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
#' @param dist.method Default euclidean. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param hclust.method Default ward.D2. The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
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
#' @param rotate.x Default 0. numeric value specifying the rotation angle. 90 for vertical x-axis text.
#'
#' @return
#' @export
#'
#' @examples
plotBarWithErr <- function(values, group, title = "", xlab = "X label", ylab = "Mean", color = "aaas", comparisons = '', method = "wilcox.test", label.y = NULL, rotate.x = 0){
  library(ggpubr)
  tmp.df <- data.frame(Value = values, Group = as.factor(group), stringsAsFactors = FALSE, check.names = F)
  colnames(tmp.df) <- c(ylab, xlab)

  p <- ggbarplot(tmp.df, x = xlab, y = ylab, add = "mean_se", fill = xlab, palette = color, title = title ) +
     rotate_x_text(angle = rotate.x)

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
                     dist(t(df), method = "euclidean"), return.score=FALSE  )

  if(return.score){
    return(data.frame(sil[,1:3]))
  }
  fviz_silhouette(sil, label = label) +
    labs(fill = class,  labels = unique(group) ) +
    scale_fill_manual(labels= unique(group), values=loonR::get.palette.color(color, length(unique(group)), alpha = alpha )  ) +
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
#' @param xlim
#' @param ylim
#'
#' @return ggplot2 object
#' @export
#'
#' @examples loonR::drawScatter(sample.info$before_reads, sample.info$mirdeep2_mapped)
drawScatter <- function(xvalue, yvalue, xlab = "X", ylab = "Y", group = NA, color = "jco", title = "", margin = TRUE, xlim = NA, ylim = NA, show.sample.name = FALSE){

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

  if (anyNA(group)){
    group = rep("Group", length(xvalue))

  }
  df = data.frame(x=xvalue, y=yvalue, Type = group)


  if (anyNA(xlim)){
    xlim = c(0,max(xvalue)*1.05)
  }
  if (anyNA(ylim)){
    ylim = c(0,max(yvalue)*1.05)
  }

  library(cowplot)
  library(ggplot2)

  # Main plot
  # pmain <- ggplot(df, aes(x = x, y = y, color = Type)) +
  #   geom_point() + theme_bw() + labs(xlab=xlab, ylab=ylab, title = title) +
  #   ggpubr::color_palette(color)  +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_blank(), axis.line = element_line(colour = "black"))


  if(show.sample.name){
    label = row.names(df_pcs)
  }

  # Main plot
  pmain <- ggpubr::ggscatter(df, x="x", y="y", color = "Type",
                             label = label, repel = show.sample.name,
                             palette = color, xlim=xlim, xlab = xlab, ylab = ylab)


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




#' Check if all values in a vector are the same
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#'
AllEqual <- structure(function(
  # from https://rdrr.io/rforge/greenbrown/src/R/AllEqual.R
  ##title<<
  ## Check if all values in a vector are the same
  ##description<<
  ## This function is used to check if all values in a vector are equal. It can be used for example to check if a time series contains only 0 or NA values.

  x
  ### numeric, character vector, or time series of type ts
) {
  res <- FALSE
  x <- na.omit(as.vector(x))
  if (length(unique(x)) == 1 | length(x) == 0) res <- TRUE
  return(res)
  ### The function returns TRUE if all values are equal and FALSE if it contains different values.
},ex=function(){
  # check if all values are equal in the following vectors:
  AllEqual(1:10)
  AllEqual(rep(0, 10))
  AllEqual(letters)
  AllEqual(rep(NA, 10))
})



#' Plot venn diagram
#'
#' @param l.list
#' @param alpha 0.5
#' @param palette aaas
#'
#' @return
#' @export
#'
#' @examples
#' l.list = list(`Up in early HCC` = early.up.genes,
#' `Down in early HCC` = early.down.genes,
#' `Up in advanced HCC` = advanced.up.genes,
#' `Down in advanced HCC` = advanced.down.genes)
#' plotVenn(l.list)
plotVenn <- function(l.list, alpha = 0.5, palette = "aaas"){
  library(VennDiagram)
  temp <- venn.diagram(l.list,
  fill = loonR::get.palette.color(palette = palette, n = length(l.list)),
  alpha = alpha,
  cex = 2, cat.fontfamily="arial",
  lty =2,  filename = NULL)
  grid.newpage()
  grid.draw(temp)
}


#' Upset plot
#'
#' @param lt
#' @param mode distinct, intersect, union
#'
#' @return
#' @export
#'
#' @examples
#' plotUpset(l.list)
plotUpset <- function(lt, mode = "intersect"){
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
  library(ComplexHeatmap)
  set.seed(123)
  m = ComplexHeatmap::make_comb_mat(lt, mode = mode)
  ss = set_size(m)
  cs = comb_size(m)


  ComplexHeatmap::UpSet(m,
    comb_order = order(cs,decreasing = T),
    left_annotation  = rowAnnotation(
      'Size' = anno_barplot(ss,
           axis_param = list(direction = "reverse"),
           border = FALSE,
           gp = gpar(fill = "black"),
           width = unit(2, "cm")
    )),
    right_annotation = NULL
    )

}

#' Convert data frame to numeric
#'
#' @param df
#'
#' @return A data.frame with numeric element in it
#' @export
#'
#' @examples
#' convertDfToNumeric(df)
convertDfToNumeric <- function(df){
  data.frame(sapply(df, function(x) as.numeric(as.character(x))), check.names = FALSE)
}


#' Scale a data.frame by column or row
#'
#' @param df
#' @param byRow Default FALSE
#' @param byColumn Default FALSE, by column
#' @param center Default TRUE. Mean = 0
#' @param scale  Default TRUE. 0-1 scale
#'
#' @return
#' @export
#'
#' @examples
#' Default by column
scaleDF <- function( df, byRow=FALSE, byColumn=FALSE, center = TRUE, scale = TRUE){
  if(byRow & byColumn){
    stop("Etheir by row or by column, can't both")
  }
  if(!byRow & !byColumn){
    warning("Set default by column")
    byColumn = TRUE
  }
  if(byRow){
    df = t(df)
    df = scale(df, center = center, scale = scale)
    df = t(df)
  }else if(byColumn){
    df = scale(df, center = center, scale = scale)
  }
  df = data.frame(df, check.names = FALSE)
  df
}


#' Draw Celventland dot plot
#'
#' @param name Names
#' @param value values
#' @param group group
#' @param palette Default aaas, if group not set, use #00AFBB
#' @param dot.size Default 2
#' @param ylab
#' @param legend.title Default "Group"
#' @param title Default ""
#'
#' @return
#' @export
#'
#' @examples
plotClevelandDot <- function(name, value, group=NA, palette = "aaas", dot.size = 2, ylab = NULL, legend.title = "Group", title = ""){

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
  library(ggpubr)
  dfm <- data.frame(Name = name,
                    Value = value,
                    Group = group,
                    stringsAsFactors = FALSE
                    )


  p <- ggdotchart(dfm, x = "Name", y = "Value",
             color = ifelse(anyNA(group), NA, "Group"),                                # Color by groups
             # c("#00AFBB", "#E7B800", "#FC4E07")
             palette = ifelse(anyNA(group)&palette=="aaas","#00AFBB", palette) , # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             rotate = TRUE,                                # Rotate vertically
             dot.size = dot.size,                          # Large dot size
             y.text.col = TRUE,                            # Color y text by groups
             ylab = ylab,
             ggtheme = theme_pubr()                        # ggplot2 theme
  )+  theme_cleveland()                                      # Add dashed grids
  p <- ggpar(p, legend.title = legend.title, main = title)



}







