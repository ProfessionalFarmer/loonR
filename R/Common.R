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
export2ppt <- function(obj,file="~/test.pptx", append=TRUE){
  if(!require(export)){
    devtools::install_github("tomwenseleers/export")
  }
  graph2ppt(obj, file=file, append=append)
}





#' Overwrite write.table with modified default option. Export data frame to file
#'
#' @param df
#' @param file file path
#' @param quote Default FALSE
#' @param sep Default \\t
#'
#' @return
#' @export
#'
#' @examples
exportTable <- function(df, file="~/test.tsv", quote = F, sep = "\t", row.names = F){
  write.table(df, file = file, quote = F, sep = "\t", row.names = row.names)
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
#' @param show.total.inTitle If to show the total number in title
#'
#' @return
#' @export
#'
#' @examples
#' plotPie(ioe.events.df$Type, title = "# of events")
#' or plotPie(ioe.events.df, col = 2, title = "# of events")
plotPie <- function(data, color = "jco", colid = 2, alpha =1 , title = "", border="white" , label = FALSE, show.total.inTitle = FALSE){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))
  if(n.color >= 9 & color != "Most" & length(color) <9 ){
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
    if(show.total.inTitle){
      title <- paste(title, " (Total ", sum(Prop),")",sep = "" )
    }else{
      title <- title
    }

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
             color = border, palette = unname(myPalette), title = title,
             legend = "right" , legend.title = "",
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
#' @param cutree Number of clusters
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' loonR::show_hcluster(t(LIRI[,3:5]), LIRI$status)
show_hcluster <- function(df, group=NULL, dist.method = "euclidean", hclust.method = "ward.D2", color.pla = "npg", main = "", cutree = 0){

  # https://www.datacamp.com/community/tutorials/hierarchical-clustering-R#howto
  # https://setscholars.net/wp-content/uploads/2019/06/How-to-visualise-Hierarchical-Clustering-agglomerative-in-R.html
  library(factoextra)

  sample.dist <- dist(t(df), method = dist.method )
  sample.dist_hc <- hclust(d = sample.dist, method =hclust.method )

  if(is.null(group)){
    group = rep(1, ncol(df))
  }

  if(cutree==0){
    p <- fviz_dend(sample.dist_hc, cex = 0.6,
                   label_cols = factor(group[sample.dist_hc$order],
                                       labels = loonR::get.palette.color(color.pla, length(unique(group)),0.7)
                   ),
                   main = main
    )
    p
  }else{
    cluster.labels = stats::cutree(sample.dist_hc, k = cutree)

    p <- fviz_dend(sample.dist_hc, cex = 0.6,
                   label_cols = factor(group[sample.dist_hc$order],
                                       labels = loonR::get.palette.color(color.pla, length(unique(group)),0.7)
                   ),
                   main = main,
                   k = cutree, k_colors = loonR::get.palette.color(color.pla, cutree, 0.7),
                   color_labels_by_k = TRUE, # color labels by groups
                   rect = TRUE # Add rectangle around groups
    )
    p
    res = list(Labels=cluster.labels, Plot=p)
    res
  }



}




#' Boxplot with jitter, barplot with mean_se, violin plot with box
#'
#' @param xvalues
#' @param yvalues
#' @param group  vector
#' @param title
#' @param xlab
#' @param ylab
#' @param color
#' @param comparisons list( c("N", "T") )
#' @param method wilcox.test or t.test
#' @param label.y
#' @param add Default jitter for boxplot, mean_se for barplot, boxplot for violin and dot plot. character vector for adding another plot element (e.g.: dot plot or error bars). Allowed values are one or the combination of: "none", "dotplot", "jitter", "boxplot", "point", "mean", "mean_se", "mean_sd", "mean_ci", "mean_range", "median", "median_iqr", "median_hilow", "median_q1q3", "median_mad", "median_range"; see ?desc_statby for more details.
#' @param alternative should be one of “two.sided”, “less”, “greater”
#' @param rotate.x Default 0. numeric value specifying the rotation angle. 90 for vertical x-axis text.
#' @param group.name Default "Group"
#' @param outlier.shape point shape of outlier. Default is 19. To hide outlier, specify outlier.shape = NA. When jitter is added, then outliers will be automatically hidden.
#' @param ylim
#' @param stat Default FALSE
#' @param barplot
#' @param violin
#' @param facet stat can work only after setting facet=TRUE
#' @param dotplot
#' @param color.by.x Group by x lab. In this way we can perform stats
#' @param shape.color.by Default "black". If you want set the color by group, please set the group name (Default "Group"). Shape color by group
#' @param fill.color.by Default by "Group". If you want set the color by group, please set the group name (Default "Group"). Shape color by group
#' @param legend.pos one of c("", "top", "bottom", "left", "right", "none")
#' @param group.position Allowed values include "identity", "stack", "dodge", "position_dodge(0.9)", position_stack(). Position adjustment, either as a string, or the result of a call to a position adjustment function
#' @param remove.element Please refer https://rpkgs.datanovia.com/ggpubr/reference/rremove.html
#'
#' @return
#' @export
#'
#' @examples
#'
#' data("LIRI")
#'
#' d.frame = LIRI[,3:6]
#' group = LIRI$status
#' liri.melt <- loonR::meltDataFrameByGroup(d.frame, group)
#'
#' xvalues=liri.melt$Gene
#' yvalues=liri.melt$value
#' group=liri.melt$Group
#'
#' loonR::plotJitterBoxplot(xvalues, yvalues, group, violin = T, facet = T)
#'
plotJitterBoxplot <- function(xvalues, yvalues, group, title = "", xlab = "", ylab = "Value",
                              group.name = "Group", color = "aaas", color.by.x = FALSE,
                              comparisons = NULL, method = "wilcox.test", label.y = NULL, add = NULL,
                              alternative = "two.sided", rotate.x = 0, outlier.shape = 19, ylim=NULL, stat = FALSE,
                              barplot = FALSE, violin = FALSE, facet=FALSE, dotplot=FALSE,
                              shape.color.by = "black", fill.color.by = NULL, legend.pos = "",
                              group.position = ggplot2::position_dodge(0.9), remove.element = NULL){

  if(!is.numeric(rotate.x)){
    rotate.x = as.numeric(rotate.x)
  }

  if(color.by.x){
    group=xvalues
  }

  if(!is.factor(group)){
    group = as.character(group)
    group = factor(group)
  }

  if(is.null(fill.color.by)){
    fill.color.by = group.name
  }

  if(length(unique(group))==2){ comparisons = list( c(unique(group)) ) }

  library(ggpubr)
  tmp.df <- data.frame(X = xvalues, Y = yvalues,
                       Group = group,
                       stringsAsFactors = FALSE,
                       check.names = F)
  colnames(tmp.df)[3] <- group.name

  if(facet){

    if(barplot){
      if(is.null(add)){add="mean_se"}
      p <- ggbarplot(tmp.df, x = group.name, y="Y", add = add,
                     color = shape.color.by, fill = fill.color.by,
                     position = group.position )
    }else if(violin){
      if(is.null(add)){add="boxplot"}
      p <- ggviolin(tmp.df, x = group.name, y="Y",
                    color = shape.color.by, fill = fill.color.by,
                    shape = group.name, add = add, add.params = list(fill = "white") )
    }else if(dotplot){
      if(is.null(add)){add="boxplot"}
      p <- ggdotplot(tmp.df, x = group.name, y="Y", add = add,
                     color = shape.color.by, fill = fill.color.by,
                     position = group.position )
    }else{
      if(is.null(add)){add="jitter"}
      p <- ggboxplot(tmp.df, x=group.name, y="Y",
                     outlier.shape = outlier.shape, add = add,
                     fill = fill.color.by, color = shape.color.by)
    }
    p = ggpar(p, legend.title = "", legend = legend.pos, palette = color, ylim = ylim)
    p = p + rotate_x_text(angle = rotate.x)
    p = facet(p, facet.by = "X", nrow = 1)

    if(stat & !is.null(comparisons)  ){# label = "p.signif"
      p <- p + stat_compare_means(method = method,
                                  comparisons = comparisons,
                                  label.y = label.y,
                                  method.args = list(alternative = alternative))
    }
    p = p  + theme(strip.background = element_blank(), strip.text.x = element_text(size = 14) )

  }else{

    if(barplot){
      if(is.null(add)){add="mean_se"}
      p <- ggbarplot(tmp.df, x="X", y="Y", add = add,
                     color = shape.color.by, fill = group.name,
                     position = group.position )

    }else if(violin){
      if(is.null(add)){add="boxplot"}
      p <- ggviolin(tmp.df, x="X", y="Y",
                    fill = fill.color.by, color = shape.color.by,
                    shape = group.name,
                    add = add, add.params = list(fill = "white") )
    }else if(dotplot){
      if(is.null(add)){add="boxplot"}
      p <- ggdotplot(tmp.df, y="Y", x= "X", add = add,
                     fill = fill.color.by, color = shape.color.by,
                     short.panel.labs = FALSE, position = group.position )
    }else{
      if(is.null(add)){add="jitter"}
      p <- ggboxplot(tmp.df, x="X", y="Y",
                     outlier.shape = outlier.shape, add = add,
                     fill = fill.color.by, color = shape.color.by)
    }
    p = ggpar(p, legend.title = "", legend = legend.pos, palette = color, ylim = ylim )
    p = p + rotate_x_text(angle = rotate.x)

    if(stat & !is.null(comparisons) & color.by.x ){# label = "p.signif"
      p <- p + stat_compare_means(method = method,
                                  comparisons = comparisons,
                                  label.y = label.y,
                                  method.args = list(alternative = alternative))
    }


  }

  p = ggpar(p, xlab = xlab, ylab = ylab, title = title)
  if(!is.null(remove.element)){
    p = p + rremove(remove.element)
  }
  p


}







#' Silhouette plot
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
plotSilhouette <- function(df, group, color = "aaas", class = "Class", label=FALSE, alpha = 0.8, return.score = FALSE){

  if (is.null(group)) {
    message("No 'group' value defined")
    stop()
  }


  library(factoextra)
  library(cluster)
  set.seed(123)

  sil <- silhouette( as.numeric(as.character(factor(group, levels = unique(group), labels = 1:length(unique(group)))   )  ),
                     dist(t(df), method = "euclidean"), return.score=return.score  )

  if(return.score){
    return(data.frame(sil[,1:3]))
  }
  fviz_silhouette(sil, label = label) +
    labs(fill = class,  labels = unique(group) ) +
    scale_fill_manual(labels= unique(group), values=loonR::get.palette.color(color, length(unique(group)), alpha = alpha )  ) +
    scale_color_manual(values=loonR::get.palette.color(color, length(unique(group)), alpha = alpha )  ) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

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
#' @param show.sample.name Default FALSE
#' @param label Default NULL
#' @param cor.coef TRUE/FALSE to show coefficient
#' @param remove.legend Default FALSE
#' @param cor.method	method for computing correlation coefficient. Allowed values are one of "pearson", "kendall", or "spearman"
#' @param add c("none", "reg.line", "loess")
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' library(ggpubr)
#' data("mtcars")
#' loonR::drawScatter(mtcars$wt, mtcars$mpg, xlab = "wt", ylab = "mpg",  remove.legend = T, cor.coef = F)
drawScatter <- function(xvalue, yvalue, xlab = "X", ylab = "Y", group = NA,
                        color = "jco", title = "", remove.legend = FALSE,
                        margin = TRUE, xlim = NULL, ylim = NULL, add = "none",
                        show.sample.name = FALSE, label = NULL, cor.coef = F, cor.method = "pearson" ){

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
                             palette = color, xlim=xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                             title = title, add = add,
                             cor.coef = cor.coef,
                             cor.coeff.args = list(method = cor.method,
                                                   label.x.npc = "left",
                                                   label.y.npc = "top")
  )

  if(remove.legend){
    pmain = pmain + ggpubr::rremove("legend")
  }

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
#'
#' @examples
#' This function is no longer used and exported
#'
#' l.list = list(`Up in early HCC` = early.up.genes,
#' `Down in early HCC` = early.down.genes,
#' `Up in advanced HCC` = advanced.up.genes,
#' `Down in advanced HCC` = advanced.down.genes)
#' plotVenn(l.list)
plotVennLegacy <- function(l.list, alpha = 0.5, palette = "aaas"){
  library(VennDiagram)
  temp <- venn.diagram(l.list,
                       fill = loonR::get.palette.color(palette = palette, n = length(l.list)),
                       alpha = alpha,
                       cex = 2, cat.fontfamily="arial",
                       lty =2,  filename = NULL)
  grid.newpage()
  grid.draw(temp)
}






#' Plot venn diagram by ggvenn
#'
#' @param l.list a list as input data
#' @param alpha Default 0.5
#' @param palette Default aaas
#' @param show_elements Show set elements instead of count/percentage.
#' @param show_percentage Show percentage for each set.
#' @param text_size Default 4, Text size for intersect contents.
#' @param set_name_size Default 6, Text color for set names.
#'
#' @return The ggplot object to print or save to file.
#' @export
#'
#' @examples plotVenn(l.list)
#' l.list = list(`Up in early HCC` = early.up.genes,
#' `Down in early HCC` = early.down.genes,
#' `Up in advanced HCC` = advanced.up.genes,
#' `Down in advanced HCC` = advanced.down.genes)
#' plotVenn(l.list)
plotVenn <- function(l.list, alpha = 0.5, palette = "aaas", show_elements = FALSE, show_percentage = TRUE, text_size = 4, set_name_size = 6){

  if(!require("ggvenn")){
    devtools::install_github("yanlinlin82/ggvenn")
  }
  ggvenn(l.list, columns = names(l.list),
         show_elements = show_elements,
         show_percentage = show_percentage,
         digits = 2,
         fill_color = loonR::get.palette.color(palette = palette, n = length(l.list)),
         fill_alpha = alpha,
         stroke_color = "black",
         stroke_alpha = alpha,
         stroke_size = 0.5,
         stroke_linetype = "solid",
         set_name_color = "black",
         set_name_size = 6,
         text_color = "black",
         text_size = text_size,
         label_sep = ","
  )


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
  row.n <- row.names(df)
  col.n <- colnames(df)
  res <- data.frame(matrix(sapply(data.frame(df, check.names = F), function(x) as.numeric(as.character(x))), nrow = nrow(df)), check.names = FALSE)
  rownames(res) <- as.character(row.n)
  colnames(res) <- as.character(col.n)
  res

}


#' Scale a data.frame by column or row
#'
#' @param df
#' @param byRow Default FALSE
#' @param byColumn Default FALSE, by column
#' @param center Default TRUE. Mean = 0
#' @param scale  Default TRUE. 0-1 scale
#' @param maxUnit Default 4
#'
#' @return
#' @export
#'
#' @examples
#' Default by column
scaleDF <- function( df, byRow=FALSE, byColumn=FALSE, center = TRUE, scale = TRUE, maxUnit = 4){
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
  df[df > maxUnit ] = maxUnit
  df[df < (-1*maxUnit) ] = -1 * maxUnit
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
#' @param cleveland
#' @param lollipop
#' @param value.name xlab, Default "Value"
#' @param sample.name ylab, Default "Name"
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' loonR::plotClevelandDot(1:10, LIRI$ANLN[1:10], lollipop = T )
#'
plotClevelandDot <- function(name, value, group=NA, palette = "aaas", dot.size = 2, ylab = NULL, legend.title = "Group", title = "",
                             cleveland = TRUE, lollipop = FALSE, value.name = "Value",  sample.name = "Name"){

  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
  library(ggpubr)
  dfm <- data.frame(Name = name,
                    Value = value,
                    Group = group,
                    stringsAsFactors = FALSE
  )
  if(lollipop){
    cleveland = FALSE
  }
  if(cleveland){
    p <- ggdotchart(dfm, x = "Name", y = "Value",
                    color = ifelse(anyNA(group), "black", "Group"),                                # Color by groups
                    # c("#00AFBB", "#E7B800", "#FC4E07")
                    palette = ifelse(anyNA(group)&palette=="aaas","#00AFBB", palette) , # Custom color palette
                    sorting = "descending",                       # Sort value in descending order
                    rotate = TRUE,                                # Rotate vertically
                    dot.size = dot.size,                          # Large dot size
                    y.text.col = TRUE,                            # Color y text by groups
                    ylab = ylab, position = position_dodge(0.3),
                    ggtheme = theme_pubr()                        # ggplot2 theme
    )  +  theme_cleveland()                                      # Add dashed grids
  }else if(lollipop){
    p <- ggdotchart(dfm, x = "Name", y = "Value",
                    color = ifelse(anyNA(group), "black", "Group"),                                # Color by groups
                    # c("#00AFBB", "#E7B800", "#FC4E07")
                    palette = ifelse(anyNA(group)&palette=="aaas","#00AFBB", palette) , # Custom color palette
                    sorting = "ascending",                        # Sort value in descending order
                    add = "segments", position = position_dodge(0.3),                             # Add segments from y = 0 to dots
                    ggtheme = theme_pubr()    # ggplot2 theme
    )


  }


  p <- ggpar(p, legend.title = legend.title, main = title, ylab = value.name, xlab = sample.name)
  p


}



#' Determine the best cluster numbers
#'
#' @param df Row is sample and column is variable
#' @param distance This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL".
#' @param dissimilarity dissimilarity matrix to be used. By default, diss=NULL, but if it is replaced by a dissimilarity matrix, distance should be "NULL".
#' @param min.nc 2 minimal number of clusters
#' @param max.nc 8 maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc.
#' @param method kmeans. the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#'
#' @return
#' @export
#'
#' @examples
determineClusterNumber <- function(df, distance = "euclidean", method = "kmeans", dissimilarity = NULL, min.nc = 2, max.nc = 8){

  # reference
  # https://towardsdatascience.com/10-tips-for-choosing-the-optimal-number-of-clusters-277e93d72d92
  # https://www.cnblogs.com/think90/p/7133753.html

  if(!require("factoextra")){
    BiocManager::install("factoextra")
  }
  if(!require("NbClust")){
    BiocManager::install("NbClust")
  }

  set.seed(1234)

  if(!is.null(dissimilarity)){
    distance = NULL
  }
  result = list()

  # nbclust
  # index: the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott",
  # "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
  # "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn",
  # "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau),
  # "alllong" (all indices with Gap, Gamma, Gplus and Tau included).

  nb_clust <- NbClust::NbClust(df,
                               diss = NULL,
                               distance = distance,
                               min.nc = min.nc,
                               max.nc = max.nc,
                               method = method,
                               index = "alllong",
                               alphaBeale = 0.1)
  barplot(table(nb_clust$Best.nc[1,]),xlab = "聚类数",ylab = "支持指标数")
  result$nb_clust = nb_clust


  #### mclust
  # https://www.cnblogs.com/think90/p/7133753.html
  if(!require("mclust")){BiocManager::install("mclust")}

  m_clust <- mclust::Mclust(as.matrix(df), G=min.nc:max.nc) #聚类数目从1一直试到20
  result$m_clust = m_clust
  plot(m_clust, "BIC")


  # Weighted Sum of squares
  #fviz_nbclust(dataset, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2)
  #km.res <- kmeans(dataset,3)
  #fviz_cluster(km.res, data = dataset)
  wssplot <- function(data, nc=15, seed=1234){
    wss <- (nrow(data)-1)*sum(apply(data,2,var))
    for (i in 2:nc){
      set.seed(seed)
      wss[i] <- sum(kmeans(data, centers=i)$withinss)
    }
    plot(1:nc, wss, type="b", xlab="Number of Clusters",
         ylab="Within groups sum of squares")
  }
  wssplot(df)

  # PAM(Partitioning Around Medoids) 围绕中心点的分割算法
  if(!require("fpc")){BiocManager::install("fpc")}

  pamk.best <- pamk(df)
  pamk.best$nc

  library(cluster)
  clusplot(pam(df, pamk.best$nc))
  result$pamk = pamk.best


  # 轮廓系数Average silhouette method
  require(cluster)
  library(factoextra)
  fviz_nbclust(dataset, kmeans, method = "silhouette")


  # Gap Statistic
  library(cluster)
  set.seed(123)
  gap_clust <- clusGap(df, kmeans, 10, B = 50, verbose = interactive())

  result$gap_clust = gap_clust

  library(factoextra)
  fviz_gap_stat(gap_clust)


  #层次聚类
  h_dist <- dist(as.matrix(df))
  h_clust<-hclust(h_dist)
  plot(h_clust, hang = -1, labels = FALSE)
  rect.hclust(h_clust)



  result
}


#' Gap statistics
#'
#' @param df Row is sample, column is variable
#' @param dist a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman"
#' @param method Default "Average". The agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
#'
#' @return
#' @export
#'
#' @examples
gapStat <- function(df, dist="spearman", method="average"){

  fun<-function(df, k) {
    y<-hclust(as.dist(1-cor(df, method="spearman")), method="average");
    clus<-cutree(y, k=k);
    return(list(cluster=clus))
  }

  gaps_default <- clusGap(t(sdat), FUNcluster=fun, K.max=8, B=50)
  factoextra::fviz_gap_stat(gap_clust)

}



#' Split group by various methods. Useful
#'
#' @param group
#' @param values
#' @param quantile.cutoff 0-1. percent value for quantile function
#' @param cut.point Default NULL. Not include the minimum and maximum value
#' @param cut.label length should be the length(quantile.cutoff), or length(cut.point) plus 1
#' @param specific.group If specified, only split within these group
#' @param group.prefix ""
#' @param fun e.g. mean, median
#' @param sample.names Default NA
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(111)
#' value = runif(100, min=0, max=100)
#' set.seed(111)
#' group = sample(c(1,2,3),100, replace = TRUE)
#' table(group)
#'
#' res <- loonR::splitGroupByCutoff(group, value, fun="mean", cut.label = c("L","H"), group.prefix = "G", specific.group = c(1,2))
#' table(res$New.Label)
splitGroupByCutoff <- function(group = NULL, values = NULL, fun = NULL, quantile.cutoff = NULL, sample.names = NULL,
                               cut.point = NULL, cut.label = NULL, specific.group = NULL, group.prefix = NULL){

  if(is.null(values)){
    stop("Please input a vector including values")
  }
  if(is.null(cut.label)){
    stop("Please input labels after spliting")
  }

  if(is.null(group) ){
    group = rep("",length(values))
  }else if(length(group)==1){
    group = rep(group,length(values))
  }else if( sum(is.null(group)|is.na(group))!=0 | sum(is.null(values)|is.na(values))!=0  ){
    stop("group/values may have NA/NULL")
  }


  data.df <- data.frame(Group = group, Value = as.numeric(values), Label = "",
                        check.names = F, stringsAsFactors = F)

  if(!is.null(sample.names)){
    row.names(data.df) <- sample.names
    data.df$Name = row.names(data.df)
  }

  if(is.null(cut.point) & is.null(fun) & is.null(quantile.cutoff)){
    stop("Cut point or fun or quantile cutoff should be set and only set one")
  }

  if(length(cut.label)!=2 & !is.null(fun) ){  stop("Error, cut.label should be two elements")  }
  if(length(cut.label)!=(length(quantile.cutoff)+1)& !is.null(quantile.cutoff)){        stop("Error, according to quantile.cutoff, cut.label should be", length(quantile.cutoff)+1,"elements")      }
  if(length(cut.label)!=(length(cut.point)+1) & !is.null(cut.point)){ stop("Error, according to cut.point, cut.label should be", length(cut.point)+1,"elements")  }


  ######### if not specify group
  if(is.null(specific.group)){

    global.cut = 0

    ###### if set a function, e.g. mean median
    if(!is.null(fun)){
      global.cut = get(fun)(data.df$Value)

      ###### if set a quantile cutoff
    }else if(!is.null(quantile.cutoff)){
      global.cut = quantile(data.df$Value, prob=c(quantile.cutoff))

      ###### if set cut point
    }else if(!is.null(cut.point)){
      global.cut = cut.point

    }

    data.df$Label <- cut(data.df$Value,
                         c(min(data.df$Value)-0.1, global.cut, max(data.df$Value)),
                         labels = cut.label )

    data.df$Label <- as.character(data.df$Label)
    cat("Cutpoint for is ", global.cut)

    ######### if specify group
  }else{

    if(length(unique(intersect(data.df$Group, specific.group)))!=length(specific.group)){
      stop("Please specify the specific.group corresponded to group")
    }

    for(g in specific.group){

      g.index = data.df$Group == g
      local.cut = 0

      g.value = data.df$Value[g.index]

      ###### if set a function, e.g. mean median
      if(!is.null(fun)){
        local.cut = get(fun)(g.value)

        ###### if set a quantile cutoff
      }else if(!is.null(quantile.cutoff)){
        local.cut = quantile(g.value, prob=c(quantile.cutoff))

        ###### if set cut point
      }else if(!is.null(cut.point)){
        local.cut = cut.point

      }

      g.label <- cut(g.value,
                     c(min(g.value)-0.1, local.cut, max(g.value)),
                     labels = cut.label )

      data.df$Label[g.index] = as.character(g.label)

      cat("Cutpoint for ",g," is ", local.cut)
    }

    data.df$Label <- as.character(data.df$Label)

  }


  if(!is.null(group.prefix)){
    data.df$Group = paste(group.prefix, data.df$Group, sep="")
  }

  if(loonR::AllEqual(group)){ # if only one group
    data.df$Label[data.df$Label!=""] = paste("", data.df$Label[data.df$Label!=""], sep="")
  }else{
    data.df$Label[data.df$Label!=""] = paste("-", data.df$Label[data.df$Label!=""], sep="")
  }


  data.df$New.Label = paste(data.df$Group, data.df$Label, sep="")
  data.df$Label = stringr::str_remove_all(data.df$Label, "^-")


  data.df

}



#' Join a table list
#'
#' @param list.df List of data.frame
#' @param sfx If not specified, use list element names instead
#' @param by 	A character vector of variables to join by.
#' @param full If TRUE, use full_join
#' @param inner If TRUE, use inner_join
#' @param left  If TRUE, use left_join
#' @param right If TRUE, use right_join
#'
#' @return
#' @export
#'
#' @examples
reduce_join_list <- function(list.df, sfx=NULL, by = NULL, full = FALSE, inner = FALSE, left = FALSE, right = FALSE){


  if(full){
    f = dplyr::full_join
  }else if(inner){
    f = dplyr::inner_join
  }else if(left){
    f = dplyr::left_join
  }else if(right){
    f = dplyr::right_join
  }else{
    warning("Pls specify join type: full, left, right or inner join?
            Default will use full_join")
    f = dplyr::full_join
  }



  if(is.null(sfx)){
    sfx <- names(list.df)
  }

  if(is.null(by)){
    stop("Pls specify by option")
  }

  res <- list.df[[1]]
  for(i in head(seq_along(list.df), -1)) {

    res <- eval(f)(
      res, list.df[[i+1]],
      all = TRUE,
      suffix = c(paste0('.',sfx[i]), paste0('.',sfx[i+1])),
      by = by
    )
  }

  res

}





#' Get current directory path
#'
#' @return
#' @export
#'
#' @examples
thisPath <- function() {
  # from 生信开发者 wechat
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {    stop("Cannot find file path")
  }
}

#' Get current running program name
#'
#' @param arguments
#'
#' @return
#' @export
#'
#' @examples
getProgramName <- function(arguments){
  # from 生信开发者 wechat
  args <- commandArgs(trailingOnly = FALSE)
  sub("--file=", "", args[grep("--file=", args)])
}


#' Generic is.nan function
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' # This function can be applied to data.frame
#' # Just use is.nan(data.frame)
is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}

#' Gnenerate combinations with specified size
#'
#' @param panel panel which will be used to make combinations
#' @param size size of a single combination. Can be a signle value or a vector value
#' @param repeats Default FALSE, not allow repeats
#' @param vector
#' @param sep If return a vector, sep is needed
#'
#' @return
#' @export
#'
#' @examples
generateCombinations <- function(panel=NULL, size = 0, repeats=FALSE, vector=FALSE, sep ="-"){
  if(length(size)>1){

    res <- lapply(size, function(x){
      combs = gtools::combinations(length(panel), x, panel, repeats=repeats)
      if(vector==TRUE){
        if(x==1){
          combs = unclass(unlist(combs[,1]))
        }else{
          combs = apply(combs, 1, function(x) paste0(sort(x), sep="", collapse = sep ))
        }
        combs
      }
    })
    names(res) <- size

  }else{
    if(is.null(panel)){
      stop("Please set a panel which you want to make a combination")
    }
    if(size == 0){
      size = length(panel) - 1
    }
    if(size > length(panel)){
      stop("Pls set size less than the length of panel")
    }
    res <- gtools::combinations(length(panel), size, panel, repeats=repeats)

    if(vector==TRUE){
      res <- apply(res, 1, function(x) paste0(sort(x), sep="", collapse = "-"))
    }
  }
  data.frame(res)
}




#' Generate label based on group
#'
#' @param label Number of elements in label vector should be the same with total group number
#' @param ... one or more group vector
#'
#' @return
#' @export
#'
#' @examples
#' n.group = c("N1", "N2", "N3", "N4")
#' t.group = c("T1", "T2")
#' m.group = c("M1", "M2", "M1")
#'
#' labels <- genereateLabelsByGroup(c("N","T", "M"), n.group, t.group, m.group)
#' table(labels)
#'
genereateLabelsByGroup <- function(label=NULL,...){

  l = list(...)

  if(is.null(label)){
    stop("Please set labels which will be used as label")
  }

  if(length(l)!=length(label)){
    stop("Please input vectors in order to generate labels. Number of vectors should be the same as element in labels")
  }

  label.list = lapply(1:length(l), function(x) {
    g.labels <- rep( label[x], length(l[[x]])  )
    names(g.labels) = l[[x]]
    g.labels
  }
  )


  return( unlist(label.list) )

}



#' Melt data frame by group
#'
#' @param d.frame Row is sample, column is feature.
#' @param group
#' @param na.rm Default TRUE. Should NA values be removed from the data set?
#' @param variable_name Default "Gene"
#' @param group2 The second group
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' d.frame = LIRI[,3:6]
#' group = LIRI$status
#' head( loonR::meltDataFrameByGroup(d.frame, group) )
meltDataFrameByGroup <- function(d.frame=NULL, group=NULL, na.rm = TRUE, variable_name="Gene", group2=NULL){

  if(is.null(d.frame)|is.null(group)){
    stop("Please set data.frame and group")
  }

  if(!is.null(group2)){
    if(length(group)!=length(group)) stop("Group 1 and 2 not the same length")
    group = data.frame(
      cond1=group,
      cond2=group2
    )
  }

  if(is.vector(group)|is.factor(group)){
    melt.df <- data.frame(d.frame, Group = group,
                          check.names = F, stringsAsFactors = F)

    melted.df <- reshape2::melt(melt.df, id.vars="Group", na.rm = na.rm, variable.name = variable_name )

  }else if(is.data.frame(group)|is.matrix(group)){
     if(nrow(d.frame) != nrow(group)){
       stop("Group number is not the same with data.frame")
     }

    melt.df <- data.frame(d.frame, group,
                          check.names = F, stringsAsFactors = F)
    melted.df <- reshape2::melt(melt.df, id.vars=colnames(group), na.rm = na.rm, variable.name = variable_name )

  }
  melted.df

}



#' Radar plot (aka spider plot)
#'
#' @param df Column is variable, Row is samples/class
#' @param palette Default aaas
#' @param min min of the axis
#' @param max max of the axis
#' @param fill.color Default is FALSE, no fill color
#' @param axistype Default 0. The type of axes, specified by any of 0:5. 0 means no axis label. 1 means center axis label only. 2 means around-the-chart label only. 3 means both center and around-the-chart (peripheral) labels. 4 is *.** format of 1, 5 is *.** format of 3. Default is 0.
#' @param cglcol Default 'navy'. color of the net
#' @param vlabels Character vector for the names for variables. If NULL, the names of the variables as colnames(df) are used. Default NULL.
#'
#' @return
#' @export
#'
#' @examples
#' # https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html
#' # Create data: note in High school for several students
#' set.seed(99)
#' df <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
#' colnames(df) <- c("math" , "english" , "biology" , "music" , "R-coding" )
#' rownames(df) <- paste("mister" , letters[1:3] , sep="-")
#' loonR::radarSpiderPlot(df)
#'
radarSpiderPlot <- function(df, palette = "aaas", min = 0, max = NULL, fill.color=FALSE, axistype = 0, cglcol = "navy", vlabels = NULL){
  # https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html
  # Use ggradar
  # https://github.com/ricardo-bion/ggradar

  if(!require("fmsb")){
    BiocManager::install("fmsb")
    library("fmsb")
  }

  df = as.data.frame(df)

  if(min==0){
    warning("Pls note the minimum tick for plot is 0")
  }
  if(is.null(max)){
    max = max(df)
  }


  df <- rbind(rep(max,ncol(df)) , rep(min,ncol(df)) , df)

  colors_border = loonR::get.palette.color(palette, ncol(df) )

  if(fill.color){
    colors_in = loonR::get.palette.color(palette, ncol(df), alpha = 0.2)
  }else{
    colors_in = NA
  }


  radarchart(
    df, axistype = axistype,
    #custom polygon
    pcol = colors_border , # line color
    pfcol = colors_in , # fill color
    plwd = 4 ,
    plty = 1,
    vlabels = vlabels,
    #custom the grid
    cglcol = "navy", # color of the net
    cglty = 3, # net line type # https://www.r-graph-gallery.com/6-graph-parameters-reminder.html
    axislabcol="grey", # color of axis labels
    # caxislabels=seq(0,20,5), #vector of axis labels to display
    cglwd=0.8, # net width
    #custom labels
    vlcex=0.8 #  group labels size
  )

  # Add a legend
  legend(x=0.7, y=1.4, legend = rownames(df[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)

}


#' Convert a text table to figure
#'
#' @param data
#' @param rowname Default TRUE
#' @param ttheme Default blank. character string the table style/theme. The available themes are illustrated in the ggtexttable-theme.pdf file. Allowed values include one of c("default", "blank", "classic", "minimal", "light", "lBlack", "lBlue", "lRed", "lGreen", "lViolet", "lCyan", "lOrange", "lBlackWhite", "lBlueWhite", "lRedWhite", "lGreenWhite", "lVioletWhite", "lCyanWhite", "lOrangeWhite", "mBlack", "mBlue", "mRed", "mGreen", "mViolet", "mCyan", "mOrange", "mBlackWhite", "mBlueWhite", "mRedWhite", "mGreenWhite", "mVioletWhite", "mCyanWhite", "mOrangeWhite" ). Note that, l = "light"; m = "medium"
#' @param top.black.line Default 1:2. Black line. From top
#' @param bottom.black.line Default 1. Black line From bottom
#' @param subtitle Subtitle if you want to show
#' @param main.title Main title if you want to show
#' @param footnote Italic footnote at botton if you want to show
#' @param vline
#'
#' @return
#' @export
#'
#' @examples
#' data(iris)
#' loonR::table2Figure( head(iris), main.title = "IRIS",  footnote = "Here is footnote")
table2Figure <- function(data, rowname = TRUE, ttheme = "blank", top.black.line = 1:2, bottom.black.line =1, subtitle = NULL, main.title = NULL, footnote = NULL, vline = NULL){

  if(rowname){
    rowname = row.names(data)
  }else{
    rowname=NULL
  }

  library(ggpubr)
  library(dplyr)

  colname = colnames(data)

  bottom.black.line = nrow(data) + 2 - bottom.black.line

  ggtexttable(data, rows = rowname,
              cols = colname,
              theme = ttheme(ttheme)
              ) %>%
    tab_add_hline(at.row = top.black.line, row.side = "top", linewidth = 1.5)  %>%
    tab_add_hline(at.row = bottom.black.line, row.side = "bottom", linewidth = 1.5) %>%
    tab_add_title(text = subtitle, face = "plain", size = 10) %>%
    tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line")) %>%
    tab_add_footnote(text = footnote, size = 10, face = "italic") %>%
    tab_add_vline(at.column = vline, column.side = "left", from.row = 1, linetype = 2)

}




#' Find the column name of maximum/minimum value for each row
#'
#' @param df A dataframe
#' @param max Default FALSE
#' @param min Default FALSE
#' @param ties.method "random", "first", "last"
#' @param specified.column Numeric vector. Which columns want to be checked
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' df = LIRI
#' res <- loonR::findMaxMinColumnNamesForEachRow(df, min = T, specified.column = 3:6)
#' res
findMaxMinColumnNamesForEachRow <- function(df, max = FALSE, min = FALSE, ties.method = c("random", "first", "last"), specified.column = NULL){
  tmp.df = df

  ties.method = match.arg(ties.method)
  res <- data.frame(tmp.df, stringsAsFactors = FALSE, check.names = FALSE)

  if(!is.null(specified.column)){
    tmp.df = tmp.df[,specified.column]
  }else{
    specified.column = 1:ncol(tmp.df)
  }

  if(max){
    res$Max.ColID = max.col(tmp.df, ties.method = ties.method )
    res$Max.ColName = colnames(tmp.df)[res$Max.ColID]

    #use raw ID
    res$Max.ColID = specified.column[res$Max.ColID]
  }

  if(min){
    res$Min.ColID = max.col(tmp.df * -1, ties.method = ties.method )
    res$Min.ColName = colnames(tmp.df)[res$Min.ColID]

    #use raw ID
    res$Min.ColID = specified.column[res$Min.ColID]
  }

  res

}


#' Splite vector and return the specific element
#'
#' @param vector
#' @param sep
#' @param index
#'
#' @return
#' @export
#'
#' @examples
splitCharacter <- function(vector, sep = NULL, index = 1 ){

  if(is.null(sep)){
    stop("Pls set sep option")
  }
  res <- stringr::str_split(x, sep, simplify = TRUE)[,index]

  res = list(res=as.vector( unlist(res) ), raw = data.frame(res))
  res

}




