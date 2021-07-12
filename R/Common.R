

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
  if(!require(export)){
    devtools::install_github("tomwenseleers/export")
  }
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
plotJitterBoxplot <- function(values, group, title = "", xlab = "Group", ylab = "Value", color = "aaas", comparisons = '', method = "wilcox.test", label.y = NULL, add = "jitter"){
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
#'
#' @return ggplot2 object
#' @export
#'
#' @examples loonR::drawScatter(sample.info$before_reads, sample.info$mirdeep2_mapped)
drawScatter <- function(xvalue, yvalue, xlab = "X", ylab = "Y", group = NA,
                        color = "jco", title = "",
                        margin = TRUE, xlim = NULL, ylim = NULL,
                        show.sample.name = FALSE, label = NULL, cor.coef = F ){

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
                    palette = color, xlim=xlim, xlab = xlab, ylab = ylab,
                    title = title,
                    cor.coef = cor.coef,
                    cor.coeff.args = list(method = "pearson",
                                          label.x.npc = "left",
                                          label.y.npc = "top")
            )


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
  res <- data.frame(sapply(data.frame(df, check.names = F), function(x) as.numeric(as.character(x))), check.names = FALSE)
  rownames(res) <- row.n
  colnames(res) <- col.n
  res

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
#' @param cleveland
#' @param lollipop
#' @param value.name xlab, Default "Value"
#' @param sample.name ylab, Default "Name"
#'
#' @return
#' @export
#'
#' @examples
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
               ylab = ylab,
               ggtheme = theme_pubr()                        # ggplot2 theme
    )  +  theme_cleveland()                                      # Add dashed grids
  }else if(lollipop){
    p <- ggdotchart(dfm, x = "Name", y = "Value",
                    color = ifelse(anyNA(group), "black", "Group"),                                # Color by groups
                    # c("#00AFBB", "#E7B800", "#FC4E07")
                    palette = ifelse(anyNA(group)&palette=="aaas","#00AFBB", palette) , # Custom color palette
                    sorting = "ascending",                        # Sort value in descending order
                    add = "segments",                             # Add segments from y = 0 to dots
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



#' Title
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
#' res <- splitGroupByCutoff(group, value, fun="mean", cut.label = c("L","H"), group.prefix = "G", specific.group = c(1,2))
#' table(res$New.Label)
splitGroupByCutoff <- function(group = "", values = NULL, fun = NULL, quantile.cutoff = NULL, sample.names = NA,
                               cut.point = NULL, cut.label = NULL, specific.group = NULL, group.prefix = NULL){

  if(is.null(values)){
    stop("Please input a vector including values")
  }
  if(is.null(cut.label)){
    stop("Please input labels after spliting")
  }

  data.df <- data.frame(Group = group, Value = values, Label = "",
                        check.names = F, stringsAsFactors = F)

  if(!is.na(sample.names)){
    row.names(data.df) <- sample.names
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

    }

    data.df$Label <- as.character(data.df$Label)

  }


  if(!is.null(group.prefix)){
    data.df$Group = paste(group.prefix, data.df$Group, sep="")
  }


  data.df$Label[data.df$Label!=""] = paste("-", data.df$Label[data.df$Label!=""], sep="")

  data.df$New.Label = paste(data.df$Group, data.df$Label, sep="")

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
              suffixes = c(sfx[i],sfx[i+1]),
              by = by
           )
  }

  res

}







