
#' Plot heatmap with log 2 fold change and risk probability information
#'
#' @param heatmap.df Row: miRNA, Column: Sample
#' @param label A factor with labels TURE/FALSE
#' @param risk.pro
#' @param lgfold corresponded to row name sequence
#' @param group.name Default "Cancer"
#' @param scale Default "TRUE"
#' @param ylim = c(0, 1), risk score range.
#' @param show.lgfold = TRUE. Whether to show right panel.
#' @param show.risk.pro
#' @param bar.name
#' @param height
#' @param show_column_names Default False
#' @param cluster_rows
#' @param z.score.cutoff Default 2
#' @param cluster_columns
#' @param specified.color Default c("#0c3e74","#77a8cd","white","#d86652","#7e0821") or colorRampPalette(c("navy", "white", "firebrick3"))(50)
#' @param palette Color palette for group. Default "jama_classic"
#' @param show_row_names Defaut TRUE
#' @param cluster_within_group Cluster within group
#' @param annotation.df Annotation data.frame
#'
#' @return A heatmap plot by complex heatmap
#' @export
#'
#' @examples heatmap.with.lgfold.riskpro(data.tmp[candi,],label, logfd,  risk.pro)
heatmap.with.lgfold.riskpro <- function(heatmap.df, label, risk.pro=NA, lgfold=NA, scale=TRUE, group.name="Cancer", bar.name = "Log2FC", ylim = c(0, 1),
                                        show.lgfold = TRUE, show.risk.pro = TRUE, height = 5, annotation.df = NULL,
                                        show_column_names = FALSE, show_row_names = TRUE, cluster_rows = FALSE,
                                        cluster_columns = FALSE, cluster_within_group = FALSE, z.score.cutoff = 2, specified.color = c("#0c3e74","#77a8cd","white","#d86652","#7e0821"), palette = "jama_classic" ){

  if (!require(ComplexHeatmap)) {
    BiocManager::install("ComplexHeatmap")
  }


  if (anyNA(lgfold)) {
    show.lgfold <- FALSE
    lgfold <- replicate(nrow(heatmap.df), 1)
  }

  if(!is.factor(label)){
    label <- factor(label, levels = unique(label))
  }

  if(length(label)!=ncol(heatmap.df)){
    stop("Label number not the same as data.frame col length")
  }

  label.risk.df <- data.frame(
    label = label,
    risk.pro = risk.pro,
    index = (1:length(label))
  )
  if (!is.null(annotation.df)){ # if annotation variable is null
    label.risk.df = data.frame(label.risk.df, annotation.df)
  }


  if (scale) {
    heatmap.df <- t(scale(t(heatmap.df)))
    heatmap.df[heatmap.df > z.score.cutoff] <- z.score.cutoff
    heatmap.df[heatmap.df < -z.score.cutoff] <- -z.score.cutoff
  }

  library(ComplexHeatmap)
  if (show.risk.pro) {
    # 根据label和risk score同时排序
    label.risk.df = dplyr::arrange(label.risk.df, label, risk.pro)

    heatmap.df <- heatmap.df[,label.risk.df$index]

    # 对riskscore排序
    risk.pro <- label.risk.df$risk.pro

  } else {
    # 根据分组排序
    label.risk.df = dplyr::arrange(label.risk.df, label)

    heatmap.df <- heatmap.df[,label.risk.df$index]

  }


  # heatmap和barplot一起画

  # row_ha = rowAnnotation( assign(eval(bar.name), anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))   )
  row_ha <- rowAnnotation(Log2FC = anno_barplot(lgfold, gp = gpar(fill = "black", col = "black")))

  label <- label.risk.df$label

  Tumor <- loonR::get.palette.color(palette, n = length(unique(label.risk.df$label)))
  names(Tumor) <- levels(label)


  # rename annotation names
  annotation <- data.frame(Tmp = label.risk.df$label)
  colnames(annotation) <- group.name

  ann_colors <- list(Tmp = Tumor)
  names(ann_colors) <- group.name

  ## add new annotation,不能直接用annotation.df，因为已经排序了,要用label.risk.df
  if(!is.null(annotation.df)){
    annotation = data.frame(annotation,
                            label.risk.df[,4:ncol(label.risk.df)],
                            check.names = FALSE)
  }


  if (show.risk.pro) {
    ha <- HeatmapAnnotation(
      df = annotation,
      col = ann_colors,
      Risk = anno_points(risk.pro,
        pch = 16, size = unit(1, "mm"),
        gp = gpar(col = "black"),
        ylim = ylim,
        axis_param = list(side = "left", at = ylim, labels = as.character(ylim))
      )
    )
  } else {
    ha <- HeatmapAnnotation(
      df = annotation,
      col = ann_colors
    )
  }

  if(cluster_within_group){
    cluster_columns = cluster_within_group(heatmap.df, label)
  }


  # 2022-06-15
  row_names_side = "left"
  if(cluster_rows | show.lgfold){
    row_names_side = "right"
  }

  #
  if (show.lgfold) {
    Heatmap(heatmap.df,
      col = specified.color,
      name = " ", cluster_rows = cluster_rows, cluster_columns = cluster_columns,
      show_row_names = show_row_names, show_column_names = show_column_names, height = unit(height, "cm"),
      top_annotation = ha, row_names_side = row_names_side,
      right_annotation = row_ha
    )
  } else {
    Heatmap(heatmap.df,
      col = specified.color,
      name = " ", cluster_rows = cluster_rows, cluster_columns = cluster_columns,
      show_row_names = show_row_names, show_column_names = show_column_names, height = unit(height, "cm"),
      top_annotation = ha, row_names_side = row_names_side,
    )
  }

}



#' Only draw annotation
#'
#' @param group
#' @param annotation.df Column is variable, row is sample
#' @param annotation.color A list object, the name should be corresponded to df column names
#' @param sort.group Order the group
#' @param title
#' @param group.color
#' @param group.prefix
#' @param column_split Default TRUE
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(111)
#' group = sample(c(1,2),10,replace = T)
#' annotation.df = data.frame(Gender = sample(c("Male","Female"),10,replace = T),
#'                            Age = runif(10, min=50, max=80),
#'                            stringsAsFactors = FALSE))
#' #Gender.col = loonR::get.palette.color("aaas",2)
#' #names(Gender.col) = unique(annotation.df$Gender)
#' Gender.col = c("Male"="blue", "Female"="red")
#'
#' col_fun = circlize::colorRamp2(c(50, 100), c("blue", "red")) # continuous variable
#'
#' annotation.color = list(Gender=Gender.col, Age=col_fun)
#'
#' loonR::heatmap.annotation(group = group, annotation.df = annotation.df,
#'                    annotation.color = annotation.color, sort.group = T)
heatmap.annotation <- function(group = NULL, annotation.df = NULL, annotation.color = NULL, sort.group=TRUE, title = "",
                               group.color = "aaas", group.prefix ="", column_split = TRUE, na_col ="grey" ){

  if(!require(ComplexHeatmap)){
    BiocManager::install("ComplexHeatmap")
  }


  if(is.null(group) | is.null(annotation.df)){
    stop("Please sest group or set annotation df")
  }

  group = paste(group.prefix, group, sep="")
  zero_row_mat = matrix(nrow = 0, ncol = length(group))
  colnames(zero_row_mat) = names(group)

  annotation.df = data.frame(Group = group, annotation.df, check.names = F, stringsAsFactors = F)

  if(sort.group){
    od = order(group)

    # change the order
    group = group[od]
    zero_row_mat = zero_row_mat
    annotation.df = annotation.df[od,]
  }

  if(length(group.color)>1){
    if(length(unique(group))!=length(group.color)){
      stop("Please set the same number group.color with group")
    }
    names(group.color) = unique(group)
    annotation.color$Group = group.color

  }else{

    group.color = loonR::get.palette.color(group.color, length(unique(group)))
    names(group.color) = unique(group)
    annotation.color$Group = group.color

  }


  set.seed(111)

  ha = HeatmapAnnotation(
    df = annotation.df,
    col = annotation.color,
    annotation_name_side = "left", na_col = na_col
  )


  if(column_split){
    column_split = group
  }else{
    column_split = NULL
  }

  set.seed(111)

  Heatmap(zero_row_mat,
          top_annotation = ha,
          column_title = title,
          show_column_names = FALSE,
          column_split = column_split,
          heatmap_legend_param = list(ncol=3)
  )

}





#' Plot heatmap while perform word cloud analysis
#'
#' @param df row is sample, column is gene. Will perform gene correlation analysis instead of sample-wise
#'
#' @return
#' @export
#'
#' @examples
#' # https://jokergoo.github.io/simplifyEnrichment/articles/word_cloud_anno.html
#'
#' df = data.frame(`a b c dd` = c(1,2,1,1,1), `a d c bb a c` =  c(3,2,3,1,1), `bb aa cc dd` = c(4,2,5,7,6))
#' M = cor(df)
#' # show correlation
#' corrplot(M, method = 'number', order = 'hclust')
#'
#'
#'
#'
heatmap.with.wordcloud <- function(df){

  # todo: 手动聚类，cutlabel之后，自己做词云分析然后画图，需要决定cut成几个簇
  if(is.null(term)){
    term = colnames(df)
  }

  if(!require(ComplexHeatmap)){
    BiocManager::install("ComplexHeatmap")
    require(ComplexHeatmap)
  }
  if(!require(simplifyEnrichment)){
    devtools::install_github("jokergoo/simplifyEnrichment")
    require(simplifyEnrichment)
  }


  Heatmap(mat, row_split = split,
          right_annotation = rowAnnotation(wc = anno_word_cloud(split, term))
  )


}












