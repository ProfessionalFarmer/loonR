
#' Plot heatmap with log 2 fold change and risk probability information
#'
#' @param heatmap.df Row: miRNA, Column: Sample
#' @param label A factor with labels TURE/FALSE
#' @param risk.pro
#' @param lgfold corresponded to row name sequence
#' @param group.name Default "Cancer"
#' @param scale Default "TRUE"
#' @param ylim = c(0, 1), risk score range.
#' @param show.lgfold = TRUE。 Whether to show right panel.
#' @param show.risk.pro
#' @param bar.name
#' @param height
#' @param show_column_names Default False
#' @param cluster_rows
#' @param z.score.cutoff Default 2
#' @param cluster_columns
#' @param specified.color Default c("#0c3e74","#77a8cd","white","#d86652","#7e0821") or colorRampPalette(c("navy", "white", "firebrick3"))(50)
#'
#' @return A heatmap plot by complex heatmap
#' @export
#'
#' @examples heatmap.with.lgfold.riskpro(data.tmp[candi,],label, logfd,  risk.pro)
heatmap.with.lgfold.riskpro <- function(heatmap.df, label, risk.pro, lgfold=NA, scale=TRUE, group.name="Cancer", bar.name = "Log2FC", ylim = c(0, 1),
                                        show.lgfold = TRUE, show.risk.pro = TRUE, height = 5, show_column_names = FALSE, cluster_rows = FALSE,
                                        cluster_columns = FALSE, z.score.cutoff = 2, specified.color = c("#0c3e74","#77a8cd","white","#d86652","#7e0821") ){
  if(!require(ComplexHeatmap)){
    BiocManager::install("ComplexHeatmap")
  }



  if (anyNA(lgfold)){
    show.lgfold = FALSE
    lgfold = replicate(nrow(heatmap.df),1)
  }

  label = factor(label, levels = unique(label))

  if(scale){
    heatmap.df = t(scale(t(heatmap.df)))
    heatmap.df[ heatmap.df > z.score.cutoff] <- z.score.cutoff
    heatmap.df[ heatmap.df < -z.score.cutoff] <- -z.score.cutoff
  }

  library(ComplexHeatmap)
  if(show.risk.pro){
    # 根据risk score排序
    heatmap.df <- cbind(
      heatmap.df[, label==levels(label)[1] ][, order(risk.pro[ label==levels(label)[1] ] ) ],
      heatmap.df[, label==levels(label)[2] ][, order(risk.pro[ label==levels(label)[2] ] ) ]  )

    # 对riskscore排序
    risk.pro <- c(
      risk.pro[label==levels(label)[1] ][order(risk.pro[label==levels(label)[1] ] ) ],
      risk.pro[label==levels(label)[2] ][order(risk.pro[label==levels(label)[2] ] ) ]  )
  }else{
    # 根据分组排序
    heatmap.df <- heatmap.df[,c(which(label==levels(label)[1]), which(label==levels(label)[2]))]
    label <- label[c(which(label==levels(label)[1]), which(label==levels(label)[2]))]
  }


  # heatmap和barplot一起画

  #row_ha = rowAnnotation( assign(eval(bar.name), anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))   )
  row_ha = rowAnnotation( Log2FC = anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))

  label = factor(label)
  Tumor = loonR::get.palette.color("jama_classic", n = 2)
  names(Tumor) = levels(label)



  # rename annotation names
  annotation <- data.frame(Tmp = label[ c(which(label==levels(label)[1]),
                                          which(label==levels(label)[2]) )
  ]
  )
  colnames(annotation) <- group.name

  ann_colors = list(Tmp = Tumor)
  names(ann_colors) = group.name


  if(show.risk.pro){
    ha = HeatmapAnnotation(df = annotation,
                           col = ann_colors,
                           Risk = anno_points(risk.pro, pch = 16, size = unit(1, "mm"),
                                              gp = gpar(col = "black"),
                                              ylim = ylim,
                                              axis_param = list( side = "left", at = ylim, labels = as.character(ylim) )
                           )
    )
  }else{
    ha = HeatmapAnnotation(df = annotation,
                           col = ann_colors    )
  }

  #
  if(show.lgfold){

    Heatmap(heatmap.df, col = specified.color,
            name = " ", cluster_rows = cluster_rows, cluster_columns = cluster_columns,
            show_row_names = TRUE, show_column_names = show_column_names, height = unit(height, "cm"),
            top_annotation = ha,
            right_annotation = row_ha  )

  }else{

    Heatmap(heatmap.df, col = specified.color,
            name = " ", cluster_rows = cluster_rows, cluster_columns = cluster_columns,
            show_row_names = TRUE, show_column_names = show_column_names, height = unit(height, "cm"),
            top_annotation = ha  )

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
#'                            stringsAsFactors = FALSE)
#' #Gender.col = loonR::get.palette.color("aaas",2)
#' #names(Gender.col) = unique(annotation.df$Gender)
#' Gender.col = c("Male"="blue", "Female"="red")
#'
#' col_fun = circlize::colorRamp2(c(50, 100), c("blue", "red")) # continuous variable
#'
#' annotation.color = list(Gender=Gender.col, Age=col_fun)
#'
#' heatmap.annotation(group = group, annotation.df = annotation.df,
#'                    annotation.color = annotation.color, sort.group = T)
heatmap.annotation <- function(group = NULL, annotation.df = NULL, annotation.color = NULL, sort.group=TRUE, title = "",
                               group.color = "aaas", group.prefix ="", column_split = TRUE ){

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
    annotation_name_side = "left"
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



#' Perform hypergeometric test
#'
#' @param row.group
#' @param col.group
#' @param row.prefix
#' @param col.prefix
#' @param lower.tail Default FALSE
#' @param title
#' @param log10.lowest Default 3. Minimux or maximum log10 value. Useful when meet inf or draw heatmap
#' @param print.fig Default print the figure
#'
#' @return list(pval.df.log, pval.df, plot)
#' @export
#'
#' @examples
#' g1 = sample(c("G1","G2","G3"),10, replace = T)
#' g2 = sample(c("1","2","3"),10, replace = T)
#' loonR::hyperGeoTest(g1, g2, col.prefix = "E")
hyperGeoTest <- function(row.group, col.group, row.prefix = "", col.prefix = "", lower.tail = FALSE, title = "", log10.lowest = 5, print.fig = TRUE){

  # https://www.omicsclass.com/article/324
  # 1-phyper(抽取样本中属于“特定类别”的数量-1,总样本中“特定类别”的数量, 总样本数-总样本中“特定类别”的数量, 从总样本中随机抽取的数量,)

  # 超几何检验，与原来的分组比较
  geomatrix  = unclass(table(row.group, col.group))
  colnames(geomatrix) = paste(col.prefix, colnames(geomatrix),sep="" )
  rownames(geomatrix) = paste(row.prefix, rownames(geomatrix),sep="" )

  # perform geometrix,把p值放在相同矩阵的数据框中
  tmpgeo = matrix(nrow=length(row.names(geomatrix)),ncol=length(colnames(geomatrix)))
  colnames(tmpgeo) = colnames(geomatrix)
  rownames(tmpgeo) = rownames(geomatrix)


  for(i in 1:length(row.names(tmpgeo))  ){ # row
    for(j in 1:length(colnames(tmpgeo))){  # column
      # 白球的个数，白球的总个数，黑球的总个数，抽的球（不是黑球，是球）个数
      p = phyper(geomatrix[i,j]-1, sum(geomatrix[i,]), sum(geomatrix)-sum(geomatrix[i,]), sum(geomatrix[,j]), lower.tail = lower.tail   )
      tmpgeo[i,j] = p
    }
  }
  tmpgeo = as.data.frame(tmpgeo)

  tmpgeo.log = -log10(tmpgeo)
  tmpgeo.log[tmpgeo.log > log10.lowest ] = log10.lowest - 0.1


  # reture formated value
  res <- list()
  res$pval.df.log <- tmpgeo.log

  tmpgeo <- loonR::convertDfToNumeric(tmpgeo)
  res$pval.df.notFormated <- tmpgeo

  tmpgeo <- apply(tmpgeo, 2, function(x) {
    scales::scientific(x, digits = 3)
  })


  res$pval.df <- format(tmpgeo, trim = TRUE, digits = 3, scientific = 3)
  res$table <- geomatrix

  # # # # # # # # # # # # # # # # # # Option 1
  # pheatmap::pheatmap(tmpgeo.log,
  #                    cluster_rows = F,
  #                    cluster_cols = F,
  #                    color = c (rep("#FFFFFF", 26), colorRampPalette(c("#FFFFFF", "#0269A4" ))(78) ),
  #                    breaks=unique(c(seq(0,5, length=100-1  ))),
  #                    display_numbers = format(tmpgeo, trim = TRUE, digits = 3, scientific = 3),
  #                    main = title
  # )

  # # # # # # # # # # # # # # # # # # Option 2
  col_fun = circlize::colorRamp2(c(0, -log10(0.05)-0.1, log10.lowest-0.5),
                                 c("#FFFFFF", "#FFFFFF", "#0269A4")  )

  library(grid)
  res$plot = ComplexHeatmap::Heatmap(tmpgeo.log,
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     col = col_fun,
                                     heatmap_legend_param = list(title = "-Log10(P)", color_bar = c("continuous")  ),
                                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                       grid.text(res$pval.df[i, j], x, y)
                                     },
                                     row_names_side = "left",
                                     column_names_rot = 45,
                                     column_names_side = "top",
                                     rect_gp = grid::gpar(col = "#c1c1c1", lwd = 2),
                                     column_title = title,
                                     row_title = character(0)

  )

  if(print.fig){
    print(res$plot)
  }


  # return res
  res

}












