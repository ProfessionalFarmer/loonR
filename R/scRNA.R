#' Read data: raw to PCA
#'
#' @param min.cells 10
#' @param min.features 200
#' @param gene.column 1
#' @param paths
#' @param sample_names
#' @param max_nCount_RNA
#' @param max_nFeature_RNA 8000
#' @param max_percent.mt 20
#' @param integrate_CCA Default FALSE
#' @param top_variable_features 2000
#' @param remove_doublet
#' @param remove_cc
#' @param remove_mt
#' @param integrate_harmony
#' @param batch Default 1
#' @param harmony_batch_colName Column names from batch correction
#' @param find_marker Specify detail DE method: e.g. wilcox_limma
#'
#' @return res
#' @export
#'
#' @examples
load10X <- function(paths = NA, sample_names = NA, min.cells = 10,
          min.features = 200, gene.column = 1, remove_doublet = T, batch = NULL,
          max_nCount_RNA = 10000, max_nFeature_RNA = 8000, max_percent.mt = 20,
          integrate_CCA = FALSE, top_variable_features = 2000, remove_cc = T, remove_mt = F,
          integrate_harmony = F, harmony_batch_colName = "orig.ident",
          find_marker = NULL, dims = 1:30){

  library(Seurat)
  library(dplyr)

  Object_list_raw = list()
  Object_list_filter = list()

  vars.to.regress = NULL
  if(remove_cc){
    vars.to.regress = c("S.Score","G2M.Score")

  }
  if(remove_mt){
       if(is.null(vars.to.regress)){
         vars.to.regress = c("percent.mt")
       }else{
         vars.to.regress = c(vars.to.regress,"percent.mt" )
       }
  }

  if(is.null(batch)){
    batch = rep(1, length(sample_names))
  }


  if(sum(is.na(paths)|is.na(sample_names))!=0){
    stop("Should provide file path and sample name")
  }
  if(length(paths) != length(sample_names)){
    stop("Should equal length")
  }

  for (i in 1:length(paths)){

    cat("\n", as.character( Sys.time() ),"\n--------------------------\nload", sample_names[i],"\n")

    data <- Read10X(paths[i], gene.column = gene.column)

    data <- CreateSeuratObject(data, min.cells = min.cells, min.features = min.features)

    data[["sample"]] = sample_names[i]
    data[["orig.ident"]] = sample_names[i]
    data[["batch"]] = batch[i]

    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")

    ############### filtering
    data.filt <- subset(data,
                        subset= nCount_RNA > min.features & nCount_RNA < max_nCount_RNA &
                          nFeature_RNA > min.features & nFeature_RNA < max_nFeature_RNA &
                          percent.mt < max_percent.mt )




    data.filt = loonR::norm_find_scale(data.filt, top_variable_features = top_variable_features)

    #### calculate cell cycle
    cat("\n--------------------------\nCalculate cell cycle" ,"\n")
    data.filt <- CellCycleScoring( object = data.filt, g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes )

    # data.filt <- RunPCA(data.filt) %>% RunUMAP(dims = 1:30) %>%  RunTSNE(dims = 1:30)
    # not run tsne to save time
    data.filt <- RunPCA(data.filt) %>% RunUMAP(dims = dims)


    if(remove_doublet){
      cat("\n--------------------------\nremove doublet", sample_names[i],"\n")

      data.filt = loonR::Find_doublet(data.filt, sct = F)
      # 提取判定为单胞的细胞进行下游分析
      data.filt <- subset(data.filt, subset=doublet_info=="Singlet")
    }

    Object_list_raw[[sample_names[i]]] = data

    Object_list_filter[[sample_names[i]]] = data.filt

  }
  raw.data = merge(Object_list_raw[[1]], Object_list_raw[2:length(Object_list_raw)])

  last.assay.name = "RNA" # CCA之后assay为integrated，harmony还是RNA
  reduction = "pca"

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ###### integrate Sample based on Object_list_filter
  if(integrate_CCA){
    cat("\n--------------------------\nPerform CCA integration\n")

    # normalize and identify variable features for each dataset independently
    filter.data <- lapply(X = Object_list_filter, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x,
                                selection.method = "vst",
                                nfeatures = top_variable_features)
    })


    features <- SelectIntegrationFeatures(object.list = filter.data)
    anchors <- FindIntegrationAnchors(object.list = filter.data, anchor.features = features)
    filter.data <- IntegrateData(anchorset = anchors)
    DefaultAssay(filter.data) <- "RNA"


    filter.data = FindVariableFeatures(filter.data,
                                       selection.method = "vst",
                                       nfeatures = top_variable_features)
    last.assay.name = "integrated"


  }else if(integrate_harmony){
    cat("\n", as.character( Sys.time() ),"\n--------------------------\nPerform SCTransfomr integration\n")

    filter.data = merge(Object_list_filter[[1]], Object_list_filter[2:length(Object_list_filter)])

    filter.data <- SCTransform(filter.data, vars.to.regress = vars.to.regress, verbose = FALSE)

    last.assay.name = "SCT"

    filter.data <- RunPCA(filter.data, assay = last.assay.name)

    #filter.data <- RunPCA(filter.data, features = VariableFeatures(filter.data) )
    cat("\n", as.character( Sys.time() ),"\n--------------------------\nPerform harmony integration\n")

    filter.data <- harmony::RunHarmony(filter.data, group.by.vars = harmony_batch_colName,
                                       reduction = "pca", assay.use = "SCT",
                                       reduction.save = "harmony", lambda = rep(1, length(harmony_batch_colName)) )

    last.assay.name = "SCT"
    reduction = "harmony"

  }else{
    cat("\n", as.character( Sys.time() ),"\n--------------------------\nMerge samples\n")

    filter.data = merge(Object_list_filter[[1]], Object_list_filter[2:length(Object_list_filter)])
    last.assay.name = "RNA"

  }

  rm(Object_list_filter, Object_list_raw)


  ########### 重新跑normalize  findvariable scale sct 的流程
  cat("\n--------------------------\nRe-scaling data \n")
  if(integrate_harmony){
    # joinlayer 不适用SCT
    filter.data <- JoinLayers(filter.data, assay = "RNA")
  }
  if(!integrate_harmony){

    filter.data <- JoinLayers(filter.data)
    # !integrate_harmony  SCT和harmony后，不用再做了
    cat("\n--------------------------\nScale data", "\n")
    filter.data = loonR::norm_find_scale(filter.data, assay.name = last.assay.name, vars.to.regress = vars.to.regress,
                                         top_variable_features = top_variable_features)

    # harmony之后不用做PCA了
    filter.data <- RunPCA(filter.data, assay = last.assay.name)
  }


  cat("\n", as.character( Sys.time() ),"\n--------------------------\nFind Cluster\n")

  DefaultAssay(filter.data) = last.assay.name

  filter.data <- FindNeighbors(filter.data, dims = dims, reduction = reduction, assay = last.assay.name )

  #计算SNN
  filter.data <- FindClusters(
    object = filter.data,
    resolution = c(seq(.1, 1, .2))
  )

  ##### umap
  cat("\n", as.character( Sys.time() ),"\n--------------------------\nRun UMAP \n")
  filter.data <- RunUMAP(filter.data, dims = dims, reduction = reduction, assay = last.assay.name )
  ##### tsne
  cat("\n", as.character( Sys.time() ),"\n--------------------------\nRun TSNE \n")
  filter.data <- RunTSNE(filter.data, dims = dims, reduction = reduction, assay = last.assay.name)


  # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
  res = list(raw = raw.data, filter = filter.data)
  DefaultAssay(filter.data) = "RNA"

  if(!is.null(find_marker)){
    Idents(filter.data) = as.vector(unlist(filter.data[[paste(last.assay.name,"_snn_res.0.1", sep = "")]]))
    if(is.logical(find_marker)){
      stop("Pls specify detailed method for FindAllMarkers")
    }
    cat("\n", as.character( Sys.time() ),"\n--------------------------\nRun Find Markers \n")
    #future::plan(strategy = "multicore", workers = 30)
    res$markers = FindAllMarkers(filter.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox_limma", assay = "RNA")
    #future::plan(strategy = "multicore", workers = 1)
  }

  cat("\n", as.character( Sys.time() ),"\n--------------------------\nDone \n")

  res
}




#' Find doublet
#'
#' @param data
#' @param PCs
#' @param sct
#' @param Doubletrate
#'
#' @return
#' @export
#'
#' @examples
Find_doublet <- function(data, PCs = 1:20, sct =FALSE, Doubletrate = 0.05, base_cluster = "RNA_snn_res.0.1"){
  library(DoubletFinder)
  # 寻找最优pk值
  sweep.res.list <- paramSweep(data, PCs = PCs, sct = sct, num.cores = 10) # 若使用SCT方法 标准化则'sct=T'
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))

  # 期望doublet数量
  homotypic.prop <- modelHomotypic(data@meta.data[[base_cluster]]) #可使用注释好的细胞类型
  Doubletrate <- Doubletrate
  nExp_poi <- round(Doubletrate*ncol(data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # 鉴定doublets
  data <- doubletFinder(data, PCs = PCs, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}

#' NormalizeData -> FindVariableFeatures -> ScaleData
#'
#' @param seurat.obj
#' @param assay.name
#' @param top_variable_features
#'
#' @return
#' @export
#'
#' @examples
norm_find_scale <- function(seurat.obj, assay.name = "RNA", top_variable_features = 2000, vars.to.regress = NULL){

  seurat.obj <- NormalizeData(seurat.obj, assay = assay.name)
  seurat.obj <- FindVariableFeatures(seurat.obj,
                                      selection.method = "vst",
                                      nfeatures = top_variable_features,
                                      assay = assay.name )

  seurat.obj <- ScaleData(seurat.obj, vars.to.regress = vars.to.regress,
                           features = VariableFeatures(seurat.obj),
                           assay = assay.name )

  seurat.obj
}


#' SCTranform HarmonyIntegration FindCluster
#'
#' @param seurat.obj
#' @param vars.to.regress colnames to remove covariate
#' @param group.by.vars colnames to remove batch effect
#'
#' @return
#' @export
#'
#' @examples
SCT_Harmony_cluster <- function(seurat.obj, vars.to.regress = NULL, group.by.vars = "orig.ident", dims = 1:20, runSCT = T, runHarmony = T){
  library(Seurat)
  library(dplyr)

  #future::plan(strategy = "multicore", workers = 50)
  #options(future.globals.maxSize = 3000 * 1024^2) # 3G

  assay = "RNA"
  reduction = "pca"
  if(runSCT){
    seurat.obj = seurat.obj %>% SCTransform(vars.to.regress = vars.to.regress) %>%
      RunPCA(assay = "SCT")
      assay = "SCT"
  }
  if(runHarmony){
    seurat.obj = seurat.obj %>%
      harmony::RunHarmony(group.by.vars = group.by.vars,
                          reduction = "pca", assay.use = "SCT",
                          reduction.save = "harmony", lambda = rep(1, length(group.by.vars)) )
    reduction = "harmony"
  }
  seurat.obj = seurat.obj %>%
    FindNeighbors(dims = dims, reduction = reduction, assay = assay ) %>%
    FindClusters( resolution = c(seq(.1, 1, .2))) %>%  #计算SNN
    RunUMAP(dims = dims, reduction = reduction, assay = assay )  %>%
    RunTSNE(dims = dims, reduction = reduction, assay = assay)

  #future::plan(strategy = "multicore", workers = 1)

  seurat.obj
}




#' Similar to dotplot, use violin plot
#'
#' @param object Seurat obj
#' @param groupBy Column names of identity
#' @param MarkerSelected marker gene names. vector
#' @param marker.group marker gene groups. Correspond with vector
#' @param color palatte
#'
#' @return
#' @export
#'
#' @examples
scViolinPlot <- function(object, groupBy, MarkerSelected, marker.group, color = NULL) {

  if(is.null(color)){
    color = c(pal_npg(alpha = 0.6)(10),pal_jama(alpha = 0.6)(7), pal_lancet(alpha = 0.6)(9))
  }
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  # https://zhuanlan.zhihu.com/p/536774748
  # (1)获取绘图数据1
  plot_data = FetchData(object = object,
                        vars = c(MarkerSelected, groupBy),
                        layer = 'data') %>%
    dplyr::rename(group = as.name(groupBy)) %>%
    tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')

  # (2)获取绘图数据2
  ident_plot = data.frame(gene = MarkerSelected, cluster = marker.group)

  # (3)绘图
  figure_1 = ggplot(data = plot_data, mapping = aes(x = Expr,
                                                    y = fct_rev(factor(x = Feat,
                                                                       levels = MarkerSelected)),
                                                    fill = group,
                                                    label = group)) +
    geom_violin(scale = 'width', adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
    facet_grid(cols = vars(group), scales = 'free') +
    cowplot::theme_cowplot(font_family = 'Arial') +
    scale_fill_manual(values = color) +
    xlab('Expression Level') +
    ylab('') +
    theme(legend.position = 'none',
          panel.spacing = unit(x = 0, units = 'lines'),
          axis.line = element_blank(), #去除x和y轴坐标线(不包括axis tick)；
          panel.background = element_rect(fill = NA, color = 'black'),
          strip.background = element_blank(), #去除分页题头背景；
          strip.text = element_text(color = 'black', size = 10, family = 'Arial', face = 'bold'),
          axis.text.x = element_text(color = 'black', family = 'Arial', size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(color = 'black', family = 'Arial', size = 15),
          axis.ticks.x = element_line(color = 'black', lineend = 'round'),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(x = 0.1, units = 'cm'))

  figure_2 = ggplot(data = ident_plot, aes(x = 1,
                                           y = fct_rev(factor(x = gene, levels = MarkerSelected)),
                                           fill = cluster)) +
    geom_tile() +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = color) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    guides(fill = guide_legend(direction = 'vertical',
                               label.position = 'right',
                               title.theme = element_blank(),
                               keyheight = 0.5,
                               nrow = 2)) +
    xlab('Feature') +
    theme(legend.text = element_text(family = 'Arial', color = 'black', size = 11),
          legend.position = 'bottom',
          legend.justification = 'left',
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,05,0,0),
          panel.spacing = unit(0, 'lines'),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(x = c(0,0,0,0), units = 'cm'),
          axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, color = 'black', family = 'Arial'),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

  figure_2 + figure_1 + patchwork::plot_layout(nrow = 1, widths = c(0.03, 0.97))
}


#' convert seurat to monocle 3 obj
#'
#' @param seu.obj
#' @param partition.col.name
#' @param use_partition
#' @param auto.root if automatically identify root
#' @param root.cluster.name if auto.root is TRUE, set root cluster name
#' @param cluster.col.name if auto.root is TRUE, set cluster column name
#'
#' @return
#' @export
#'
#' @examples
Seurat2Monocle3 = function(seu.obj, partition.col.name = NULL, use_partition = F, auto.root = T, root.cluster.name = NULL, cluster.col.name = NULL){

  if(!require(monocle3)){
    devtools::install_github('cole-trapnell-lab/monocle3')
    require(monocle3)
  }
  if(!require(SeuratWrappers)){
    remotes::install_github('satijalab/seurat-wrappers')
  }

  seu.obj.mnc = SeuratWrappers::as.cell_data_set(seu.obj)

  # https://www.jianshu.com/p/afbef525a03b
  #mapping cluster
  #cds对象的cds2@clusters$UMAP$clusters是一个named vector。
  seu.obj.mnc@clusters$UMAP$clusters <- Idents(seu.obj)[rownames(colData(seu.obj.mnc))]

  #parition的数目根据自己的实际情况处理，如果后面的 learn_graph(cds2, use_partition = F)的use_partition为F，则不用partition来做不同pseudotime推测，没有没有必要设置partition；
  if(is.null(partition.col.name)){
    seu.obj.mnc@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(seu.obj.mnc)))), levels = 1)
  }else{
    seu.obj.mnc@clusters$UMAP$partitions <- factor(seu.obj@meta.data[rownames(colData(seu.obj.mnc)), partition.col.name])

  }


  names(seu.obj.mnc@clusters$UMAP$partitions) <- rownames(colData(seu.obj.mnc))

  # 从seurat obj transfer过去的cds对象没有做estimate size factor参数
  ## Calculate size factors using built-in function in monocle3
  seu.obj.mnc <- estimate_size_factors(seu.obj.mnc)

  #创建cds对象的时候，cds要求gene_meta这gedf必须要有一列为gene_short_name
  ## Add gene names into CDS

  seu.obj.mnc@rowRanges@elementMetadata@listData$gene_short_name <- rownames(seu.obj.mnc)

  #如果出现了order_cells(cds)报错"object 'V1' not found"，否则不执行
  #来源：https://github.com/cole-trapnell-lab/monocle3/issues/130](https://github.com/cole-trapnell-lab/monocle3/issues/130

  # # # # # # # # # # # # # # # # # # # #
  rownames(seu.obj.mnc@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
  # colnames(cds@reducedDims$UMAP) <- NULL
  colnames(seu.obj.mnc@int_colData@listData$reducedDims@listData$UMAP) <- NULL
  # # # # # # # # # # # # # # # # # # # #

  seu.obj.mnc <- learn_graph(seu.obj.mnc, use_partition = use_partition)


  if(auto.root){
    if(is.null(root.cluster.name) | is.null(cluster.col.name)){
      stop("Pls set root name and column name when automaticlly find root")
    }

    get_earliest_principal_node <- function(cds, time_bin = root.cluster.name){
      cell_ids <- which(colData(cds)[, cluster.col.name] == time_bin)

      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]

      root_pr_nodes
    }
    seu.obj.mnc <- order_cells(seu.obj.mnc, root_pr_nodes=get_earliest_principal_node(seu.obj.mnc))

  }else{
    seu.obj.mnc <- order_cells(seu.obj.mnc, reduction_method = "UMAP")
  }

  # https://cloud.tencent.com/developer/article/1819550
  ##寻找拟时轨迹差异基因
  #graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
  #空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
  Track_genes <- graph_test(seu.obj.mnc, neighbor_graph="principal_graph", cores=10)
  #挑选top10画图展示
  Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
    pull(gene_short_name) %>% as.character()
  #基因表达趋势图
  plot_genes_in_pseudotime(seu.obj.mnc[Track_genes_sig,], color_cells_by="predicted.id",
                           min_expr=0.5, ncol = 2)
  #FeaturePlot图
  plot_cells(seu.obj.mnc, genes=Track_genes_sig, show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,  label_leaves=FALSE)
  ##寻找共表达模块
  genelist <- pull(Track_genes, gene_short_name) %>% as.character()
  gene_module <- find_gene_modules(seu.obj.mnc[genelist,], resolution=1e-2, cores = 10)
  cell_group <- tibble::tibble(cell=row.names(colData(seu.obj.mnc)),
                               cell_group=colData(seu.obj.mnc)$predicted.id)
  agg_mat <- aggregate_gene_expression(seu.obj.mnc, gene_module, cell_group)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")

  # plot_cells(t.cell.mnc, color_cells_by = "pseudotime", label_cell_groups=F, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5)
  res = list(monocle.obj = seu.obj.mnc,
             Track_genes = Track_genes,
             gene_module, cell_group, agg_mat
  )
  res
}



export2CellphoneDB <- function(seu.obj, cell.type.column= NULL, dir = NULL, diff.gene = NULL){

  if(is.null(cell.type.column)){
      stop("Input cell type column")
  }

  if(is.null(dirx)){
    stop("Input direcotry")
  }

  if(!dir.exists(dir)){
      dir.create(dir)
  }

    # Save normalised counts - NOT scaled!
  #Matrix::writeMM(seu.obj@assays$RNA$data, file = paste0(dir,'/matrix.mtx') )
  write.table(seu.obj@assays$RNA$data, sep = "\t", quote = F, file = paste0(dir,'/matrix.tsv') )
  # save gene and cell names
  write(x = rownames(seu.obj@assays$RNA$data), file = paste0(dir, "/features.tsv")  )
  write(x = colnames(seu.obj@assays$RNA$data), file = paste0(dir, "/barcodes.tsv")  )

  meta <-  data.frame(Cell = rownames(seu.obj@meta.data),
                     cell_type = myeloid.cell[[cell.type.column]] )
  write.table(meta, file = paste0(dir, 'meta.tsv'), sep = '\t', quote = F, row.names = F)

  if(!is.null(diff.gene)){
    fDEGs = diff.gene[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')]
    colnames(fDEGs)= c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_logFC', 'pct.1', 'pct.2')
    write.table(fDEGs, file = paste0(dir, '/DEGs.tsv'), sep = '\t', quote = F, row.names = F)
  }

}



#' Pseudo differential analysis
#'
#' @param mono.cds.obj
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
pre_pseudotime_matrix <- function(mono.cds.obj, genes = NULL){
  # https://github.com/cole-trapnell-lab/monocle-release/issues/295
  if(is.null(genes)){
    modulated_genes <- graph_test(mono.cds.obj,
                                  neighbor_graph = "principal_graph",
                                  cores = 40)
    genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
  }

  pt.matrix <- exprs(mono.cds.obj)[match(genes,rownames(rowData(mono.cds.obj))),order(pseudotime(mono.cds.obj))]
  #Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
  #normalized_counts(cds, norm_method = "log")

  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

  rownames(pt.matrix) <- genes
  pt.matrix
}


#' Pseudo heatmap
#'
#' @param mt from pre_pseudotime_matrix function
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
pseudoHeatmap <- function(mt, genes = NULL){
  library(ClusterGVis)

  # kmeans
  ck <- clusterData(exp = pt.matrix,
                    cluster.method = "kmeans",
                    cluster.num = 5)

  # add line annotation
  visCluster(object = ck,
             plot.type = "both",
             add.sampleanno = F,
             markGenes =  genes )

}



#' Export seurat object for SCENIC analysis
#'
#' @param seurat.obj
#' @param dir Directory for export object
#' @param prefix Default scenic. Prefix in dir
#'
#' @return
#' @export
#'
#' @examples
exportScenicLoom <- function(seurat.obj, dir = "scenic", prefix = "scenic"){
  library(SCENIC)
  library(Seurat)

    ## Get data from sce object:
  exprMat <- GetAssayData(object = seurat.obj, assay = "RNA", slot = "counts")
  # To use Seurat clusters as cell annotation (e.g. for visualization):
  cellInfo <- data.frame(seuratCluster=Idents(seurat.obj))

  if(!dir.exists(dir)){
    dir.create(dir)
  }

  if(!require("SCopeLoomR")){
    devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
  }
  require("SCopeLoomR")
  loom <- SCopeLoomR::build_loom(file.path(dir, paste0(prefix, ".loom") ), dgem = exprMat)
  loom <- SCENIC::add_cell_annotation(loom, cellInfo)
  SCopeLoomR::close_loom(loom)

  cat("Resource: ","https://github.com/aertslab/pySCENIC/tree/master/resources", "\n")
  cat("TF:", "https://resources.aertslab.org/cistarget/tf_lists/", "\n")
  cat("motif–TF: ", "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", "\n")
  cat("Database: ","https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/", "\n")




}



#' Run inferCNV
#'
#' @param seurat.obj
#' @param ref_group_names Reference cell type
#' @param gene_pos gene order file
#' @param save_local_path Path to save rds
#'
#' @return
#' @export
#'
#' @examples
runInferCNV <- function(seurat.obj, ref_group_names = NULL,
                        gene_pos = "~/ref/hg38/gencode.v44/human_genes_pos.txt",
                        save_local_path = NULL){

  if(is.null(ref_group_names)){
    stop("Provide ref group")
  }

  if(!file.exists(gene_pos)){
    stop("Gene order/position file not exist. \n")
  }


  library(Seurat)
  library(infercnv)

  levels(seurat.obj)
  ##筛选分析需要的细胞类型
  #seurat.obj <- subset(seurat.obj, idents=c('Epithelial cells', 'Myeloid cells', 'T cells'))

  #抽样，仅用于该教程
  # counts <- GetAssayData(seurat.obj, slot = 'counts')
  counts <- GetAssayData(seurat.obj, assay = "RNA", "counts")

  # 保持gene order和count的行名一致
  # 位置信息可以来自ensembl，也可以来自refSeq （https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/）
  # gene_order = read.delim(gene_pos, sep = "\t", header = F)
  # overlap_gene = intersect(gene_order$V1, rownames(counts))
  # # setbset counts
  # counts = counts[overlap_gene,]
  # gene_order = gene_order[match(overlap_gene, gene_order$V1), ]
  # tmp.file = tempfile()
  # gene_pos = tmp.file
  # write.table(gene_order, file = tmp.file, sep = "\t", col.names = F, row.names = F, quote = F)

  # get annotation
  anno <- data.frame(Idents(seurat.obj))

  #先把脚本下载下来，再处理gtf文件，生成human_genes_pos.txt，https://github.com/broadinstitute/infercnv/tree/master/scripts
  #python ./scripts/gtf_to_position_file.py --attribute_name "gene_name" your_reference.gtf human_genes_pos.txt

  # 创建inferCNV对象
  # https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                      annotations_file = anno,
                                      delim="\t",
                                      gene_order_file = gene_pos,
                                      min_max_counts_per_cell = c(100, +Inf),
                                      ref_group_names = ref_group_names)

  # run inferCNV
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = "./infercnv",
                               cluster_by_groups = T,
                               HMM = T,
                               denoise = TRUE,
                               num_threads = 40)

  if(!is.null(save_local_path)){
    warning("save data to local file: ", save_local_path)
    saveRDS(infercnv_obj, file = save_local_path)
  }else{
    infercnv_obj
  }

}




#' Run cytoTRACE 1
#'
#' @param seurat.obj
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
runCytoTrace1 <- function(seurat.obj, cores = 10){

  library(CytoTRACE)

  counts <- GetAssayData(seurat.obj, assay = "RNA", "counts")
  results <- CytoTRACE(data.frame(counts), ncores = cores, subsamplesize = 1000)
  seurat.obj@meta.data$CytoTRACE = results$CytoTRACE[Cells(seurat.obj)]
  seurat.obj
}

#' Run cytoTRACE 2
#'
#' @param seurat.obj
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
runCytoTrace2 <- function(seurat.obj, cores = 20){
  # https://github.com/digitalcytometry/cytotrace2/tree/main?tab=readme-ov-file
  library(CytoTRACE2)

  seurat.obj <- cytotrace2(seurat.obj, is_seurat = TRUE,
                           slot_type = "counts", species = "human", seed = 14,
                           batch_size = 10000, parallelize_models = T,
                           ncores = cores)
  seurat.obj
}


