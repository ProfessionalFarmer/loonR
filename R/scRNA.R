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
#' @param find_marker Default NULL. Specify detail DE method from SeuratWrappers::RunPrestoAll: e.g. wilcox_limma
#' @param dims
#' @param max_precent.hb Default 10
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
          find_marker = NULL, dims = 1:30, max_precent.hb = 10){

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

    data[["percent.mt"]]   <- PercentageFeatureSet(data, pattern = "^MT-")
    data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
    data[["percent.hb"]]   <- PercentageFeatureSet(data, pattern = "^HB[^(P|E|S)]")

    ############### filtering
    data.filt <- subset(data,
                        subset= nCount_RNA > min.features & nCount_RNA < max_nCount_RNA &
                          nFeature_RNA > min.features & nFeature_RNA < max_nFeature_RNA &
                          percent.mt < max_percent.mt & percent.hb < max_precent.hb )



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
    res$markers = SeuratWrappers::RunPrestoAll(filter.data, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = find_marker, assay = "RNA")
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
Seurat2Monocle3 = function(seu.obj, use_partition = F, auto.root = T, root.cluster.name = NULL, cluster.col.name = NULL){

  if(!require(monocle3)){
    devtools::install_github('cole-trapnell-lab/monocle3')
    require(monocle3)
  }
  if(!require(SeuratWrappers)){
    remotes::install_github('satijalab/seurat-wrappers')
  }

  seu.obj.mnc = SeuratWrappers::as.cell_data_set(seu.obj)

  # https://www.jianshu.com/p/afbef525a03b
  # 自己做clustering，则不用mapping
  seu.obj.mnc <- cluster_cells(cds = seu.obj.mnc, reduction_method = "UMAP")

  #mapping cluster
  seu.obj.mnc@clusters$UMAP$clusters <- Idents(seu.obj)[rownames(colData(seu.obj.mnc))]

  # 从seurat obj transfer过去的cds对象没有做estimate size factor参数
  ## Calculate size factors using built-in function in monocle3
  seu.obj.mnc <- estimate_size_factors(seu.obj.mnc)

  #创建cds对象的时候，cds要求gene_meta这gedf必须要有一列为gene_short_name
  ## Add gene names into CDS

  seu.obj.mnc@rowRanges@elementMetadata@listData$gene_short_name <- rownames(seu.obj.mnc)

  # use_parition这里一直是true，但具体的partition在前面指定
  seu.obj.mnc <- learn_graph(seu.obj.mnc, use_partition = T)


  if(auto.root){
    if(is.null(cluster.col.name)){
      stop("Pls set root name and column name when automaticlly find root")
    }

    get_earliest_principal_node <- function(cds, time_bin = root.cluster.name){
      cell_ids <- which(colData(cds)[, "ident"] == time_bin)

      closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]

      root_pr_nodes
    }
    seu.obj.mnc <- order_cells(seu.obj.mnc, root_pr_nodes=get_earliest_principal_node(seu.obj.mnc, root.cluster.name))

  }else{
    seu.obj.mnc <- order_cells(seu.obj.mnc, reduction_method = "UMAP")
  }

  plot_cells(
    cds = seu.obj.mnc,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE
  )

  # https://cloud.tencent.com/developer/article/1819550
  ##寻找拟时轨迹差异基因
  #graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
  #空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
  Track_genes <- graph_test(seu.obj.mnc, neighbor_graph="principal_graph", cores=5)
  #挑选top10画图展示
  Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
    pull(gene_short_name) %>% as.character()
  #基因表达趋势图

  plot_genes_in_pseudotime(seu.obj.mnc[Track_genes_sig,], color_cells_by="ident",
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


#' Run pesudo analysis
#'
#' @param seu.obj Seurat object
#' @param start.clus NULL
#' @param approx_points 150
#'
#' @return slingshot object with multiple function variable
#' @export
#'
#' @examples
#' # Pls run runSlingshot first
#' # Next run f1, f2 function
#' # Then run runTradeSeq in the backend and save to continue
#' # Final run f4 and f5
runSlingshot <- function(seu.obj, start.clus = NULL, approx_points = 150){

  # 用于构建 sce 的对象
  library(Seurat)
  library(slingshot)

  sce <- as.SingleCellExperiment(seu.obj, assay = "RNA")

  sce_slingshot <- slingshot(sce,      #输入单细胞对象
                              reducedDim = 'UMAP',  #降维方式
                              clusterLabels = sce$ident,  #cell类型
                              reweight = FALSE , # 让细胞不会在多个轨迹间重复评估
                              start.clus = start.clus,     #轨迹起点,也可以不定义
                              approx_points = approx_points)

  res = list(slingshot = sce_slingshot)


  ###### 画图参考 https://www.jianshu.com/p/e85d23a25a43

  plot_trajectory <- function(sce_slingshot = NULL, seu.obj = NULL, color = NULL, lineage.ind = NULL){

    if(is.null(sce_slingshot)){
      stop("sligshot obj is required")
    }
    if(is.null(seu.obj)){
      stop("Seurat obj is required")
    }
    if(is.null(color)){
      stop("Pls set named color vector")
    }

    #### 如果lineage ind为NULL，则画所有的
    if(is.null(lineage.ind)){
        cell_colors <- color[sce_slingshot$ident]
        plot(reducedDims(sce_slingshot)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.8)
        lines(SlingshotDataSet(sce_slingshot), lwd=2, col='black')

        ###################### 1 轨迹图
        #计算celltype坐标位置，用于图中标记
        celltype_label <- seu.obj@reductions$umap@cell.embeddings %>%
          as.data.frame() %>%
          cbind(celltype = Idents(seu.obj) ) %>%
          group_by(celltype) %>%
          summarise(UMAP1 = median(umap_1),
                    UMAP2 = median(umap_2))

        for (i in 1:length(unique(sce_slingshot$ident))) {
          text(celltype_label$celltype[i], x=celltype_label$UMAP1[i]-1, y=celltype_label$UMAP2[i])
        }
    }else{
        plot(reducedDims(sce_slingshot)$UMAP, asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
             col = hcl.colors(100, alpha = 0.5)[cut(slingPseudotime(sce_slingshot)[,lineage.ind], breaks = 100)])

        lines(SlingshotDataSet(sce_slingshot), linInd = lineage.ind, lwd = 2, col = 'black')
    }

  }
  res$f1_plot_trajectory = plot_trajectory


  ####################### 2 密度图
  plot_density = function(sce_slingshot = NULL, color = NULL, lineage.id = NULL){

    # 查看多少lineage SlingshotDataSet(sce_slingshot)

    if(is.null(sce_slingshot)){
      stop("sligshot obj is required")
    }
    if(is.null(color)){
      stop("Pls set named color vector")
    }
    if(is.null(lineage.id)){
      stop("Pls set lineage ID")
    }
    density.list = lapply(unique(sce_slingshot$ident), function(x){
      density(slingPseudotime(sce_slingshot)[colData(sce_slingshot)$ident == x, lineage.id], na.rm = T)
    })
    names(density.list) = unique(sce_slingshot$ident)

    #作图范围
    xlim <- range( unlist( lapply(density.list, function(x) x$x) ) )
    ylim <- range( unlist( lapply(density.list, function(x) x$y) ))

    par(mar=c(6, 6, 6, 6), xpd = TRUE)
    plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "")  #基本作图
    #添加密度图
    for(i in names(density.list)) {
      polygon(c(min(density.list[[i]]$x),density.list[[i]]$x, max(density.list[[i]]$x)), c(0, density.list[[i]]$y, 0),
              col = color[i], alpha = 0.5)
    }

    # legend("right",
    #        inset= -0.2,#legend位于画框右侧正中
    #        pch=15,
    #        legend= names(density.list),
    #        bty="n",
    #        col = color,
    #        border="black",
    #        title = "celltype",
    #        cex = 0.5)

  }
  res$f2_plot_density = plot_density


  # # # ## #########运行 tradeseq时间长，可以单独跑然后保存
  # 微信文章：【单细胞高级分析. 26】slingshot：赏心悦目的轨迹分析图
  # https://blog.csdn.net/qq_43022495/article/details/132708255
  runTradeSeq <- function(sce_slingshot = NULL, variable.features = NULL, nknots = 6){

    # 整体比较 https://www.jianshu.com/p/e85d23a25a43
    # 谱系内比较 https://blog.csdn.net/qq_43022495/article/details/132708255
    if(is.null(sce_slingshot)){
      stop("sligshot obj is required")
    }
    if(is.null(variable.features)){
      stop("Pls set variable features by VariableFeatures(seu.obj)")
    }
    library(tradeSeq)

    # Fit negative binomial model
    counts <- sce_slingshot@assays@data$counts[variable.features,]
    crv <- SlingshotDataSet(sce_slingshot)
    pseudotime <- slingPseudotime(sce_slingshot, na = FALSE)
    cellWeights <- slingCurveWeights(sce_slingshot)


    # 在使用NB-GAM模型之前，需要确定用于构建模型的基函数的数量。
    # 这些基函数被称为节点（knots），它们在模型的拟合过程中起到关键作用。
    # 使用evaluateK函数。该函数会运行一段时间

    # 2k cells ~16min(如果你的轨迹较多的话，时间更长)
    # 运行之后会自动出现图
    # icMat <- evaluateK(counts = counts,
    #                    sds = crv, nGenes = 500,
    #                    k = 3:10,verbose = T, plot = TRUE)
    # https://statomics.github.io/tradeSeq/articles/fitGAM.html
    # 拟合

    library(BiocParallel)
    set.seed(111)
    tradeseq.sce <- fitGAM(counts = counts,
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = nknots, verbose = TRUE,
                  BPPARAM = MulticoreParam(20),parallel=T)

    table(rowData(tradeseq.sce)$tradeSeq$converged)#查看收敛基因个数

    # 整体比较
    assocRes <- associationTest(tradeseq.sce, lineages = TRUE, l2fc = log2(2))
    rowData(tradeseq.sce)$assocRes <- assocRes[order(assocRes$waldStat, decreasing = TRUE),]

    # https://blog.csdn.net/qq_43022495/article/details/132708255
    # 整体比较 https://www.jianshu.com/p/e85d23a25a43
    # 谱系内比较 https://blog.csdn.net/qq_43022495/article/details/132708255
    startRes <- startVsEndTest(tradeseq.sce, lineages=T)
    rowData(tradeseq.sce)$assocRes.start = startRes[ order(startRes$waldStat, decreasing = TRUE), ]

    #### 谱系间比较
    endRes <- diffEndTest(tradeseq.sce,pairwise=T)
    rowData(tradeseq.sce)$assocRes.end <- endRes[order(endRes$waldStat, decreasing = TRUE),]

    tradeseq.sce
  }
  res$f3_runTradeSeq = runTradeSeq

  ############# 基因表达和伪时间
  plot_gene_pesudo_exp = function(gene = NULL, tradeseq.sce = NULL){
    if(is.null(gene) | is.null(tradeseq.sce)){
      stop("Pls set gene name and trade seq objc")
    }
    library(ggpubr)
    p1 = plotSmoothers(tradeseq.sce, assays(tradeseq.sce)$counts,
                      gene = gene, alpha = 0.6, border = T, lwd = 2)+
      ggtitle(gene)
    p1
  }
  res$f4_plot_gene_pesudo_exp = plot_gene_pesudo_exp


  add_pesudoTime_to_seuratMeat <- function(seurat.obj, sce_slingshot.obj){
    t = slingPseudotime(sce_slingshot) %>% data.frame(check.names = F)
    seurat.obj = Seurat::AddMetaData(
      seurat.obj, t
    )
    seurat.obj
  }
  ###################### 拟时间热图

  plot_pesudo_heatmap <- function(seurat.obj = NULL, lineage.name = NULL, genes = NULL,
                                       n_bins = 100,#拟时需要分割的区间，将相似的拟时区间合并，这类似于我们monocle3中的方式
                                       min_exp = 0.2, cutree_rows = 5)
  {
    if(is.null(seurat.obj)| is.null(lineage.name)){
      stop("Pls set Seurat object and lineage.name")
    }
    if(is.null(gene)){
      stop("Pls set candidates genes after test\n\n
           intersect( rownames( rowData(tradeseq.sce)$assocRes.end)[1:50], rownames( rowData(tradeseq.sce)$assocRes.start)[1:50] )
           \n\n")
    }
    seurat_meta = seurat.obj@meta.data
    seurat_meta = as_tibble(cbind(cell.id = as.character(rownames(seurat_meta)), seurat_meta))
    warning("Not consider the cell belong to other lineages")
    seurat_meta = seurat_meta[!is.na(seurat_meta[[lineage.name]]), ]

    seurat_meta = seurat_meta[order(seurat_meta[[lineage.name]]),]


    pl_cells = as.character(seurat_meta$cell.id)

    #提取表达矩阵,并将cell id的exp排序与前面排序好的cell id一致
    expr_mat = as.matrix( seurat.obj@assays$RNA$data[genes,pl_cells] )

    clust_expr_mat = matrix(nrow = nrow(expr_mat),
                            ncol = n_bins, dimnames = list(rownames(expr_mat), 1:n_bins))

    max_pseudotime = max(seurat_meta[[lineage.name]])
    pseudotime_bin_size = max_pseudotime/n_bins

    pseudotime_cluster_stat = NULL
    seurat_meta$pseudotime_bin = NA_integer_

    for (i in 1 : n_bins){
      # 属于特定bin范围内的细胞
      bin_cells = seurat_meta$cell.id[(seurat_meta[[lineage.name]] > (i-1)*pseudotime_bin_size &
                                         seurat_meta[[lineage.name]] <= i*pseudotime_bin_size)]

      # 对细胞归到某个组
      seurat_meta$pseudotime_bin[seurat_meta$cell.id %in% bin_cells] = i

      #计算基因平均表达量
      if (length(bin_cells)>10){
        m2 = expr_mat[,colnames(expr_mat) %in% bin_cells]
        clust_expr_mat[,i] = apply(m2, 1, mean, na.rm = TRUE)
      }

    }

    #数据缩放一下，为了更好的展现热图，并删除低表达基因
    mm = clust_expr_mat - apply(clust_expr_mat, 1, mean, na.rm = TRUE)
    mm = mm[apply(abs(mm),1, max, na.rm = TRUE)> min_exp,]

    ####### 画图
    max_range = max(range(is.finite(mm)))
    lim = c(-max_range, max_range)

    library(pheatmap)
    heatmap1 = pheatmap(mm, show_rownames=F, cluster_rows = TRUE,
                        cluster_cols = FALSE, show_colnames = FALSE,
                        clustering_distance_rows = "euclidean",
                        clustering_method = "ward.D2",
                        treeheight_row = 50,
                        cutree_rows = cutree_rows,
                        color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "PRGn")))(250),
                        breaks = seq(lim[1]/4, lim[2]/4, length.out = 251),
                        border_color = NA)
    heatmap1

   ########## 获取模块的基因
    tree_module = cutree(heatmap1$tree_row, k=cutree_rows)
    tree_module = tibble(gene = names(tree_module), module = as.character(tree_module))
    tree_module = tree_module[heatmap1$tree_row$order,]
    tree_module
  }
  res$f5_plot_pesudo_heatmap = plot_pesudo_heatmap

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

#' Run cellchat
#'
#' @param seurat.obj
#' @param type
#' @param min.cells minimum 100 for analysis
#' @param cores Default 20
#'
#' @return cellChat.Object
#' @export
#'
#' @examples
run_cellChat = function(seurat.obj, type = "triMean", min.cells = 100, cores = 20){


  if(length(intersect("samples", colnames(seurat.obj@meta.data) ) ) == 0 ){
    stop("Plse set @meta.data$samples")
  }

  library(CellChat)

  ####### step create cellchat object
  # Starting from a Seurat object
  if(class(seurat.obj) == "Seurat"){
    cellchat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
  }else if(class(seurat.obj) == "CellChat"){
    cellchat = seurat.obj
  }

  # https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#part-ii-inference-of-cell-cell-communication-network

  ############ Step Select database
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling).
  rm(CellChatDB)

  # set the used database in the object
  cellchat@DB <- CellChatDB.use

  # This step is necessary even if using the whole database
  cellchat <- subsetData(cellchat)

  ########## step Preprocessing the expression data for cell-cell communication analysis
  # https://developer.aliyun.com/article/1256526
  future::plan("multisession", workers = cores) # do parallel
  # #识别过表达基因
  print("identifyOverExpressedGenes")
  cellchat <- identifyOverExpressedGenes(cellchat)
  # #识别过表达配体受体对
  print("identifyOverExpressedInteractions")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  #> The number of highly variable ligand-receptor pairs used for signaling inference is 692

  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)

  options(future.globals.maxSize = 150000 * 1024^2)  # 30G

  ####### Part II: Inference of cell-cell communication network
  # The key parameter for this analysis is type, the method for computing
  # the average gene expression per cell group. By default type = "triMean",
  # producing fewer but stronger interactions. When setting type = "truncatedMean",
  # a value should be assigned to trim, producing more interactions.
  # Please check above in detail on the method for calculating the average gene expression per cell group.
  print("computeCommunProb")
  cellchat <- computeCommunProb(cellchat, type = type)

  # Users can filter out the cell-cell communication if there are only
  # few cells in certain cell groups. By default, the minimum number
  # of cells required in each cell group for cell-cell communication is 10.
  print("filterCommunication")
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)


  ####  step Infer the cell-cell communication at a signaling pathway level
  # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is
  # stored in the slot ‘net’ and ‘netP’, respectively.
  print("computeCommunProbPathway")
  cellchat <- computeCommunProbPathway(cellchat)

  ######## step Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)

  cellchat
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
  exprMat <- GetAssayData(object = seurat.obj, assay = "RNA", layer = "counts")
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

  cat("Resource: ","https://github.com/aertslab/pySCENIC/tree/master/resources/hs_hgnc_tfs.txt", "\n")
  cat("TF:", "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt", "\n")
  cat("motif–TF: ", "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", "\n")
  cat("Database: ","https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/", "\n")

}

#' AUC heatmap for scenic
#'
#' @param auc.file
#' @param seu.obj
#'
#' @return
#' @export
#'
#' @examples
#' loonR::scenicHeatmap(auc.file.path, seu.obj)
scenicHeatmap = function(auc.file =NULL, seu.obj = NULL, color = NULL){
  if(is.null(auc.file)|is.null(seu.obj)){
    stop("set auc.file and seu.obj")
  }
  if(is.null(color)){
    stop("A named vector should be set")
  }
  # load("/master/zhu_zhong_xu/work/TACE-AH/2024-07-09/2024-07-07-allsamples.noRMcc.rdata-hcc.Rdata")
  # auc.file= "/master/zhu_zhong_xu/work/TACE-AH/scenic/hcc.auc_mtx.csv"
  # seu.obj = hcc.cell
  # library(ggsci)
  # pal_npg(alpha = 0.6)(10)
  # color = c(pal_npg(alpha = 0.6)(10),pal_jama(alpha = 0.6)(7), pal_lancet(alpha = 0.6)(9))
  # color = color[c(2:12, 22, 18, 26, 19:25, 13:17)]

  library(Seurat)

  auc.df = read.csv(auc.file, check.names = F, row.names = 1)
  auc.df = auc.df[rownames(seu.obj@meta.data),]

  seu.obj = Seurat::AddMetaData(
    seu.obj, auc.df, colnames(auc.df)
  )
  co2_vector <- c(as.matrix(auc.df))
  hist(co2_vector)

  loonR::heatmap.with.lgfold.riskpro(t(auc.df), Idents(seu.obj),palette = color,
                                     show.lgfold = F, show.risk.pro = F, show_row_names = F,
                                     cluster_rows = T)


  # Seurat::DotPlot(seu.obj, colnames(auc.df))

  auc.df.melted = loonR::meltDataFrameByGroup(auc.df,Idents(seu.obj))
  auc.df$Cell.Type = Idents(seu.obj)

  library(dplyr)
  auc.df.melted = auc.df.melted %>%
    group_by(Group, Gene) %>%
    summarise(mean= mean(value), max = max(value), min = min(value))

  res = list(matrix = auc.df, melt.df = auc.df.melted)
  res

}

#' Run inferCNV
#'
#' @param seurat.obj
#' @param ref_group_names Reference cell type
#' @param gene_pos gene order file
#' @param save_local_path Path to save rds
#' @param num_threads Default: 40
#' @param HMM when set to True, runs HMM to predict CNV level (default: FALSE)
#' @param denoise If True, turns on denoising according to options below (default: FALSE)
#' @param no_plot
#'
#' @return
#' @export
#'
#' @examples
runInferCNV <- function(seurat.obj, ref_group_names = NULL,HMM = F, denoise = F,
                        gene_pos = "~/ref/hg38/gencode.v44/human_genes_pos.txt",
                        save_local_path = NULL, num_threads = 40, no_plot = FALSE){

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
  anno <- data.frame(Idents(seurat.obj), check.names = F)

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
  # https://github.com/broadinstitute/infercnv/issues/418
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = "./infercnv",
                               cluster_by_groups = T, write_expr_matrix = T,
                               HMM = HMM, useRaster=FALSE,
                               denoise = denoise, no_plot = no_plot,
                               num_threads = num_threads)

  if(!is.null(save_local_path)){
    warning("save data to local file: ", save_local_path)
    saveRDS(infercnv_obj, file = save_local_path)
  }else{
    infercnv_obj
  }

}


#' Cal_CNV score
#'
#' @param file_path infercnv.observations.txt
#' @param save_path save path for cnv score
#'
#' @return
#' @export
#'
#' @examples
#' cal_cnvScore("infercnv.observations.txt",  "cnv.score.rdata")
cal_cnvScore <- function(file_path = NULL, save_path = NULL){

  library(dplyr)
  if(is.null(file_path)|is.null(save_path)){
     stop("Pls set file input and output path")
  }

  data = data.table::fread(file_path, nThread = 20, data.table = F)  %>% tibble::column_to_rownames("V1")
  library(dplyr)
  data <- data %>% as.matrix() %>%
    t() %>%
    scale() %>%
    scales::rescale(to=c(-1, 1)) %>%
    t()

  cnv_score <- as.data.frame(colSums(data * data) )
  rownames(cnv_score) = colnames(data)
  colnames(cnv_score) = "cnv.score"

  save(cnv_score, file = save_path)
}




#' CopyKAT
#'
#' @param seurat.obj
#' @param sam.name
#' @param n.cores
#'
#' @return
#' @export
#'
#' @examples
runCopyKAT <- function(seurat.obj = NULL, sam.name = NULL,n.cores =40){

  if(!require(copykat)){
    devtools::install_github("navinlabcode/copykat")
  }
  library(copykat)
  if(is.null(seurat.obj)|is.null(sam.name)){
    stop("Pls set seurat object and sam name")
  }
  copykat.res <- copykat::copykat(rawmat = seurat.obj@assays$RNA$counts, sam.name = sam.name, n.cores = n.cores)
  copykat.res

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




#' Plot dim plot and cell count bar
#'
#' @param object
#' @param color
#' @param dim
#' @param rm.axis
#' @param cell.id
#' @param bar.width
#' @param point.size
#' @param rm.barplot
#' @param legend.psize
#' @param arrow.len
#'
#' @return
#' @export
#'
#' @examples
#' # modified from scRNAtoolVis::scatterCellPlot
#'
dimPlotWithCountBar <- function (object = NULL, color = NULL, dim = "umap", rm.axis = FALSE,
          cell.id = NULL, bar.width = 0.2, point.size = 1, rm.barplot = FALSE,
          legend.psize = 1.5, arrow.len = 0.2)
{
  library(Seurat)
  library(ggplot2)
  library(grid)

  reduc <- data.frame(Seurat::Embeddings(object, reduction = dim))
  meta <- object@meta.data
  pc12 <- cbind(reduc, meta)
  pc12$idents <- as.character(Seurat::Idents(object))
  if (is.null(cell.id)) {
    cell_num <- dplyr::arrange(dplyr::summarise(dplyr::group_by(pc12,
                                                                idents), n = dplyr::n()), n)
  }
  else {
    cell_num <- dplyr::arrange(dplyr::summarise(dplyr::group_by(pc12,
                                                                idents, .data[[cell.id]]), n = dplyr::n()), n)
  }
  if (rm.axis == FALSE) {
    lab.shift = unit(-2.5, "lines")
  }
  else {
    lab.shift = unit(-1, "lines")
  }
  grid.newpage()
  pushViewport(viewport(x = unit(0.1, "npc"), y = unit(0.5,
                                                       "npc"), width = unit(0.5, "npc"), height = unit(0.7,
                                                                                                       "npc"), just = "left", xscale = grDevices::extendrange(range(pc12[,
                                                                                                                                                                         1]), f = 0.05), yscale = grDevices::extendrange(range(pc12[,
                                                                                                                                                                                                                                    2]), f = 0.05), ))
  grid.rect()
  if (rm.axis == FALSE) {
    jjPlot::grid.xaxis2(label.space = 0.5)
    jjPlot::grid.yaxis2(label.space = 0.25)
  }
  celltype <- cell_num$idents
  if (is.null(color)) {
    cols <- circlize::rand_color(n = length(celltype))
  }
  else {
    cols <- color
  }
  ################# modified by myself 2024-08-01
  if(!is.null(names(cols))){
  cols = cols[celltype]
  }
  # # ## ##### ##### #### ## ###

  for (i in seq_along(celltype)) {
    tmp <- pc12[which(pc12$idents %in% celltype[i]), ]
    grid.points(x = tmp[, 1], y = tmp[, 2], pch = 19, size = unit(point.size,
                                                                  "pt"), gp = gpar(col = cols[i]))
  }
  if (rm.axis == TRUE) {
    grid.segments(x0 = 0.025, x1 = arrow.len, y0 = 0.05,
                  y1 = 0.05, arrow = arrow(length = unit(2, "mm"),
                                           type = "closed"), gp = gpar(fill = "black"))
    grid.text(label = paste0(toupper(dim), " 1"), x = (arrow.len +
                                                         0.025)/2, y = 0.025, gp = gpar(fontsize = 6, fontface = "bold.italic"))
    grid.segments(x0 = 0.05, x1 = 0.05, y0 = 0.025, y1 = arrow.len,
                  arrow = arrow(length = unit(2, "mm"), type = "closed"),
                  gp = gpar(fill = "black"))
    grid.text(label = paste0(toupper(dim), " 2"), x = 0.025,
              y = (arrow.len + 0.025)/2, rot = 90, gp = gpar(fontsize = 6,
                                                             fontface = "bold.italic"))
  }
  else {
    grid.text(label = paste0(toupper(dim), " dimension 1"),
              x = 0.5, y = lab.shift)
    grid.text(label = paste0(toupper(dim), " dimension 2"),
              x = lab.shift, y = 0.5, rot = 90)
  }
  popViewport()
  if (isFALSE(rm.barplot)) {
    pushViewport(viewport(x = unit(0.61, "npc"), y = unit(0.5,
                                                          "npc"), width = unit(bar.width, "npc"), height = unit(0.7,
                                                                                                                "npc"), just = "left", yscale = c(0, nrow(cell_num) +
                                                                                                                                                    0.75), xscale = c(0, max(cell_num$n) + 0.1 * max(cell_num$n))))
    if (rm.axis == FALSE) {
      jjPlot::grid.xaxis2(label.space = 0.5, at = c(0,
                                                    max(cell_num$n)), labels = as.character(c(0,
                                                                                              max(cell_num$n))))
    }
    grid.rect(x = rep(0, nrow(cell_num)), y = 1:nrow(cell_num),
              width = cell_num$n, height = unit(0.08, "npc"),
              just = "left", gp = gpar(fill = cols, col = NA),
              default.units = "native")
    grid.rect(gp = gpar(fill = "transparent"))
    grid.text(label = "Number of cells", x = 0.5, y = lab.shift)
    popViewport()
  }
  if (isTRUE(rm.barplot)) {
    bar.width = 0
  }
  pushViewport(viewport(x = unit(0.61 + bar.width, "npc"),
                        y = unit(0.5, "npc"), width = unit(0.2, "npc"), height = unit(0.7,
                                                                                      "npc"), just = "left", yscale = c(0, nrow(cell_num) +
                                                                                                                          0.75)))
  grid.points(x = rep(0.1, nrow(cell_num)), y = 1:nrow(cell_num),
              pch = 19, gp = gpar(col = cols), size = unit(legend.psize,
                                                           "char"))
  if (!is.null(cell.id)) {
    grid.text(label = as.character(unlist(cell_num[, cell.id])),
              x = 0.1, y = 1:nrow(cell_num), default.units = "native")
  }
  grid.text(label = cell_num$idents, x = 0.2, y = 1:nrow(cell_num),
            just = "left", gp = gpar(fontsize = 10), default.units = "native")
  popViewport()
}



#' Correlation plot based on average expression
#'
#' @param seurat.obj
#' @param group.by
#' @param features
#'
#' @return correlation plot
#' @export
#'
#' @examples
groupCorrelation <- function(seurat.obj, group.by = "ident",
                             features = NULL){

  av.expr = Seurat::AggregateExpression(
      seurat.obj,
      assays = "RNA",
      features = features,
      return.seurat = FALSE,
      group.by = group.by)

  av.expr = data.frame(av.expr, check.names = F)
  colnames(av.expr) = stringr::str_remove(colnames(av.expr),"^RNA.g")

  cg = names(head( sort( apply(av.expr, 1, mad), decreasing = T ), 1000))


  M = cor(av.expr[cg,])


  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # matrix of the p-value of the correlation
  p.mat <- cor.mtest(av.expr[cg,])


  color = colorRampPalette(colors=c("black", "gray", "#FFFFFF", "steelblue", "red"))(20)

  corrplot::corrplot(M, method="circle", type="lower",
                     order="hclust", col= color, cl.pos = 'n',
                     # bg="lightblue",
                     col.lim = c(0, 1) )



}

#' Run metaflux pipeline
#'
#' @param seurat.obj
#' @param myident
#' @param n_bootstrap
#' @param seed
#' @param medium
#'
#' @return
#' @export
#'
#' @examples
runMetaFlux <- function(seurat.obj = NULL, myident = NULL, n_bootstrap = 100, seed=1, medium = NULL){
  library(METAFlux)

  if(is.null(medium)){
    print("Medium not set, will use medium = cell_medium")
    print("Else you may set human_blood")
    medium = cell_medium
  }else if(medium == "cell_medium"){
    medium = cell_medium
  }else if(medium == "human_blood"){
    medium = human_blood
  }else{
    stop("Pls set medium correctlly: human_blood or cell medium")
  }

  if(is.null(seurat.obj)){
    stop("Pls input seurat object")
  }
  if(is.null(myident)){
    stop("Pls set idents col name in metadata")
  }

  # https://htmlpreview.github.io/?https://github.com/KChen-lab/METAFlux/blob/main/Tutorials/pipeline.html
  generate_boots <- function(celltype, n) {
    dt <- data.frame(cluster = celltype, id = 1:length(celltype))
    index <- do.call(cbind, sapply(1:n, function(x) {
      splits <- dt %>%
        group_by(cluster) %>%
        sample_n(dplyr::n(), replace = TRUE) %>%
        ungroup() %>%
        dplyr::select("id")
    }))
    return(index)
  }

  get_ave_exp <- function(i, myseurat, samples,myident) {
    meta.data=myseurat@meta.data[samples[,i],]
    sample <-GetAssayData(myseurat, assay = "RNA")[,samples[,i]]
    name <- colnames(sample)
    for (j in 1:length(name)) {
      name[j] <- paste0(name[j], "_", j)
    }
    colnames(sample) <- name
    rownames(meta.data) <- name
    SeuratObject<-suppressWarnings(
      CreateSeuratObject(count=sample, meta.data = meta.data))
    SeuratObject<-NormalizeData(SeuratObject,verbose = TRUE)
    ave<-GetAssayData(AverageExpression(SeuratObject,group.by = myident,return.seurat = T), assay = "RNA") %>% as.data.frame()
    return(ave)
  }



  #Calculate the mean expression for bootstrap samples from seurat object
  #Using "Cell_type" in seurat as my grouping variable
  #For the sake of demonstration, we only set the number of bootstraps to 3.
  #In real analysis, the number of bootstraps should be much higher(e.g. 100, 200....)
  ###### https://github.com/KChen-lab/METAFlux/issues/11
  edit_calculate_avg_exp <- function(myseurat,myident,n_bootstrap,seed) {
    set.seed(seed)
    samples=generate_boots(myseurat@meta.data[,myident],n_bootstrap)
    exp <- lapply(1:n_bootstrap,get_ave_exp,myseurat,samples,myident)
    exp <- do.call(cbind, exp)
    return(exp)
  }


  mean_exp = edit_calculate_avg_exp(myseurat = seurat.obj, myident = myident, n_bootstrap = n_bootstrap, seed = seed)

  # Calculate MRAS(Metabolic Reaction Activity Score)
  # We can calculate single sample normalized MRAS from data using Gene-protein-reaction (GPR).This step is the same with the bulk RNA-seq pipeline.
  scores <- calculate_reaction_score(data=mean_exp)

  #calculate the fractions of celltype/clusters

  fra.tb = round(table(seurat.obj[[myident]])/nrow(seurat.obj@meta.data),1)
  fra.tb = data.frame(fra.tb, check.names = F)
  fra.vec = fra.tb$Freq
  names(fra.vec) = fra.tb$groups
  rm(fra.tb)
  ## re-order
  fra.vec = fra.vec[colnames(mean_exp)[1:length(unique(unlist(seurat.obj[[myident]])))]]

  #num_cell=number of cell types/clusters, here we have 4 cell types, thus input is 4. Make sure the sum of cell percentage is 1.The order of fraction should be the same as that of "Mean_exp" data
  flux = compute_sc_flux(num_cell = length(fra.vec), fraction = fra.vec, fluxscore = scores, medium = medium)

  #optional: flux scores normalization
  cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
  }

  flux=cbrt(flux)

  data("human_gem")
  data("nutrient_lookup_files")

  res = list(MetabolicReactionActivityScore = scores, flux = flux, medium = medium, human_gem = human_gem, nutrient_lookup_files = nutrient_lookup_files)

  res

}



#' run spotlight
#'
#' @param scRNA Seurat object
#' @param spatial Seurat object
#' @param downsample Default 100
#' @param feature.auc Default 0.75 when selecting cluster marker
#' @param top_n_feature Default 50
#'
#' @return
#' @export
#'
#' @examples
run_spotlight = function(scRNA=NULL, spatial=NULL, downsample=100, feature.auc = 0.75, top_n_feature = 50){

  if(is.null(scRNA)|is.null(spatial)){
    stop("provide scRNA and spatial RNA")
  }

  library(scater)
  library(scran)
  library(dplyr)
  library(Seurat)

  library(SPOTlight)
  cat("\nConvert to sce object\n")

  scRNA$spot.cluster= Idents(scRNA)

  sce <- as.SingleCellExperiment(scRNA)

  cat("\nFeature selection\n")
  ### Feature selection
  sce <- logNormCounts(sce)

  # Get vector indicating which genes are neither ribosomal or mitochondrial
  genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

  dec <- modelGeneVar(sce, subset.row = genes)
  plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
  curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

  cat("\nTop genes\n")
  # Get the top 3000 genes.
  hvg <- getTopHVGs(dec, n = 3000)

  # 加上细胞注释信息
  colLabels(sce) <- colData(sce)$spot.cluster

  # Compute marker genes
  mgs <- scoreMarkers(sce, subset.row = genes)

  # 保留最相关的marker基因
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > feature.auc, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)

  mgs_df <- mgs_df %>% group_by(cluster) %>% top_n(top_n_feature, mean.AUC) %>% ungroup() %>% as.data.frame()

  cat("\nMarkers \n")
  cat(table(mgs_df$cluster))

  # split cell indices by identity
  idx <- split(seq(ncol(sce)), sce$spot.cluster)
  # downsample to at most 20 per identity & subset
  # We are using 5 here to speed up the process but set to 75-100 for your real
  # life analysis
  n_cells <- downsample
  cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    set.seed(111)
    sample(i, n_cells)
  })
  sce <- sce[, unlist(cs_keep)]

  cat("\nRun spotlight\n")
  res <- SPOTlight(
    x = sce,
    y = spatial@assays$Spatial$counts,
    groups = as.character(sce$spot.cluster),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

  # 提取比例矩阵
  deconv_mat <- res$mat
  spatial@meta.data <- cbind(spatial@meta.data, deconv_mat)
  spatial
}


#' run_RCTD
#'
#' @param scRNA
#' @param spatial
#' @param rctd_mode RCTD has three modes, each suited for different spatial technologies: "doublet": Assigns 1-2 cell types per pixel. Best for high-resolution technologies like Slide-seq and MERFISH. "multi": Like doublet mode but can assign more cell types (up to max_multi_types). Best for lower-resolution technologies like 100-micron Visium. "full": No restrictions on number of cell types.
#'
#' @return
#' @export
#'
#' @examples
run_RCTD = function(scRNA=NULL, spatial=NULL, rctd_mode = "doublet", max_multi_types = 10){

  if(is.null(scRNA)|is.null(spatial)){
    stop("provide scRNA and spatial RNA")
  }
  library(SpatialExperiment)
  library(SummarizedExperiment)

  library(spacexr)
  library(Seurat)
  # Create SummarizedExperiment.

  reference_se <- SummarizedExperiment(
    assays = list(counts = scRNA@assays$RNA$counts),
    colData = data.frame(cell_type = Idents(scRNA),
                         row.names = rownames(scRNA@meta.data) )
  )


  # 获取坐标信息
  coords <- GetTissueCoordinates(spatial)
  # 为 RCTD 准备坐标
  spatial_coords <- data.frame(
    row.names =  rownames(coords),
    x = coords$imagecol,
    y = coords$imagerow
  )

  # Create SpatialExperiment.
  spatial_spe <- SpatialExperiment(
      assay = spatial@assays$Spatial$counts,
      spatialCoords = as.matrix(spatial_coords)
  )

  # Preprocess data
  rctd_data <- createRctd(spatial_spe, reference_se)
  results_spe <- runRctd(rctd_data, rctd_mode = rctd_mode, max_cores = 50, max_multi_types = 10)
  resutts_weight = data.frame(assays(results_spe)[["weights"]], check.names = F)
  resutts_weight = t(resutts_weight)

  spatial@meta.data <- cbind(spatial@meta.data,
                             data.frame(resutts_weight, check.names=F)[rownames(spatial@meta.data),])
  spatial

}


