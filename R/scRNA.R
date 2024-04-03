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
#' @param harmony_batch Column names from batch correction
#'
#' @return
#' @export
#'
#' @examples
load10X <- function(paths = NA, sample_names = NA, min.cells = 10,
          min.features = 200, gene.column = 1, remove_doublet = T, batch = NULL,
          max_nCount_RNA = 10000, max_nFeature_RNA = 8000, max_percent.mt = 20,
          integrate_CCA = FALSE, top_variable_features = 2000, remove_cc = T, remove_mt = F,
          integrate_harmony = F, harmony_batch = "batch"){

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

    cat("\n--------------------------\nload", sample_names[i],"\n")

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

    data.filt <- RunPCA(data.filt) %>% RunUMAP(dims = 1:30) %>%  RunTSNE(dims = 1:30)


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
    cat("\n--------------------------\nPerform harmony integration\n")

    filter.data = merge(Object_list_filter[[1]], Object_list_filter[2:length(Object_list_filter)])

    filter.data <- SCTransform(filter.data, vars.to.regress = vars.to.regress, verbose = FALSE)

    last.assay.name = "SCT"

    filter.data <- RunPCA(filter.data, assay = last.assay.name)

    #filter.data <- RunPCA(filter.data, features = VariableFeatures(filter.data) )

    filter.data <- harmony::RunHarmony(filter.data, vars_use = harmony_batch,
                                       reduction = "pca", assay.use = "SCT",
                                       reduction.save = "harmony")
    last.assay.name = "SCT"
    reduction = "harmony"

  }else{
    cat("\n--------------------------\nMerge samples\n")

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


  filter.data <- FindNeighbors(filter.data, dims = 1:30, assay = last.assay.name )

  cat("\n--------------------------\nFind Cluster\n")

  #计算SNN
  filter.data <- FindClusters(
    object = filter.data,
    resolution = c(seq(.1,1.6,.2))
  )

  ##### umap
  cat("\n--------------------------\nRun UMAP \n")
  filter.data <- RunUMAP(filter.data, dims = 1:30, reduction = reduction, assay = last.assay.name )
  ##### tsne
  cat("\n--------------------------\nRun TSNE \n")
  filter.data <- RunTSNE(filter.data, dims = 1:30, reduction = reduction, assay = last.assay.name)

  DefaultAssay(filter.data) = "RNA"

  # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
  res = list(raw = raw.data, filter = filter.data)

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
