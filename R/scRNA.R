#' Read data: raw to FindVariableFeatures
#'
#' @param min.cells 3
#' @param min.features 200
#' @param gene.column 1
#' @param paths
#' @param sample_names
#' @param max_nCount_RNA
#' @param max_nFeature_RNA 8000
#' @param max_percent.mt 20
#' @param integrate_CCA Default FALSE
#' @param top_variable_features 2000
#'
#' @return
#' @export
#'
#' @examples
load10X <- function(paths = NA, sample_names = NA, min.cells = 3,
          min.features = 200, gene.column = 1, remove_doublet = T,
          max_nCount_RNA = 10000, max_nFeature_RNA = 8000, max_percent.mt = 20,
          integrate_CCA = FALSE, top_variable_features = 2000, remove_cc = T){

  Object_list_raw = list()
  Object_list_filter = list()

  if(sum(is.na(paths)|is.na(sample_names))!=0){
    stop("Should provide file path and sample name")
  }
  if(length(paths) != length(sample_names)){
    stop("Should equal length")
  }

  for (i in 1:length(paths)){

    cat("\n--------------------------\nlaod", sample_names[i],"\n")

    data <- Read10X(paths[i], gene.column = gene.column)

    data <- CreateSeuratObject(data, min.cells = min.cells, min.features = min.features)

    data[["sample"]] = sample_names[i]
    data[["orig.ident"]] = sample_names[i]

    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")

    ############### filtering
    data.filt <- subset(data,
                        subset= nCount_RNA > min.features & nCount_RNA < max_nCount_RNA &
                          nFeature_RNA > min.features & nFeature_RNA < max_nFeature_RNA &
                          percent.mt < max_percent.mt )

    #### calculate cell cycle
    data.filt <- CellCycleScoring(
      object = data.filt,
      g2m.features = cc.genes$g2m.genes,
      s.features = cc.genes$s.genes
    )

    if(remove_doublet){
      data.filt = NormalizeData(data.filt)
      data.filt = FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = top_variable_features)
      data.filt = ScaleData(data.filt)
      data.filt <- RunPCA(data.filt, dims = 30)
      data.filt <- FindNeighbors(data.filt, dims = 30)

      data.filt <- FindClusters(
        object = data.filt,
        resolution = c(0.1,0.3)
      )

      data.filt = loonR::Find_doublet(data.filt)
      # 提取判定为单胞的细胞进行下游分析
      data.filt <- subset(data.filt, subset=doublet_info=="Singlet")
    }

    Object_list_raw[[sample_names[i]]] = data

    Object_list_filter[[sample_names[i]]] = data.filt

  }
  raw.data = merge(Object_list_raw[[1]], Object_list_raw[2:length(Object_list_raw)])


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

  }else{
    cat("\n--------------------------\nMerge samples\n")
    filter.data = merge(Object_list_filter[[1]], Object_list_filter[2:length(Object_list_filter)])
    filter.data = NormalizeData(filter.data)
    filter.data = FindVariableFeatures(filter.data,
                                       selection.method = "vst",
                                       nfeatures = top_variable_features)
  }

  rm(Object_list_filter, Object_list_raw)

  filter.data <- JoinLayers(filter.data)

  ########### 重新跑normalize  (前面已跑) findvariable  (重跑)scale （待跑） 的流程

  if(remove_cc){
    filter.data <- ScaleData(filter.data,
                             vars.to.regress = c("S.Score","G2M.Score"),
                             features = VariableFeatures(filter.data),
                             assay = "integrated"
    )
  }else{
    filter.data <- ScaleData(filter.data,
                             features = VariableFeatures(filter.data),
                             assay = "integrated"
    )
  }

  ########### 然后是PCA
  filter.data <- RunPCA(filter.data, dims = 30, assay = "integrated")
  filter.data <- FindNeighbors(filter.data, dims = 1:30, assay = "integrated" )


  DefaultAssay(filter.data) = "integrated"
  #计算SNN
  filter.data <- FindClusters(
    object = filter.data,
    resolution = c(seq(.1,1.6,.2))
  )
  DefaultAssay(filter.data) = "RNA"

  ##### umap
  filter.data <- RunUMAP(filter.data, dims = 1:30, reduction.use = "pca", assay = "integrated" )
  filter.data <- RunTSNE(filter.data, dims = 1:30, reduction.use = "pca", assay = "integrated")


  # https://cran.r-project.org/web/packages/harmony/vignettes/Seurat.html
  res = list(raw = raw.data, filter = filter.data)

  res
}


#' simple QC figures
#'
#' @param seurat_obj
#'
#' @return
#' @export
#'
#' @examples
prelimenaryQC <- function(seurat_obj){

  p1 = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0.1)

  p2 = FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

  p3 = FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")

  p4 = FeatureScatter(seurat_obj, feature1="percent.ribo", feature2="nFeature_RNA")

  res = list(Volin = p1, Count_RNA = p2, Feature_mt = p3, ribo_Feature = p4)

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
  sweep.res.list <- paramSweep(data, PCs = PCs, sct = sct) # 若使用SCT方法 标准化则'sct=T'
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

