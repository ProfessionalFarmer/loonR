

#' Title
#'
#' @param affinityL List including different affinity matrix
#' @param evidence.type Corresponded to affinity list
#' @param group SNF subtyping result
#'
#' @return
#' @export
#'
#' @examples
#' loonR::SNF_Similairity_Hist(tmp$AffinityL,evidence.type = c("RNA","Methylation", "CNV"), group = snf.group)
SNF_Similairity_Hist <- function(affinityL=NULL, evidence.type=NULL, group = NULL){
  library(reshape2)

  affinityL.matrix <- affinityL
  # evidence.type = c("mRNA","Metylation","CNV")
  # group = snf.group

  if (is.null(affinityL)) {
    message("No 'affinityL' value defined")
    stop()
  }

  if (is.null(evidence.type)) {
    message("No 'evidence.type' value defined.")
    stop()
  }

  if (is.null(group)) {
    message("No 'group' value defined. Must supply with names")
    stop()
  }


  # 1. Melt get paired similarity
  names(affinityL.matrix) <- evidence.type
  affinityL.matrix.melt <- lapply(affinityL.matrix, function(m){
    m[upper.tri(m, diag = TRUE)] <- NA
    m = melt(m)
    m = na.omit(m)
    m$pair <- paste(m[,c(1)],m[,c(2)],sep="-")
    m <- m[,-c(1,2)]
    m
  })
  affinityL.matrix.melt <- lapply(evidence.type, function(x){
    colnames( affinityL.matrix.melt[[x]] ) <- c(x, "pair")
    affinityL.matrix.melt[[x]]
  })
  names( affinityL.matrix.melt ) <- evidence.type


  # 2. Histogram plot
  library(ggpubr)
  affinityL.similarity.hist <- lapply(names(affinityL.matrix.melt), function(x){

    affinityL.matrix.melt[[x]] <- affinityL.matrix.melt[[x]]

    p <- gghistogram(affinityL.matrix.melt[[x]], x=x, xlab = "Similarity", title = x, ylab = "# of pairs of patients")
    p

  })

  names( affinityL.similarity.hist ) <- names(affinityL.matrix.melt)


  # 3. Evidence support count
  # not related with step 1 and 2
  affinityL.evidence.support <- lapply(unique(group), function(x){
    smps = names(group[group==x])

    # similar to step 1 and 2
    single.group <- lapply(evidence.type, function(x){
      m = affinityL.matrix[[x]]
      m = m[smps,smps]
      m = melt(m)
      m = na.omit(m)
      m$pair <- paste(m[,c(1)],m[,c(2)],sep="-")
      m <- m[,-c(1,2)]
      colnames(m) <- c(x,"pair")
      m
    })

    single.group <- Reduce(function(x,y) merge(x=x, y=y, by = "pair"), single.group)
    single.group <- single.group[,-c(1)]

    # get supported evidence
    single.group.support.evidence <- apply(single.group, MARGIN=1, FUN = function(value){

      value.rank <- order(value, decreasing = T)
      support.evidence.vector <- names(value)[value.rank[1]]

      for(ind in 2:length(value.rank)){
        if ( 1.1 *  value[ value.rank[ind] ] > value[ value.rank[1] ] ){
          support.evidence.vector <- c(support.evidence.vector, names(value)[value.rank[ind]] )
        } else{
          break
        }
      }
      support.evidence.vector = paste(  sort(support.evidence.vector), collapse = "-")
      support.evidence.vector
    })


    single.group$support.evidence = single.group.support.evidence
    single.group
  })
  names( affinityL.evidence.support ) <- paste("Group",unique(group),sep=" ")


  # step 4 support evidence pie plot
  affinityL.evidence.support.pie <- lapply( names(affinityL.evidence.support), function(x){
    tmp.data <- affinityL.evidence.support[[x]]
    p = loonR::plotPie(tmp.data$support.evidence, title = x, color="Most")
    p

  } )
  names(affinityL.evidence.support.pie) <-  names(affinityL.evidence.support)

  #Plot given similarity matrix by clusters
  #displayClusters(W, snf.group)


  # Ranks each features by NMI based on their clustering assingments
  #NMI_scores <- rankFeaturesByNMI(dataL, W)

  # ind = order(snf.group)
  # consensusmap(W,color=c(rep("#191970",2),rep("white",40)),
  #              Rowv=ind,
  #              Colv=ind,
  #              main = "",  ## 1. #191970 #00FFFF
  #              annCol = annotation,
  #              annColors = ann_colors,
  #              scale = "none",legend = F,
  #              annLegend = F, cellwidth = 0.1, cellheight = 0.1)

  res<- list(paired.similarity = Reduce(function(x,y)
                                        merge(x=x, y=y, by = "pair"),
                                        affinityL.matrix.melt), # merged data.frame
    paired.similarity.hist = affinityL.similarity.hist,

    group.support.evidence = affinityL.evidence.support,
    group.support.evidence.pie = affinityL.evidence.support.pie

  )
  res
}

#' Title
#'
#' @param x
#' @param similarity
#'
#' @return
#' @export
#'
#' @examples
tanimoto <- function(x, similarity=F) {
  res<-sapply(x, function(x1){
    sapply(x, function(x2) {i=length(which(x1 & x2)) / length(which(x1 | x2)); ifelse(is.na(i), 0, i)})
  })
  if(similarity==T) return(res)
  else return(1-res)
}


#' Title Similarity network fusion
#'
#' @param dataL list( t(mRNA.snf.df), t(methylation.snf.df), t(cnv.snf.df) )
#' @param alpha Default 0.5. hyperparameter, usually (0.3~0.8)   Variance for local model
#' @param K Default 20. Number of neighbors, must be greater than 1. usually (10~30)  20
#' @param Iterations T.Default 20. Number of Iterations, usually (10~50)
#' @param dist.method Euclidean, Pearson, Spearman, mutualInfo, Manhattan, KLD
#' @param survival Must a data frame. colnames OS.time, OS.event, RFS.time, RFS.event. Rownames must be sample name. Two columns or four columns.
#'
#' @return
#' @export
#'
#' @examples run_SNF( list( t(mRNA.snf.df), t(methylation.snf.df), t(cnv.snf.df) ),  alpha = 0.5, K = 20, Iterations = 20    )
#' https://cran.r-project.org/web/packages/SNFtool/SNFtool.pdf
#' Distance reference: https://www.rdocumentation.org/packages/bioDist/versions/1.44.0
run_SNF <- function(dataL = NULL, alpha = 0.5, K = 20, Iterations = 20, dist.method ="Euclidean", survival=NA, max.cluster = 5, std.normalize = TRUE){

  library(SNFtool)
  library(bioDist)

  T = Iterations

  if (is.null(dataL)) {
    message("No 'dataL' value defined")
    stop()
  }

  ################# Normalization
  if(std.normalize){
    # normalize Normalize each column of the input data to have mean 0 and standard deviation 1.
    dataL.normalized = lapply(dataL, SNFtool::standardNormalization) #注意前面是否normalization了，注意这里是否需要log2
    # min max
  }else{
    dataL.normalized = dataL
  }


  ################ Calculate distance
  # calculate disctance   'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
  # chiDist2  Pairwise Chi-squared distances  dist 2 Euclidean Distance   dist2(x, x)^(1/2)
  # dist(data,method="euclidean") method "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  # library(bioDist)  * KLD.matrix(Continuous version of Kullback-Leibler Distance (KLD)),*  cor.dist(Pearson correlational distance),* spearman.dist(Spearman correlational distance),* tau.dist(Kendall's tau correlational distance),* man (Manhattan distance),* KLdist.matriX(Discrete version of Kullback-Leibler Distance (KLD)), *  mutualInfo(Mutual Information),*  euc(Euclidean distance)   https://www.rdocumentation.org/packages/bioDist/versions/1.44.0

  if(dist.method == "Euclidean"){
    dataL.dist = lapply(dataL.normalized, function(x) (SNFtool::dist2(x, x)^(1/2)))   # Euclidean Distance
  }else if(dist.method == "Pearson"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::cor.dist(x) )   # Pearson correlational distance
  }else if(dist.method == "Spearman"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::spearman.dist(x) ) # Spearman correlational distance
  }else if(dist.method == "mutualInfo"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::mutualInfo(x) ) # mutualInfo correlational distance
  }else if(dist.method == "Kendall"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::tau.dist(x) ) # Kendall's tau correlational distance
  }else if(dist.method == "Manhattan"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::man(x) ) # Manhattan distance
  }else if(dist.method == "KLD"){
    dataL.dist = lapply(dataL.normalized, function(x) bioDist::KLdist.matriX(x) ) # Discrete version of Kullback-Leibler Distance (KLD)
  }else if(dist.method %in% c("tanimoto","jaccard","euclidean","hamming","cont tanimoto", "MCA_coord","gower","chi.squared","cosine") ){
    dataL.dist = lapply(dataL.normalized, function(x) IntClustd::KLdist.matriX(x,distmeasure=dist.method,normalize=FALSE) )
  }else{
    stop("Distance: Euclidean, Pearson, Spearman, mutualInfo, Manhattan, KLD")
  }


  # Construct the similarity graphs
  affinityL = lapply(dataL.dist, function(x) affinityMatrix(as.matrix(x), K, alpha))
  W = SNF(affinityL, K, T)

  ## an example of how to use concordanceNetworkNMI
  ## The output, Concordance_matrix,
  ## shows the concordance between the fused network and each individual network.
  # Concordance_matrix = concordanceNetworkNMI(affinityL, 3)


  # Estimate Number Of Clusters Given Graph
  estimationResult = estimateNumberOfClustersGivenGraph(W, NUMC=2:max.cluster)
  #cat(estimationResult)

  snf.group <- lapply(2:max.cluster, function(i){
    once.group <- spectralClustering(W,i)
    names(once.group) <- colnames(W)
    once.group
  })
  names(snf.group) <- paste("Cluster",2:max.cluster,sep="")


  #Plot given similarity matrix by clusters
  #displayClusters(W, snf.group)
  p.survival = NA
  if (!is.null(survival)) {
    # 1对应的是2个分组
    tmp.surv.df <- survival[names(snf.group[[1]]),]
    best1 = estimationResult$`Eigen-gap best`
    best2 = estimationResult$`Eigen-gap 2nd best`

    surv.col.num = ncol(survival)

    library(foreach)
    p.survival <- foreach::foreach(i = 1:length(snf.group), .combine = rbind) %do% {
      tmp.surv.df$SNF <- as.numeric( snf.group[[i]] )
      diff = survdiff( Surv(OS.time, OS.event) ~ SNF , data = tmp.surv.df)
      os.pval <- round(broom::glance(diff)$p.value,3)
      if (surv.col.num > 3){
        diff = survdiff( Surv(RFS.time, RFS.event) ~ SNF , data = tmp.surv.df)
        rfs.pval <- round(broom::glance(diff)$p.value,3)
      }
      if(surv.col.num > 3){
        fe.res <- c(i+1, os.pval, rfs.pval, alpha, K, T, best1, best2)
        names(fe.res) <- c("Group", "OS.P","RFS.P","alpha","K","T", "Best1", "Best2")
      }else{
        fe.res <- c(i+1, os.pval, alpha, K, T, best1, best2)
        names(fe.res) <- c("Group", "OS.P","alpha","K","T", "Best1", "Best2")
      }
      fe.res
    }


  }



  res = list(AffinityL = affinityL, Wall = W, alpha = alpha, K = K, Iterations = Iterations,
             Data.normalized = dataL.normalized,
             Data.dist = dataL.dist,
             EstimateResult = estimationResult,
             Clustering = snf.group,
             Survival.Analysis = p.survival)

}




