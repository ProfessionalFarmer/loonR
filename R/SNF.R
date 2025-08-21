

#' Similairity count
#'
#' @param affinityL List including different affinity matrix
#' @param evidence.type Corresponded to affinity list
#' @param group SNF subtyping result
#' @param group.prefix Default Group. User can specify
#'
#' @return
#' @export
#'
#' @examples
#' loonR::SNF_Similairity_Hist(tmp$AffinityL,evidence.type = c("RNA","Methylation", "CNV"), group = snf.group)
SNF_Similairity_Hist <- function(affinityL=NULL, evidence.type=NULL, group = NULL, group.prefix = "Group"){

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

      m[upper.tri(m, diag = TRUE)] <- NA

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
        if ( value[ value.rank[ind] ] > 0.9 * value[ value.rank[1] ] ){
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
  names( affinityL.evidence.support ) <- paste(group.prefix, unique(group),sep=" ")


  # obtain all the evidence combinations
  library(foreach)
  combs.res <- foreach(c = 1:length(evidence.type), .combine = c)%do%{

    combs <- gtools::combinations(length(evidence.type), c, evidence.type)

    if(c==1){
      combs = unclass(unlist(combs[,1]))
    }else{
      combs = apply(combs, 1, function(x) paste0(sort(x), sep="", collapse = "-"))
    }

    combs
  }

  pie.colors <- loonR::get.mostDistint.color.palette(n=length(combs.res))
  names(pie.colors) <- combs.res


  # step 4 support evidence pie plot
  affinityL.evidence.support.pie <- lapply( names(affinityL.evidence.support), function(x){
    tmp.data <- affinityL.evidence.support[[x]]
    p = loonR::plotPie(tmp.data$support.evidence, title = x,
                       color=pie.colors[sort(unique(tmp.data$support.evidence))]  )
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

#' Calculate tanimoto distance
#'
#' @param x matrix
#' @param similarity Default F. If true return similarity, else return distance
#'
#' @return
#' @export
#'
#' @examples
tanimoto_distance <- function(x, similarity=F) {
  res<-sapply(x, function(x1){
    sapply(x, function(x2) {i=length(which(x1 & x2)) / length(which(x1 | x2)); ifelse(is.na(i), 0, i)})
  })
  if(similarity==T) return(res)
  else return(1-res)
}


#' Run similarity network fusion
#'
#' @param dataL list( t(mRNA.snf.df), t(methylation.snf.df), t(cnv.snf.df) )
#' @param alpha Default 0.5. hyperparameter, usually (0.3~0.8)   Variance for local model
#' @param K Default 20. Number of neighbors, must be greater than 1. usually (10~30)  20
#' @param Iterations T.Default 20. Number of Iterations, usually (10~50)
#' @param dist.method Default Euclidean. pearson, spearman, kendall
#' @param survival Must a data frame. colnames OS.time, OS.event, RFS.time, RFS.event. Rownames must be sample name. Two columns or four columns.
#' @param max.cluster
#' @param std.normalize Default TRUE
#' @param cnv.index Default 0. If CNV data index is specified, Euclidean will be used to calculate distance
#'
#' @return
#' @export
#'
#' @examples run_SNF( list( t(mRNA.snf.df), t(methylation.snf.df), t(cnv.snf.df) ),  alpha = 0.5, K = 20, Iterations = 20    )
#' https://cran.r-project.org/web/packages/SNFtool/SNFtool.pdf
#' Distance reference: https://www.rdocumentation.org/packages/bioDist/versions/1.44.0
run_SNF <- function(dataL = NULL, alpha = 0.5, K = 20, Iterations = 20, dist.method ="Euclidean", survival=NULL, max.cluster = 5, std.normalize = TRUE, cnv.index = 0){

  library(SNFtool)
  library(bioDist)
  library(survival)

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
  }else if(dist.method == "Chi-squared"){
    dataL.dist = lapply(dataL.normalized, function(x) (SNFtool::chiDist2(x, x)^(1/2)))   # Chi-squared
  }else if(dist.method %in% c("pearson", "kendall", "spearman") ){
    dataL.dist = lapply(dataL.normalized, function(x) (1-cor(t(x), method = dist.method) )   )
  }else{
    stop("Distance: Euclidean, pearson, spearman, kendall")
  }

  # If CNV data index is specified, Euclidean will be used to calculate distance
  if(cnv.index!=0){
    dataL.dist[[cnv.index]] = SNFtool::dist2(dataL.normalized[[cnv.index]], dataL.normalized[[cnv.index]])^(1/2)
  }



  # Construct the similarity graphs
  affinityL = lapply(dataL.dist, function(x) affinityMatrix(as.matrix(x), K, alpha) )
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



#' Ranks each features by NMI (normalized mutual information) based on their clustering assingments
#'
#' @param data.list A list, pls specify the name. Element is data.frame (Row is sample, column is feature.). dataL <- list( mRNA=t(mRNA.snf.df), miRNA=t(miRNA.snf.df), Met=t(methylation.snf.df), CNV=t(cnv.snf.df) ).

#' @param wall
#'
#' @return
#' @export
#'
#' @examples
calculateNMI <- function(data.list, wall){

  nmi.res <- rankFeaturesByNMI(data.list, wall)

  if( is.null(names(nmi.res[[1]]) )  ){
    stop("Pls specify names for the list")
  }

  names(data.list) = nmi.res[[1]]
  names(data.list) = nmi.res[[2]]

  res = list()
  res$NMI.Raw = nmi.res

  library(foreach)
  omic.nmi.res <- foreach(omic.name = names(data.list), .combine = rbind) %do%{
       omic.nmi  = nmi.res[[1]][[omic.name]]
       omic.rank = nmi.res[[2]][[omic.name]]

       omic.data.frame = data.frame(
         Feature = colnames(data.list[[omic.name]]),
         NMI = omic.nmi,
         Omic = omic.name,
         OmicLevelRank = omic.rank, stringsAsFactors = FALSE
       )
       omic.data.frame
  }

  omic.nmi.res$OverallRank = rank(-omic.nmi.res$NMI, ties.method = "first")
  library(dplyr)
  omic.nmi.res = arrange(omic.nmi.res,OverallRank)

  res$Overall.NMI = omic.nmi.res

  rm(omic.nmi.res)

  res

}



