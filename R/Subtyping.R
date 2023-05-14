

#' Fill NA in symmetric table
#'
#' @param symmetric.df
#'
#' @return
#' @export
#'
#' @examples
fillSymmetricNATable <- function(symmetric.df){

  if(sum(colnames(symmetric.df)!=rownames(symmetric.df))!=0){
    stop("Row name is not the same with colname")
  }

  for(r in 1:nrow(symmetric.df)){
    for(c in 1:ncol(symmetric.df)){
      value = symmetric.df[r,c]
      rev_value = symmetric.df[c,r]

      if(is.na(value) & !is.na(rev_value)){
        symmetric.df[r,c] = rev_value
      }
      if(!is.na(value) & is.na(rev_value)){
        symmetric.df[c,r] = value
      }

    }

  }
  data.frame(symmetric.df, check.names = F)

}

#' Hyper geometrix test for contigency table
#'
#' @param df data.frame with row and col names, and numeric value
#' @param lower.tail
#'
#' @return
#' @export
#'
#' @examples
#' df = data.frame(a=c(1,4),b=c(2,19))
#' rownames(df) = c("c","d")
#' loonR::hyperGeoTest4ContingencyTable(df)
hyperGeoTest4ContingencyTable <- function(df, lower.tail = FALSE){

  geomatrix = loonR::convertDfToNumeric(df)
  # perform geometrix,把p值放在相同矩阵的数据框中
  tmpgeo = matrix(nrow=length(row.names(geomatrix)),ncol=length(colnames(geomatrix)))
  colnames(tmpgeo) = colnames(geomatrix)
  rownames(tmpgeo) = rownames(geomatrix)


  for(i in row.names(tmpgeo)  ){ # row
    for(j in colnames(tmpgeo) ){  # column
      # 白球的个数，白球的总个数，黑球的总个数，抽的球（不是黑球，是球）个数
      p = phyper(geomatrix[i,j]-1, sum(geomatrix[i,]), sum(geomatrix)-sum(geomatrix[i,]), sum(geomatrix[,j]), lower.tail = lower.tail   )
      tmpgeo[i,j] = p
    }
  }
  tmpgeo = as.data.frame(tmpgeo)

  tmpgeo

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
#' @param adjusted.p If to adjust P value c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' @param remove.na If to remove NA
#' @param bar.color palatte for barplot
#' @param chiq.square.test
#' @param not.consider A vector that not consider the inside element
#'
#' @return list(barplot, pval.df.log, pval.df, plot)
#' @export
#'
#' @examples
#' g1 = sample(c("G1","G2","G3"),10, replace = T)
#' g2 = sample(c("1","2","3"),10, replace = T)
#' loonR::hyperGeoTest(g1, g2, col.prefix = "E")
hyperGeoTest <- function(row.group, col.group, row.prefix = "", col.prefix = "", lower.tail = FALSE, title = "", log10.lowest = 5,
                         print.fig = TRUE, adjusted.p=NULL, remove.na=FALSE, bar.color="npg", chiq.square.test = T, not.consider = NA){

  # https://www.omicsclass.com/article/324
  # 1-phyper(抽取样本中属于“特定类别”的数量-1,总样本中“特定类别”的数量, 总样本数-总样本中“特定类别”的数量, 从总样本中随机抽取的数量,)
  # 也是直接用黑白球的思路来说的：“袋子里有黑球n个和白球m个，不放回的取出k个球，其中有白球q个”
  # phyper(q,m,n,k)

  library(dplyr)

  if(remove.na){
      na.index = is.na(row.group) | is.null(row.group) | row.group == "" | is.na(col.group) | is.null(col.group) | col.group == ""
      na.index = na.index | row.group %in% not.consider | col.group %in% not.consider

      row.group = row.group[!na.index]
      col.group = col.group[!na.index]
      rm(na.index)
  }

  message("Total samples: ", length(col.group))

  # 超几何检验，与原来的分组比较
  geomatrix  = unclass(table(row.group, col.group))
  colnames(geomatrix) = as.character(paste(col.prefix, colnames(geomatrix),sep="" ))
  rownames(geomatrix) = as.character(paste(row.prefix, rownames(geomatrix),sep="" ))


  # 20220616 add barplot
  geomatrix.melted = geomatrix
  # chisq test
  geomatrix.melted.chiq.res = chisq.test(geomatrix.melted)
  if(chiq.square.test){
    main = paste0(geomatrix.melted.chiq.res$method,": ",
           formatC(geomatrix.melted.chiq.res$p.value,digits = 2, format = "E") )
    message(main)
  }else{
    main = ""
  }
  # covert proportion
  geomatrix.melted = prop.table(geomatrix.melted, margin = 2)
  # melt
  geomatrix.melted = reshape2::melt(geomatrix.melted)

  library(ggplot2)
  library(ggpubr)

  barplot = ggbarplot(geomatrix.melted, x = "col.group", y = "value",  # color = "row.group",
                      fill = "row.group", position = position_stack(),
                      palette = loonR::get.palette.color(bar.color, alpha = 0.7 ), legend = "right"  ) +
    ylab("Proportion")

  barplot = ggpar(barplot, legend.title = row.prefix, xlab = col.prefix,
                  main = main,
                  font.main = c(8, "plain", "black")
  )


  #用在这个函数替代 reduce wheel, cheers
  tmpgeo = loonR::hyperGeoTest4ContingencyTable(geomatrix)
  # # perform geometrix,把p值放在相同矩阵的数据框中
  # tmpgeo = matrix(nrow=length(row.names(geomatrix)),ncol=length(colnames(geomatrix)))
  # colnames(tmpgeo) = colnames(geomatrix)
  # rownames(tmpgeo) = rownames(geomatrix)
  #
  #
  # for(i in row.names(tmpgeo)  ){ # row
  #   for(j in colnames(tmpgeo) ){  # column
  #     # 白球的个数，白球的总个数，黑球的总个数，抽的球（不是黑球，是球）个数
  #     p = phyper(geomatrix[i,j]-1, sum(geomatrix[i,]), sum(geomatrix)-sum(geomatrix[i,]), sum(geomatrix[,j]), lower.tail = lower.tail   )
  #     tmpgeo[i,j] = p
  #   }
  # }
  # tmpgeo = as.data.frame(tmpgeo)




  res <- list()
  res$barplot = barplot


  if(!is.null(adjusted.p)){
    res$raw.pval.df = tmpgeo

    longvector = unlist(tmpgeo)
    longvector = p.adjust(longvector, method = adjusted.p)
    longvector = as.data.frame(matrix(longvector, ncol = ncol(tmpgeo)))
    row.names(longvector) = row.names(tmpgeo)
    colnames(longvector) = colnames(tmpgeo)
    tmpgeo = longvector
  }


  tmpgeo.log = -log10(tmpgeo)
  tmpgeo.log[tmpgeo.log > log10.lowest ] = log10.lowest - 0.1


  # reture formated value
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
  res$plot = ComplexHeatmap::Heatmap(as.matrix(tmpgeo.log),
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


#' Multiple hyper geomatrix test, include 2xn, 2x2
#'
#' @param main.row.group
#' @param minor.col.group
#' @param row.prefix
#' @param col.prefix
#' @param lower.tail
#' @param title
#' @param log10.lowest
#' @param print.fig
#' @param adjusted.p
#' @param remove.na
#' @param bar.color
#' @param chiq.square.test
#' @param not.consider
#'
#' @return
#' @export
#'
#' @examples
#' g1 = sample(c("G1","G2","G3"),10, replace = T)
#' g2 = sample(c("1","2","3"),10, replace = T)
#' res = loonR::multipleHyperGeoTest(g1, g2, col.prefix = "E")
multipleHyperGeoTest <- function(main.row.group, minor.col.group, row.prefix = "", col.prefix = "", lower.tail = FALSE, title = "", log10.lowest = 5,
                                print.fig = FALSE, adjusted.p=NULL, remove.na=FALSE, bar.color="npg", chiq.square.test = T, not.consider = NA){

  twoXtwo = list()
  twoXN = list()

  #  main is the group you want to test
  main.row.group = paste0(row.prefix,main.row.group)
  minor.col.group = paste0(col.prefix,minor.col.group)

  # obtain class level
  main.groups = unique(main.row.group)
  minor.groups = unique(minor.col.group)




  # 2023-04-28 convert to 2x2 matrix
  for (main.group in main.groups) {
     for(minor.group in minor.groups){
          row.group = rep(paste0("non-",main.group), length(main.row.group))
          col.group = rep(paste0("non-",minor.group),length(minor.col.group))

          row.group[main.row.group  == main.group ] = main.group
          col.group[minor.col.group == minor.group] = minor.group

          tmp.res = loonR::hyperGeoTest(row.group, col.group,
                              lower.tail = lower.tail, title = title, log10.lowest = 5,
                              print.fig = print.fig, adjusted.p=adjusted.p, remove.na=remove.na,
                              bar.color=bar.color, chiq.square.test = chiq.square.test, not.consider = not.consider)
          tmp.name = paste0(main.group, "--", minor.group)
          twoXtwo[[tmp.name]] <- tmp.res
     }
  }

  # convert to 2xn
  for (main.group in main.groups) {
      row.group = rep(paste0("non-",main.group), length(main.row.group))
      col.group = minor.col.group

      row.group[main.row.group  == main.group ] = main.group

      tmp.res = loonR::hyperGeoTest(row.group, col.group,
                                    lower.tail = lower.tail, title = title, log10.lowest = 5,
                                    print.fig = print.fig, adjusted.p=adjusted.p, remove.na=remove.na,
                                    bar.color=bar.color, chiq.square.test = chiq.square.test, not.consider = not.consider)
      tmp.name = paste0(main.group)
      twoXN[[tmp.name]] <- tmp.res

  }


  res = list(twoXtwo = twoXtwo,
             twoXN = twoXN)
  res
}




#' Subtype association: Calculate Jaccard similarity coefficient and perform hypergeometrix test
#'
#' @param df row is sample, column is study subtype. Plsease note rownames is sample ID
#' @param concateStudy If to add study name (from colnames) into result
#' @param adjusted.hypergeometrix.p Default FALSE, if to adjust hypergeometrix P value by BH method.
#' @param cut.edge.byPval Default 0.05
#' @param print.message If print message
#'
#' @return
#' @export
#'
#' @examples
#' data(Subtype.matrix)
#' res <- loonR::subtypeAssociationAnalysis(Subtype.matrix)
#' res$hyperGeoTest.Analysis
subtypeAssociationAnalysis <- function(df, concateStudy = F, adjusted.hypergeometrix.p = F, cut.edge.byPval = 0.05, print.message = T){

  # https://www.nature.com/articles/nm.3967#Sec9
  # we employed a network-based approach. The network encodes on nodes the information of subtype
  # prevalence and on edges their association calculated on the basis of Jaccard similarity coefficient,
  # which is defined by the size of the intersection between two sample sets over the size of their union.
  # To quantify the statistical significance of subtype associations, we performed hypergeometric tests
  # for overrepresentation of samples classified to one subtype in another. The resulting P values were
  # adjusted for multiple hypotheses testing using the Benjamini-Hochberg (BH) method.
  # Using this approach, we built a network consisting of the total 27 subtypes defined in the six different subtyping systems, interconnected by 96 significant (BH-corrected, P value < 0.001) edges.

  # Process: Jaccard analysis -> hypergeometrix analysis -> combine -> filter edge by p value

  res = list()
  res$raw = df

  ####################################### Jaccard analysis
  if(print.message){cat("Perform jaccard analysis\n")}

  study.Com <- loonR::generateCombinations(colnames(df), size=2, repeats = T)

  library(foreach)

  jaccardAnalysis.res <- foreach(i=1:nrow(study.Com), .combine = rbind) %do% {
      study1.name = study.Com[i,1]
      study2.name = study.Com[i,2]

      study1.vector = as.vector(unlist(df[study1.name]))
      study2.vector = as.vector(unlist(df[study2.name]))

      foreach(subtype1 = unique(study1.vector), .combine = rbind) %do%{

        foreach(subtype2 = unique(study2.vector), .combine = rbind) %do%{

          ind = study1.vector != 0 | study2.vector != 0

          j.coef = jaccard::jaccard(study1.vector[ind]==subtype1, study2.vector[ind]==subtype2)

          if(concateStudy){
            c( paste0(study1.name,"-",subtype1), paste0(study2.name, "-", subtype2), j.coef)
          }else{
            c( subtype1, subtype2, j.coef )
          }

        }

      }

  }
  colnames(jaccardAnalysis.res) <- c("Subtype1","Subtype2","Jaccard.index")
  jaccardAnalysis.res <- data.frame(jaccardAnalysis.res, check.names = F)

  jaccardAnalysis.res = reshape2::dcast(jaccardAnalysis.res, Subtype1~Subtype2, value.var = 'Jaccard.index')
  rownames(jaccardAnalysis.res) <- jaccardAnalysis.res$Subtype1
  jaccardAnalysis.res <- data.frame( jaccardAnalysis.res[,-c(1)], check.names = F )
  # reorder
  jaccardAnalysis.res <- jaccardAnalysis.res[colnames(jaccardAnalysis.res),]
  jaccardAnalysis.res <- loonR::fillSymmetricNATable(jaccardAnalysis.res)

  res$JaccardAnalysis = jaccardAnalysis.res
  res$JaccardAnalysis = loonR::convertDfToNumeric(res$JaccardAnalysis)


  ####################################### Hypergeometrix analysis
  if(print.message){cat("Perform Hypergeometrix analysis\n")}
  study.Com <- loonR::generateCombinations(colnames(df), size=2, repeats = T)

  library(foreach)

  hyperGeoTest.res <- foreach(i=1:nrow(study.Com), .combine = rbind) %do% {

    study1.name = study.Com[i,1]
    study2.name = study.Com[i,2]

    study1.vector = as.vector(unlist(df[study1.name]))
    study2.vector = as.vector(unlist(df[study2.name]))

    if(adjusted.hypergeometrix.p){
      hyp.test = loonR::hyperGeoTest( study1.vector, study2.vector, print.fig = F, adjusted.p = "BH" )
    }else{
      hyp.test = loonR::hyperGeoTest( study1.vector, study2.vector, print.fig = F, adjusted.p = NULL )
    }


    foreach(subtype1 = unique(study1.vector), .combine = rbind) %do%{

      foreach(subtype2 = unique(study2.vector), .combine = rbind) %do%{

        hyper.adjusted.p = hyp.test$pval.df.notFormated[subtype1,subtype2]

        if(concateStudy){
          c( paste0(study1.name,"-",subtype1), paste0(study2.name, "-", subtype2), hyper.adjusted.p)
        }else{
          c( subtype1, subtype2, hyper.adjusted.p )
        }

      }
    }

  }


  colnames(hyperGeoTest.res) <- c("Subtype1","Subtype2","Hypergeometrix.P")

  hyperGeoTest.res <- data.frame(hyperGeoTest.res, check.names = F)

  hyperGeoTest.res = reshape2::dcast(hyperGeoTest.res, Subtype1~Subtype2, value.var = 'Hypergeometrix.P')
  rownames(hyperGeoTest.res) <- hyperGeoTest.res$Subtype1
  hyperGeoTest.res <- data.frame( hyperGeoTest.res[,-c(1)], check.names = F )
  # reorder
  hyperGeoTest.res <- hyperGeoTest.res[colnames(hyperGeoTest.res),]
  hyperGeoTest.res <- loonR::fillSymmetricNATable(hyperGeoTest.res)

  # diagonal value changed to 1
  hyperGeoTest.res[col(hyperGeoTest.res)==row(hyperGeoTest.res)] = 1

  # the same order
  res$hyperGeoTest.Analysis = hyperGeoTest.res[rownames(res$JaccardAnalysis), colnames(res$JaccardAnalysis)]
  res$hyperGeoTest.Analysis = loonR::convertDfToNumeric(res$hyperGeoTest.Analysis)


  ########################### Another style data frame combined with hypergeometrix test and jaccard similarity
  if(print.message){cat("Combining\n")}

  hyperGeoTest.res[lower.tri(hyperGeoTest.res, diag = T)] = NA
  jaccardAnalysis.res[lower.tri(jaccardAnalysis.res, diag = T)] = NA

  hyperGeoTest.res$Rowname = rownames(hyperGeoTest.res)
  jaccardAnalysis.res$Rowname = rownames(jaccardAnalysis.res)

  hyperGeoTest.res.melt = reshape2::melt(hyperGeoTest.res, id.vars="Rowname", value.name = "Hypergeometrix.P", na.rm = T)
  jaccardAnalysis.res.melt = reshape2::melt(jaccardAnalysis.res, id.vars="Rowname", value.name = "Jaccard.similarity", na.rm = T)

  hyperGeoTest.res.melt$Hypergeometrix.P <- as.numeric(hyperGeoTest.res.melt$Hypergeometrix.P)
  jaccardAnalysis.res.melt$Jaccard.similarity <- as.numeric(jaccardAnalysis.res.melt$Jaccard.similarity)

  combined.matrix = dplyr::full_join(jaccardAnalysis.res.melt, hyperGeoTest.res.melt, by = c("Rowname"="Rowname", "variable"="variable") )
  colnames(combined.matrix)[1:2] = c("Subtype1","Subtype2")
  combined.matrix$Jaccard.similarity = as.numeric(as.vector(combined.matrix$Jaccard.similarity))
  combined.matrix$Hypergeometrix.P = as.numeric(as.vector(combined.matrix$Hypergeometrix.P))
  combined.matrix$Subtype1 <- as.vector(combined.matrix$Subtype1)
  combined.matrix$Subtype2 <- as.vector(combined.matrix$Subtype2)


  res$Combined = combined.matrix



  ########################### Genereate clean network

  cut.res = list()
  cut.res$combined.matrix = dplyr::filter(combined.matrix, Hypergeometrix.P < cut.edge.byPval)

  cut.res$JaccardAnalysis = res$JaccardAnalysis
  cut.res$hyperGeoTest.Analysis = res$hyperGeoTest.Analysis

  cut.res$JaccardAnalysis[ !is.na(cut.res$JaccardAnalysis) ] = NA
  cut.res$hyperGeoTest.Analysis[ !is.na(cut.res$hyperGeoTest.Analysis) ] = NA

  for(i in 1:nrow(cut.res$combined.matrix)){

    subtype1 = cut.res$combined.matrix[i,c("Subtype1")]
    subtype2 = cut.res$combined.matrix[i,c("Subtype2")]

    cut.res$JaccardAnalysis[subtype1,subtype2] = cut.res$combined.matrix[i,c("Jaccard.similarity")]
    cut.res$JaccardAnalysis[subtype2,subtype1] = cut.res$combined.matrix[i,c("Jaccard.similarity")]

    cut.res$hyperGeoTest.Analysis[subtype1,subtype2] = cut.res$combined.matrix[i,c("Hypergeometrix.P")]
    cut.res$hyperGeoTest.Analysis[subtype2,subtype1] = cut.res$combined.matrix[i,c("Hypergeometrix.P")]

  }

  cut.res$adjacencyMatrix = cut.res$JaccardAnalysis
  cut.res$adjacencyMatrix[is.na(cut.res$adjacencyMatrix)] = 0
  cut.res$adjacencyMatrix[cut.res$adjacencyMatrix!=0] = 1

  res$cut.res = cut.res

  res
}



#' Identification of consensus subtypes
#'
#' @param df row is sample, column is study subtype. Plsease note rownames is sample ID
#' @param replicate Default 100
#' @param seed Default 1
#' @param adjusted.hypergeometrix.p Default 0.05
#' @param inflation Default 2
#' @param proportion Default 0.8
#' @param adjacencyMatrixCutoff Default NULL. Cutoff for adjacenty matrix.
#' @param subtype.prefix
#' @param pOradjacent Default "adjacent". Could be p or adjacent. Using p or adjacent to perform mcl cluster
#' @param inflationConsensus Default equal to inflation. inflation is for one boot while inflationConsensus is for subtyping
#'
#' @return
#' @export
#'
#' @examples
#' data("Subtype.matrix")
#' res = loonR::consensusSubtyping(Subtype.matrix, replicate = 50)
#' res$consensus.map
consensusSubtyping <- function(df, replicate=100, seed=1, proportion = 0.8, adjusted.hypergeometrix.p = F, inflation = 2, adjacencyMatrixCutoff = NULL, subtype.prefix="CMS", pOradjacent = "adjacent", inflationConsensus = 0){

  # https://www.nature.com/articles/nm.3967#Sec9
  # Identification of consensus subtypes. To identify consensus groups from the network of subtype association,
  # we used a consensus clustering approach involving the following steps.
  # (a) Network construction: using the approach described above, 80% of patient samples are randomly selected
  # to generate a network of subtype association.
  # (b) Network clustering: the network generated is partitioned into clusters using MCL (Markov cluster algorithm)11,12,
  # which is a scalable and efficient unsupervised cluster algorithm for networks.
  # (c) Cluster evaluation: steps (a) and (b) are repeated for n = 1,000 times.
  # From all clustering results, we calculated a 27 × 27 consensus matrix,
  # defined by the frequency that each pair of subtypes is partitioned into the same cluster.
  # On the basis of the consensus matrix, we assessed the robustness of each subtype with a stability score,
  # which is the average frequency that its within-cluster association with other subtypes is the same as predicted by MCL on the network generated with all samples.
  # For evaluation of clustering performance, we employed weighted Silhouette width (R package 'WeightedCluster'),
  # which extends Silhouette width by giving more weights to subtypes that are more representative of their assigned clusters.
  # Here, we used stability scores as weights to calculate weighted Silhouette width and took the median over all subtypes
  # as a measure of clustering performance, which was used to evaluate the optimal number of clusters.

  # Process: a) --> b) --> c) --> Consensus matrix - > Identify consensus subtype --> stability --> Identification of core consensus samples

  library(foreach)
  library(dplyr)

  library(doParallel)
  registerDoParallel(40)
  set.seed(100)

  library(utils)
  pb <- utils::txtProgressBar(style = 3)

  ########################################## a) --> b) --> c)
  outter.each.res <- foreach::foreach(i=1:replicate, .combine = rbind) %dopar% {
    set.seed(seed+i)
    new_df <- df %>% sample_n(  as.integer( nrow(df) * proportion)  )
    new_df_analysis <- loonR::subtypeAssociationAnalysis(new_df, adjusted.hypergeometrix.p = adjusted.hypergeometrix.p, print.message = F)

    if(pOradjacent=="p"){ # 20221220
      t = new_df_analysis$cut.res$hyperGeoTest.Analysis
      t[is.na(t)] = 1
      t = -1 * log10(t)
    }else if(pOradjacent=="adjacent"){
      t = new_df_analysis$cut.res$adjacencyMatrix
    }else{
      stop("Can not find ", pOradjacent)
    }


    mcl.res <- MCL::mcl(t, addLoops = T, inflation = inflation)

    subtype.names <- colnames(t)

    inner.each.res <- foreach::foreach(cluster=unique(mcl.res$Cluster[mcl.res$Cluster!=0]), .combine = rbind) %do% {
        # 20211018: Cluster 0 is not right
        cluster.subtypes = subtype.names[mcl.res$Cluster==cluster]
        if(length(cluster.subtypes)==1){
          c(cluster.subtypes,cluster.subtypes,1)
        }else{
          comb = loonR::generateCombinations(cluster.subtypes,size = 2, repeats = T)
          comb$Connected = 1
          comb
        }

    }
    colnames(inner.each.res) = c("Subtype1","Subtype2","Connected")
    inner.each.res$Boot = i

    setTxtProgressBar(pb, i/replicate)

    data.frame(inner.each.res)
  }

  close(pb)

  res = list()
  # 20230322 save parameter
  res$input = list(data.frame = df, replicate = replicate, seed = seed, proportion = proportion,
                   adjusted.hypergeometrix.p = adjusted.hypergeometrix.p,
                   inflation = inflation, adjacencyMatrixCutoff = adjacencyMatrixCutoff,
                   subtype.prefix = subtype.prefix, pOradjacent = pOradjacent, inflationConsensus = 0)

  res$raw = outter.each.res

  res$group.df = outter.each.res %>% dplyr::group_by(Subtype1,Subtype2) %>% dplyr::summarise( Freq= (n()/replicate) )
  rm(outter.each.res)

  ############################### construct consensus matrix from grouped dataframe
  consensusMatrix = reshape2::dcast(res$group.df, Subtype1~Subtype2, value.var = c("Freq"))
  rownames(consensusMatrix) = consensusMatrix$Subtype1
  consensusMatrix = consensusMatrix[,-c(1)]

  consensusMatrix = loonR::fillSymmetricNATable(consensusMatrix)
  consensusMatrix[is.na(consensusMatrix)] = 0

  # set diag value, pls note this comment.这个地方决定在做MCL之前是否要设置对角线为1
  # consensusMatrix = as.matrix(consensusMatrix)
  # diag(consensusMatrix) = 1
  consensusMatrix = data.frame( consensusMatrix, check.names = F )
  res$consensusMatrix = consensusMatrix


  ############################## adjacenty matrix
  res$adjacencyMatrix = consensusMatrix

  if(inflationConsensus == 0){ # default to inflation
    inflationConsensus = inflation
  }


  # Identify consensus subtype
  if(is.null(adjacencyMatrixCutoff)){
    subtyping.res <- loonR::identifySubtypeFromMatrix(consensusMatrix, usingRawDf = T, adjacency.cutoff = adjacencyMatrixCutoff, clusterPrefix = subtype.prefix, inflation = inflationConsensus)
    res$adjacencyMatrix.proceed = consensusMatrix
  }else{
    subtyping.res <- loonR::identifySubtypeFromMatrix(consensusMatrix, usingRawDf = F, adjacency.cutoff = adjacencyMatrixCutoff, clusterPrefix = subtype.prefix, inflation = inflationConsensus)
    res$adjacencyMatrix.proceed = consensusMatrix
    res$adjacencyMatrix.proceed[consensusMatrix > adjacencyMatrixCutoff] = 1
    res$adjacencyMatrix.proceed[consensusMatrix <= adjacencyMatrixCutoff] = 0

  }

  res$ConsensusSubtype.res = subtyping.res
  res$ConsensusSubtype.clean = subtyping.res$cluster.df

  # include study information
  t.df = t(df)
  t.df.melt = loonR::meltDataFrameByGroup(t.df,row.names(t.df))[,c(1,3)] %>% unique()
  colnames(t.df.melt) = c("Study","Subtype")

  res$ConsensusSubtype.clean = dplyr::full_join(res$ConsensusSubtype.clean, t.df.melt, by=c("Sample"="Subtype"))
  rm(t.df, t.df.melt)

  # 统计每个study的亚型的样本数目
  subtype.count = loonR::countClassByColumn(df)
  res$ConsensusSubtype.clean = dplyr::full_join(res$ConsensusSubtype.clean, subtype.count, by=c("Sample"="Class","Study"="Column"))
  rm(subtype.count)


  # 检查每个study的亚型是否分配的CMS，如果没有则报错终止
  if(sum(is.na(res$ConsensusSubtype.clean$Cluster))!=0){
    warning("Not found CMS subtype for the following type: ",res$ConsensusSubtype.clean$Sample[is.na(res$ConsensusSubtype.clean$Cluster)] )
    stop("Please check the study and stype")
  }

  res$CMSCount <- res$ConsensusSubtype.clean %>% dplyr::group_by(Cluster) %>% dplyr::summarise(SubtypeCount=n())

  res$CMSCount <- data.frame(res$CMSCount, row.names = res$CMSCount$Cluster) %>% dplyr::rename(Subtype=Cluster)



  ############################################## robustness of each subtype with a stability score
  res$group.df$ConsensusSubtype1 = subtyping.res$cluster.df$Cluster[match(res$group.df$Subtype1, subtyping.res$cluster.df$Sample)]
  res$group.df$ConsensusSubtype2 = subtyping.res$cluster.df$Cluster[match(res$group.df$Subtype2, subtyping.res$cluster.df$Sample)]

  res$group.df$SameGroup = res$group.df$ConsensusSubtype1 == res$group.df$ConsensusSubtype2

  res$stability = res$group.df %>% filter(SameGroup) %>% dplyr::rename(Subtype=ConsensusSubtype1) %>%
                            dplyr::group_by(Subtype) %>% dplyr::summarise(Stability=mean(Freq))



  ############################################# Identify core samples
  newlabels <- foreach(sample=rownames(df), .combine = rbind) %do%{

    sample.raw.subtype = unlist(df[sample,])
    # clean$sample其实是亚型
    sample.new.subtype = res$ConsensusSubtype.clean$Cluster[match(sample.raw.subtype, res$ConsensusSubtype.clean$Sample)]
    sample.new.subtype = as.character(sample.new.subtype)

    sig.count = 0
    core.sample = FALSE
    cms.type = NA

    uniq.cms = as.character(row.names(res$CMSCount))
    # 20230321 add p value for each CMS
    p.cms = c()

    for(t.cms in uniq.cms){

      # 也是直接用黑白球的思路来说的：“袋子里有黑球n个和白球m个，不放回的取出k个球，其中有白球q个”
      # phyper(q,m,n,k)
      # p = phyper(geomatrix[i,j]-1, sum(geomatrix[i,]), sum(geomatrix)-sum(geomatrix[i,]), sum(geomatrix[,j]), lower.tail = lower.tail   )

      t.pval = phyper(sum(sample.new.subtype==t.cms)-1,
                      as.numeric(  unclass(res$CMSCount[t.cms, c("SubtypeCount")])  ),
                      as.numeric(  unclass(sum(res$CMSCount$SubtypeCount)-res$CMSCount[t.cms, c("SubtypeCount")])  ),
                      length(sample.new.subtype), lower.tail = F   )

      if(t.pval <= 0.1){
         sig.count = sig.count + 1
         core.sample = TRUE
         cms.type = t.cms
      }

      # 20230321 add p value for each CMS
      p.cms = c(p.cms, t.pval)

    }
    if(sig.count>1){
      cms.type = "Confusing"
      core.sample = FALSE
    }

    # Count appreance frequency
    subtype.count = unlist(table(sample.new.subtype))
    if(sum(subtype.count == subtype.count[which.max(subtype.count)])!=1){
      HighFrequencySubtype = "Confusing"
    }else{
      HighFrequencySubtype = names(subtype.count)[which.max(subtype.count)]
    }
    p.cms = as.numeric( format(p.cms, digits = 3) )
    c(sample.new.subtype, core.sample, cms.type, HighFrequencySubtype, p.cms)
  }


  newlabels = data.frame(newlabels, row.names = rownames(df))
  colnames(newlabels) <- c( paste0("New",colnames(df)),
                            "CoreSample","CMSSubtype", "HighFrequencySubtype",
                            paste0("P.", uniq.cms) )

  # Select CMS subtype based on p value 20230321
  ttt = loonR::findMaxMinColumnNamesForEachRow(newlabels, ties.method = "first",
                                               min = T, specified.column = seq(ncol(newlabels)-3, ncol(newlabels) )  )
  newlabels$CMS.minP = stringr::str_remove_all(ttt$Min.ColName, "P.")
  newlabels$CMS.minP.value = ttt$Min.Value
  rm(ttt)

  newlabels = cbind(df, newlabels)
  res$Samples <- newlabels





  ############################################# END调整一下列名和其他地方
  colnames(res$ConsensusSubtype.clean)[1] <- c("Subtype")

  # 此时对角线一定要为1，表示节点自己的关联
  diag(res$consensusMatrix) = 1
  res$consensusMatrix = data.frame( res$consensusMatrix, check.names = F )


  ############################################## 20211018 add plot
  core.annotation.df = res$Samples %>% filter(HighFrequencySubtype!="Confusing")
  core.annotation.df = core.annotation.df[,c(colnames(df),"HighFrequencySubtype")] # select by HighFrequencySubtype
  colnames(core.annotation.df) = c(colnames(df),"CMS")

  res$core.annotation.plot = loonR::heatmap.annotation(
                                      group = core.annotation.df$CMS,
                                      annotation.df = core.annotation.df[,colnames(df)],
                                      sort.group = T)
  rm(core.annotation.df)

  ############################################# 20211019 add igraph and consensus map

  cluster.info = res$ConsensusSubtype.clean

  consensus.map = res$consensusMatrix
  consensus.map[lower.tri(consensus.map,diag = T)] = NA
  # filter by value control the edge
  res$consensus.map.for.cytoscape = loonR::meltDataFrameByGroup(consensus.map, rownames(consensus.map),variable_name = "Subtype2") %>% filter(value>0.0)

res$ConsensusSubtype.clean$Subtype

  plotNetwork <- function(){

    library(igraph)

    net <- graph.data.frame(res$consensus.map.for.cytoscape,
                            directed=FALSE,
                            vertices=cluster.info)

    # set node color点颜色
    # V(net)$color <- loonR::get.palette.color()[as.numeric(factor(V(net)$Study))]
    V(net)$color <- loonR::get.palette.color("aaas")[as.numeric(stringr::str_remove_all(V(net)$Cluster, subtype.prefix ) )]
    # 点border颜色
    V(net)$frame.color = loonR::get.palette.color("aaas")[as.numeric(stringr::str_remove_all(V(net)$Cluster, subtype.prefix ) )]


    # Set node size点大小
    V(net)$size <- 25

    # 字体颜色
    V(net)$label.color <- "black"
    # Setting them to NA will render no labels:
    # V(net)$label <- NA

    # Width
    E(net)$width <- log10(E(net)$value * 500)

    set.seed(1)
    plot(net, e=TRUE, v=TRUE)

  }
  res$plotNetwork = plotNetwork

  ######################################## 20211019 Consensus map
  sort.ind = order(res$ConsensusSubtype.clean$Cluster[ match(colnames(res$consensusMatrix), res$ConsensusSubtype.clean$Subtype)])

  res$consensus.map.plot = loonR::heatmap.with.lgfold.riskpro(
    res$consensusMatrix[sort.ind,sort.ind],
    res$ConsensusSubtype.clean$Cluster[ match(colnames(res$consensusMatrix), res$ConsensusSubtype.clean$Subtype)][sort.ind],
    show.lgfold = F, show.risk.pro = F, scale = F,
    specified.color = c("white","orange"), show_column_names = T, group.name = subtype.prefix
  )

  res

}



#' A matrix
#'
#' @param df n*n row and column are samples
#' @param adjacency.cutoff
#' @param inflation Default is 2 for mcl clustering
#' @param usingRawDf Default FALSE. If use raw data frame instead of adjacency matrix
#' @param clusterPrefix
#'
#' @return
#' @export
#'
#' @examples
#' data(Subtype.matrix)
#' sub.asso.ana.res <- loonR::subtypeAssociationAnalysis(Subtype.matrix)
#' cluster.res = loonR::identifySubtypeFromMatrix(sub.asso.ana.res$JaccardAnalysis, adjacency.cutoff = 0.8)
#' cluster.res
identifySubtypeFromMatrix <- function(df, adjacency.cutoff = 0.5, inflation = 2, usingRawDf=F, clusterPrefix=NULL){

  adjacency.matrix = df
  if(!usingRawDf){
     adjacency.matrix[adjacency.matrix > adjacency.cutoff] = 1
     adjacency.matrix[adjacency.matrix <= adjacency.cutoff] = 0
  }

   mcl.res <- MCL::mcl(adjacency.matrix, addLoops = F, inflation = inflation)

   cluster = mcl.res$Cluster
   sample = row.names(adjacency.matrix)

   k = mcl.res$K

   unique.cluster = unique(cluster)

   for(i in 1:length(unique.cluster)){
     cluster[cluster==unique.cluster[i]] = paste0("Cluster",i)
   }

   cluster = as.numeric( stringr::str_remove_all(cluster,"Cluster") )
   if(!is.null(clusterPrefix)){
     cluster = paste0(clusterPrefix, cluster)
   }

   names(cluster) = sample
   rm(unique.cluster)

   cluster.df = data.frame(Sample = sample, row.names = sample, Cluster = cluster)

   ##################################### use cluster label instead of 0/1
   # cluster network is similar to adjacency matrix
   cluster.network = adjacency.matrix
   cluster.network[!is.na(cluster.network)] = 0

   for(clu in unique(cluster)){
      samples <- names(cluster)[cluster==clu]
      if(length(samples)==1){
        cluster.network[samples,samples] = clu
      }else{
        sample.comb = loonR::generateCombinations(samples, size = 2, repeats = T)
        for(i in 1:nrow(sample.comb)){
          sample1 = sample.comb[i,1]
          sample2 = sample.comb[i,2]
          cluster.network[sample1, sample2] = clu
          cluster.network[sample2, sample1] = clu
        }
      }

   }


   res = list(
     mcl.res = mcl.res,
     K = k,
     sample = sample,
     cluster = cluster,
     cluster.df = cluster.df,
     adjacency.matrix = adjacency.matrix,
     cluster.network = cluster.network
   )

   res

}



#' Count the class by each column
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
#' data(Subtype.matrix)
#' loonR::countClassByColumn(Subtype.matrix)
countClassByColumn <- function(df){

  library(foreach)
  res <- foreach(cname=colnames(df), .combine = rbind ) %do%{

    c.count <- unclass(table(unlist(df[,cname])))
    c.count = data.frame(c.count, Class = names(c.count), row.names = names(c.count))
    colnames(c.count) = c("Count","Class")
    c.count$Column = cname
    c.count

  }

  res
}



#' Convert data frame to GenePattern gct file
#'
#' @param df Row is gene
#' @param file Default ./expression.gcf
#'
#' @return
#' @export
#'
#' @examples
dataframe2gct_file = function(df, file = "./expression.gcf", description = NULL, version = "#1.2" ){

  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3893799/

  # https://github.com/drmjc/metaGSEA/blob/master/DESCRIPTION
  # https://rdrr.io/github/drmjc/metaGSEA/src/R/export.gsea.gct.R
  message("Pls note, the output gct file path is ", file)

  if(is.null(description)){
    gct <- data.frame(Name=rownames(df), Description=rownames(df), df, stringsAsFactors=FALSE, check.names=FALSE)
  }else{
    gct <- data.frame(Name=rownames(df), Description=description, df, stringsAsFactors=FALSE, check.names=FALSE)
  }


  OUT <- file(file, "w")
  writeLines(version, OUT)
  writeLines(paste(nrow(df), ncol(df), sep="\t"), OUT)
  close(OUT)
  readr::write_delim(gct, file, delim = "\t", col_names = T, append = T)

  invisible( gct )

}


#' Convert a vector to GenePattern class file
#'
#' @param group A vector
#' @param file Default ./class.cls
#'
#' @return
#' @export
#'
#' @examples
group2class_file = function(group, file = "./class.cls"){

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3893799/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2065909/

# The first line of the CLS file contains three values: the number of samples,
# the number of classes, and the version number of file format (always 1).
# The second line begins with a pound sign (#) followed by a name for each class.
# The last line contains a class label for each sample. The number and order
# of the labels must match the number and order of the samples in the expression dataset.
# The class labels are sequential numbers (0, 1, . . .)
# assigned to each class listed in the second line.

  message("Pls note, the output class file path is ", file)

  version = 1

  OUT <- file(file, "w")
  writeLines( stringr::str_c( length(group), length(unique(group)), version, sep =" "), OUT)
  writeLines( stringr::str_c( "#", paste0( unique(group), collapse = " " ), sep = " "), OUT)
  writeLines( stringr::str_c( paste0(  as.numeric(as.factor(group)) - 1 , collapse = " " ) , sep = " "), OUT)

  close(OUT)

}

