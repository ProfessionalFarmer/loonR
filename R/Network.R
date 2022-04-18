

#' WGCNA, Weighted correlation network analysis
#'
#' @param datExpr row is sample
#' @param datTraits 0 1 to indicate group. Numeric to indicate trait
#' @param power Default by sft$powerEstimate
#' @param minModuleSize We like large modules, so we set the minimum module size relatively high:
#' @param MEDissThres We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
#' @param traitColumnIndex Trait we are interested in
#' @param nSelectHeatmap n (default 1000) genes to show in heatmap (a way to visualize network)
#'
#' @return
#' @export
#'
#' @examples
runWGCNA <- function(datExpr, datTraits, traitColumnIndex = 1,power=NULL, minModuleSize = 30, MEDissThres = 0.25, nSelectHeatmap = 10000){

  ## 查看是否有离群样品
  sampleTree = hclust(dist(datExpr), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


  res = list()
  # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

  if(!require(WGCNA)){
    BiocManager::install("WGCNA")
    require(WGCNA)
  }
  #################################### 1
  # Step-by-step network construction and module detection
  #################################### 1
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



  # 2.b.2 Co-expression similarity and adjacency
  ## Detect outlier sample
  #sampleTree = hclust(dist(datExpr), method = "average")
  # plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  #
  if(is.null(power)){
    softPower = sft$powerEstimate
    if(is.null(softPower)){
      softPower = 6
    }
  }else{
    softPower = power
  }
  adjacency = adjacency(datExpr, power = softPower);

  # 2.b.3
  # Topological Overlap Matrix (TOM)
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM

  # 2.b.4 Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);

  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = minModuleSize;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  cat("Label 0 is reserved for unassigned genes.")
  table(dynamicMods)

  # We now plot the module assignment under the gene dendrogram
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")

  # 2.b.5 Merging of modules whose expression profiles are very similar
  # The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
  # such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
  # calculate their eigengenes and cluster them on their correlation
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")

  # We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
  # MEDissThres = 0.25
  MEDissThres = MEDissThres
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;

  # To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
  sizeGrWindow(12, 9)
  #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)

  # In the subsequent analysis, we will use the merged module colors in mergedColors. We save the relevant variables for
  # use in subsequent parts of the tutorial:
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;



  ######################################################### 2
  # Relating modules to external clinical traits and identifying important genes
  ######################################################### 2
  # 3.a Quantifying module–trait associations
  # we simply correlate eigengenes with external traits and look for the most significant associations
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

  res$`MEs(moduleEigengenes)` = MEs
  res$moduleTraitCor = moduleTraitCor
  res$moduleTraitPvalue = moduleTraitPvalue

  # Since we have a moderately large number of modules and traits, a suitable graphical representation will help in
  # reading the table. We color code each association by the correlation value
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))


  # 3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
  # We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as
  # (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
  # measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This
  # allows us to quantify the similarity of all genes on the array to every module.

  # Define variable weight containing the weight column of datTrait
  # trait.df is made by myself
  trait.df = as.data.frame(datTraits[,traitColumnIndex]);
  names(trait.df) = colnames(datTraits)[traitColumnIndex]
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, trait.df, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="");
  names(GSPvalue) = paste("p.GS.", names(trait.df), sep="");

  # 3.c Intramodular analysis: identifying genes with high GS and MM
  # Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
  # membership in interesting modules. As an example, we look at the brown module that has the highest association
  # with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
  # module membership MM, GS gene significance

  res$modNames = modNames
  res$moduleColors = moduleColors

  GSvsMM.list = lapply(modNames, function(module){
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    tmp.df = data.frame(
      geneModuleMembership  = abs(geneModuleMembership[moduleGenes, column]),
      geneTraitSignificance = abs(geneTraitSignificance[moduleGenes, 1]),
      gene = rownames(geneTraitSignificance)[moduleGenes],
      row.names = rownames(geneTraitSignificance)[moduleGenes]
    )
    tmp.df
  })
  names(GSvsMM.list) = modNames
  res$GSvsMM.list = GSvsMM.list

  ############################################### 3
  ## Visualization of networks within R
  ############################################### 3
  # 5.a Visualizing the gene network by heatmap
  nSelect = nSelectHeatmap
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

  # 5.b Visualizing the network of eigengenes
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  # Isolate weight from the clinical traits
  trait.df = as.data.frame(datTraits[,traitColumnIndex]);
  names(trait.df) = colnames(datTraits)[traitColumnIndex]
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, trait.df))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

  res$MET = MET
  # Export of networks to external software
  # Exporting to Cytoscape
  # Recalculate topological overlap if needed
  TOM = TOMsimilarityFromExpr(datExpr, power = softPower);

  ctyscape.module.list = lapply(modNames, function(module){

    # Select module probes
    probes = colnames(datExpr)
    inModule = is.finite(match(moduleColors, module));
    modProbes = probes[inModule];
    # Select the corresponding Topological Overlap
    modTOM = TOM[inModule, inModule];

    dimnames(modTOM) = list(modProbes, modProbes)
    modTOM = reshape2::melt(modTOM, value.name = "Edge", na.rm = TRUE)
    modTOM
  })
  names(ctyscape.module.list) = modNames
  res$ctyscape.module.list = ctyscape.module.list

  res

}



#' RTN Reconstruction of Transcriptional regulatory Networks analysis
#'
#' @param TFs transcription factors. vector
#' @param expData row is gene, include TFs and genes
#' @param group
#' @param cores Default 1
#' @param nPermutations Default 1000
#' @param phenotype A named vector, logFC, include TFs and genes
#' @param hits
#'
#' @return
#' @export
#'
#' @examples
RTN.analysis <- function(TFs = NULL, expData = NULL, group = NULL, cores=1, nPermutations = 1000, phenotype = NULL, hits = NULL){
  # https://www.bioconductor.org/packages/release/bioc/vignettes/RTN/inst/doc/RTN.html#overview

  if(!require(RTN)){
    BiocManager::install("RTN")
    library(RTN)
  }
  library(parallel)
  library(doParallel)

  if(is.null(TFs) | is.null(expData) | is.null(group) | is.null(hits) ){
    stop("Pls specify all required option")
  }

  if(is.null(phenotype)){
    warning("User didn't provide phenotype. We will use limma to calculate logFC")
    diff.res = loonR::limma_differential(expData, group, pre.filter = 0, cal.AUC = F)
    diff.res = diff.res[rownames(expData),]
    phenotype = diff.res$logFC
    rm(diff.res)
    names(phenotype) = diff.res$REF
  }


  rowAnnotation = data.frame(
    ID = row.names(expData),
    GENEID = rownames(expData),
    SYMBOL = rownames(expData),
    stringsAsFactors = F
  )

  colAnnotation = data.frame(
    row.names = colnames(expData),
    IDs = colnames(expData),
    Group = group
  )

  ####################### Transcriptional Network Inference (TNI)
  ####################### Transcriptional Network Inference (TNI)
  rtni <- tni.constructor(expData = as.matrix(expData),
                          regulatoryElements = TFs,
                          rowAnnotation = rowAnnotation,
                          colAnnotation = colAnnotation)



  registerDoParallel(cores=cores)
  parallel::mcaffinity(c(1:cores))
  ################################## run
  rtni <- tni.permutation(rtni, nPermutations = nPermutations, pValueCutoff = 0.05)
  #Unstable interactions are subsequently removed by bootstrap analysis
  parallel::mcaffinity(c(1:cores))
  rtni <- tni.bootstrap(rtni)
  #scans all triplets formed by two regulators and one target and removes the edge with the smallest MI value of each triplet, which is regarded as a redundant association.
  rtni <- tni.dpi.filter(rtni)
  ############################# stop

  # a list with regulons, including the weight assigned for each interaction
  regulons <- tni.get(rtni, what = "regulons.and.mode", idkey = "SYMBOL")


  # Compute regulon activity for individual samples
  rtni1st <- tni.gsea2(rtni, regulatoryElements = TFs, minRegulonSize = 5)
  regact <- tni.get(rtni1st, what = "regulonActivity")

  # Get sample attributes from the 'rtni1st' dataset
  col_annot <- tni.get(rtni1st, "colAnnotation")
  # Get ER+/- and PAM50 attributes for pheatmap
  col_annot <- data.frame(
    Group = as.factor( col_annot[,c("Group")] ),
    row.names = col_annot$IDs
  )

  regact$differential = as.data.frame(lapply(regact$differential,
                                             function(x) as.numeric(as.character(x))))

  # # Plot regulon activity profiles
  # regulon.activity.profiles =
  # pheatmap::pheatmap(
  #     data.frame( t(regact$differential) ),
  #     main="",
  #     annotation_col = col_annot,
  #     show_colnames = FALSE, annotation_legend = F,
  #     clustering_method = "ward.D2", fontsize_row = 6,
  #     clustering_distance_rows = "correlation",
  #     clustering_distance_cols = "correlation")


  ##################### Transcriptional Network Analysis (TNA)
  ##################### Transcriptional Network Analysis (TNA)
  # Input 1: 'object', a TNI object with regulons
  # Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
  # Input 3: 'hits', a character vector, usually a set of differentially expressed genes
  # Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype

  phenoIDs = data.frame(
    PROBEID = names(phenotype),
    ENTREZ  = names(phenotype),
    SYMBOL  = names(phenotype),
    stringsAsFactors = FALSE
  )

  rtna <- tni2tna.preprocess(object = rtni,
                             phenotype = phenotype,
                             hits = hits,
                             phenoIDs = phenoIDs)

  # Run the MRA method
  rtna <- tna.mra(rtna)
  #..setting 'ntop = -1' will return all results, regardless of a threshold
  mra <- tna.get(rtna, what="mra", ntop = -1)

  # Run the GSEA-2T method
  rtna <- tna.gsea2(rtna, nPermutations = nPermutations, minRegulonSize = 5)

  # Get GSEA-2T results
  gsea2 <- tna.get(rtna, what = "gsea2", ntop = -1)

  res = list(
    rtni = rtni,
    regulons = regulons,
    rtna = rtna,
    mra = mra,
    gsea2 = gsea2,
    regact = regact
  )

  res
}

#' Plot RTN network
#'
#' @param target.TFs which TFs to show. A vector
#' @param TF.fc named vector with TF fold change. Signaficantly differential expression
#' @param gene.fc named vector with gene fold change. Signaficantly differential expression
#' @param rtni.obj
#' @param hits Gene of interest, E.g. EMT gene
#'
#' @return
#' @export
#'
#' @examples
RTN.network.plot <- function(rtni.obj, target.TFs, TF.fc, gene.fc, hits, font.size=20){

  # https://www.bioconductor.org/packages/devel/bioc/vignettes/RTN/inst/doc/RTN.html#regulon-activity-profiles

  if(missing(rtni.obj) | missing(target.TFs) | missing(TF.fc) | missing(gene.fc) |  missing(hits) ) {
    stop("All options required")
  }

  library(RedeR)
  library(igraph)
  library(RTN)

  g <- tni.graph(rtni.obj, regulatoryElements = target.TFs )
  # The next chunk shows how to plot the igraph-class object using RedeR (Figure 1).

  ######################################### Start modification
  TF.fc <- TF.fc[intersect(V(g)$name, target.TFs)] # only show target.TFs
  TF.cols <- c(colorRampPalette(c("mediumseagreen", "white", "OrangeRed"))(10))
  names(TF.cols) <- c(1:10)

  TF.fc[TF.fc > 1] = 0.9999
  TF.fc[TF.fc < -1] = -0.9999

  TF.labs <- cut(TF.fc, breaks= seq(-1,1,0.2), labels=1:10)
  V(g)$nodeColor[match(names(TF.fc), V(g)$name)] <- TF.cols[TF.labs]

  ### gene color
  gene.fc <- gene.fc[intersect(V(g)$name, names(gene.fc))]
  gene.cols <- c(colorRampPalette(c("DeepSkyBlue", "white", "Orange"))(10))
  names(gene.cols) <- c(1:10)

  gene.fc[gene.fc > 1] = 0.9999
  gene.fc[gene.fc < -1] = -0.9999

  gene.col.labs <- cut(gene.fc,
                     breaks = seq(-1,1,0.2),
                     labels = 1:10)
  V(g)$nodeColor[match(names(gene.fc), V(g)$name)] <- gene.cols[gene.col.labs]

  ### node size
  V(g)$nodeSize <- rep(20, length(V(g)$name))
  V(g)$nodeSize[which(V(g)$name %in% target.TFs)] <- 56
  V(g)$nodeSize[which(V(g)$name %in% hits)] <- 40

  ### node shape
  V(g)$nodeShape[which(V(g)$name %in% hits)] <- "TRIANGLE"
  V(g)$nodeShape[which(V(g)$name %in% target.TFs)] <- "DIAMOND"


  ### only show hits and TF node name
  V(g)$nodeFontSize[! V(g)$name %in% c(target.TFs, hits)] = 0
  V(g)$nodeFontSize[  V(g)$name %in% c(hits)] = font.size
  V(g)$nodeFontSize[  V(g)$name %in% c(target.TFs)] = font.size * 2


  ################################### End modification

  rdp <- RedPort()
  calld(rdp)
  addGraph(rdp, g, layout=NULL)
  addLegend.color(rdp, g, type="edge")
  addLegend.shape(rdp, g)
  relax(rdp, ps = TRUE)


}








