
#' Download dataset from GEO by accession ID and platform
#'
#' @param geo.accession.id GEO Accession ID
#' @param platform Platform
#' @param destdir Default tempdir()
#' @param platform.available If GLP platform not available
#' @param normalizeBetweenArrays Default FALSE. If perform normlize between arrays. Check boxplot first
#' @param id.mapping Default TRUE in order to select the maximum expression value if a gene has muiltiple probes
#'
#' @return list(expression, phenotype, probe.annotation)
#' @export
#'
#' @examples
download.geo.dataset <- function(geo.accession.id, platform = NULL, destdir = "~/GSE", platform.available = TRUE, normalizeBetweenArrays = FALSE, id.mapping = TRUE) {
  if (missing("geo.accession.id")) {
    stop("Please provide geo.accession.id")
  }

  if (!require(AnnoProbe)) {
    install.packages("AnnoProbe")
  }
  library(GEOquery)


  # load series and platform data from GEO
  gset <- getGEO(geo.accession.id, GSEMatrix = TRUE, AnnotGPL = platform.available, getGPL = TRUE, destdir = destdir)

  if (!is.null(platform)) {
    if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  } else {
    if (length(gset) > 1) stop("More than 1 dataset in this GSE")
    idx <- 1
  }

  gset <- gset[[idx]]

  platform <- gset@annotation
  warning(
    "-------------------\n",
    "Pls check paltform is, ", platform, "\n",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo.accession.id,
    "\n-------------------\n"
  )

  gpl.annotation <- fData(gset)

  # show(gset)

  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  phenotype <- pData(gset)

  # eliminate samples
  exp.df <- exprs(gset)

  # check if need to perform normalize
  boxplot(exp.df, outline = FALSE, notch = T)

  if (normalizeBetweenArrays) {
    exp.df <- limma::normalizeBetweenArrays(exp.df)
    boxplot(exp.df, outline = FALSE, notch = T)
  }

  # if log2 transform.    VALUE	Normalized log2 data represent microRNA expression levels
  qx <- as.numeric(quantile(exp.df, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    exp.df[which(exp.df <= 0)] <- NaN
    exp.df <- log2(exp.df)
    print("Perform log2 transformation")
  } else {
    print("Note: here not perform log2 transformation")
  }
  rm(qx, LogC)


  result <- list(
    expression = exp.df,
    rawExpression = exprs(gset),
    phenotype = phenotype,
    probe.annotation = gpl.annotation
  )




  # 2022-06-10 updated, great
  if (id.mapping) {
    ids <- AnnoProbe::idmap(platform, type = "soft")
    # only keep the max value one. Down by AnnoProbe
    result$geneExpressionMaxValue <- AnnoProbe::filterEM(exp.df, ids)

    result$id.mapping <- ids
  }

  return(result)
}




#' GEO dataset analysis by limma.
#'
#' @param exp
#' @param group
#' @param check.log2
#' @param normalizeBetweenArrays
#'
#' @return Normalized expression data.frame and differential analysis result
#' @export
#'
#' @examples
#' # Code from GEO2E. Thanks
#' library(dplyr)
#' GSE32575 <- loonR::download.geo.dataset("GSE32575", "GPL6102")
#' control <- GSE32575$phenotype %>%
#'   filter(characteristics_ch1.2 == "disease state: obese before bariatric surgery") %>%
#'   pull(geo_accession)
#' experiment <- GSE32575$phenotype %>%
#'   filter(characteristics_ch1.2 == "disease state: obese after bariatric surgery") %>%
#'   pull(geo_accession)
#' exp.df <- GSE32575$expression[, c(control, experiment)]
#' group <- rep(c("Before", "After"), c(length(control), length(experiment)))
#' loonR::microarray_limma_differential(exp.df, group)
microarray_limma_differential <- function(exp, group, check.log2 = TRUE, normalizeBetweenArrays = TRUE) {
  print(Sys.time())
  # -------- # -------- # -------- # -------- # -------- #
  # Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
  ###
  #   Differential expression analysis with limma
  library(GEOquery)
  library(limma)
  library(umap)

  group <- factor(group, levels = unique(group), labels = c("Ctrl", "Exp"))

  # load series and platform data from GEO
  # group membership for all samples
  sml <- group
  gset <- exp

  # box-and-whisker plot
  ord <- order(sml) # order samples by group
  par(mar = c(7, 4, 2, 1))
  boxplot(gset[, ord], boxwex = 0.6, notch = T, main = "Expression before normalization", outline = FALSE, las = 2, col = factor(sml, levels = unique(sml))[ord])
  legend("topleft", legend = unique(sml), fill = palette(), bty = "n")


  # log2 transformation
  if (check.log2) {
    print("- Check log2")
    ex <- gset
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
    LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

    if (LogC) {
      print("- Perform log2 transformation")
      ex[which(ex <= 0)] <- NaN # 小于等于的标为NaN
      gset <- log2(ex)
    } else {
      print("- No need to perform log2 transformation")
    }
  }

  if (normalizeBetweenArrays) {
    print("- Normalization, if you don't want, pls set normalizeBetweenArrays=FALSE")
    gset <- normalizeBetweenArrays(gset) # normalize data
  } else {
    print("- No normalization")
  }


  # assign samples to groups and set up design matrix
  print("- assign samples to groups and set up design matrix")
  gs <- factor(sml, levels = unique(sml))

  print("- limma")
  design <- model.matrix(~ gs + 0)
  colnames(design) <- levels(gs)

  fit <- lmFit(gset, design) # fit linear model

  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(unique(sml)[2], "-", unique(sml)[1], sep = "")
  cont.matrix <- makeContrasts(contrasts = cts, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)

  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust = "BH", sort.by = "p", number = Inf, coef = cts)
  tT$REF <- row.names(tT)

  # -------- # -------- # -------- # -------- # -------- #
  # Visualize and quality control test results.
  # Build histogram of P-values for all genes. Normal test
  # assumption is that most genes are not differentially expressed.
  tT2 <- topTable(fit2, adjust = "BH", sort.by = "B", number = Inf)
  hist(tT2$adj.P.Val,
    col = "grey", border = "white", xlab = "P-adj",
    ylab = "Number of genes", main = "P-adj value distribution"
  )


  # summarize test results as "up", "down" or "not expressed"
  dT <- decideTests(fit2, adjust.method = "BH", p.value = 0.05)
  # Venn diagram of results
  vennDiagram(dT, circle.col = palette())

  # create Q-Q plot for t-statistic
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")

  # volcano plot (log P-value vs log fold change)
  colnames(fit2) # list contrast names
  ct <- 1 # choose contrast of interest
  volcanoplot(fit2,
    coef = ct, main = colnames(fit2)[ct], pch = 20,
    highlight = length(which(dT[, ct] != 0)), names = rep("+", nrow(fit2))
  )

  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  plotMD(fit2, column = ct, status = dT[, ct], legend = F, pch = 20, cex = 1)
  abline(h = 0)

  ####
  # General expression data analysis
  ex <- gset

  # box-and-whisker plot
  ord <- order(gs) # order samples by group

  par(mar = c(7, 4, 2, 1))

  boxplot(ex[, ord], boxwex = 0.6, notch = T, main = "Expression after log2 and normalization", outline = FALSE, las = 2, col = gs[ord])
  legend("topleft", legend = unique(sml), fill = palette(), bty = "n")

  # expression value distribution
  par(mar = c(4, 4, 2, 1))
  plotDensities(ex, group = gs, main = "Expression value distribution", legend = "topright")

  # UMAP plot (dimensionality reduction)
  ex <- na.omit(ex) # eliminate rows with NAs
  ex <- ex[!duplicated(ex), ] # remove duplicates
  ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
  par(mar = c(3, 3, 2, 6), xpd = TRUE)
  plot(ump$layout, main = "UMAP plot, nbrs=3", xlab = "", ylab = "", col = gs, pch = 20, cex = 1.5)
  legend("topright",
    inset = c(-0.15, 0), legend = levels(gs), pch = 20,
    col = 1:nlevels(gs), title = "Group", pt.cex = 1.5
  )
  library("maptools") # point labels without overlaps
  pointLabel(ump$layout, labels = rownames(ump$layout), method = "SANN", cex = 0.6)

  # mean-variance trend, helps to see if precision weights are needed
  plotSA(fit2, main = "Mean variance trend")

  print(Sys.time())

  list(
    expr.df = gset,
    diff.res = tT
  )
}
