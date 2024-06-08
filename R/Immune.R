#' Run estiamte
#'
#' @param expression.df   the row name in an input data must be gene symbols
#' @param platform c("affymetrix", "agilent", "illumina")
#'
#' @return Data frame with immune infiltration information
#' @export
#'
#' @examples
runEstimate <- function(expression.df, platform = c("affymetrix", "agilent", "illumina") ){
  # https://bioinformatics.mdanderson.org/estimate/rpackage.html
  #

  platform <- match.arg(platform)

  if(!require(estimate)){
    library(utils)
    rforge <- "http://r-forge.r-project.org"
    install.packages("estimate", repos=rforge, dependencies=TRUE)
    help(package="estimate")
  }


  library(estimate)
  tmpid = stringi::stri_rand_strings(1, 10)

  raw.file.path = paste0("./", tmpid,".rna.abundance.txt")
  filter.file.path = paste0("./", tmpid,".genes.gct")
  estimate.file = paste0("./", tmpid,".estimate.gct")

  write.table(expression.df, file = raw.file.path, quote = F, sep = "\t")

  # filter common genes
  filterCommonGenes(
    input.f = raw.file.path,
    output.f= filter.file.path,
    id="GeneSymbol"
  )
  file.remove(raw.file.path)


  # run estimate
  # input.ds character string specifying name of input GCT file containing stromal, immune, and estimate scores for each sample
  # output.ds character string specifying name of output file
  estimateScore(input.ds = filter.file.path,
                output.ds= estimate.file,
                platform = platform)
  file.remove(filter.file.path)

  # load data
  estimate.scores <- read.table(estimate.file, skip = 2, header = T)
  file.remove(estimate.file)

  estimate.scores = t(estimate.scores)
  colnames(estimate.scores) = estimate.scores[c(1),]
  estimate.scores = estimate.scores[-c(1,2),]

  rownames(estimate.scores) = colnames(expression.df)
  estimate.scores = data.frame(estimate.scores, check.names = F, stringsAsFactors = F)


  # StromalScorenumeric scalar specifying the presence of stromal cells in tumor tissue
  # ImmuneScorenumeric scalar specifying the level of infiltrating immune cells in tumor tissue
  # ESTIMATEScorenumeric scalar specifying tumor cellularity
  # TumorPuritynumeric scalar specifying ESTIMATE-based tumor purity with value in range[0,1]
  estimate.scores

}


#' run EPIC
#'
#' @param expression.df a matrix of the TPM (or RPKM) gene expression from the samples for which to estimate cell proportions
#'
#' @return Data frame with immune infiltration information
#' @export
#'
#' @examples
runEPIC <- function(expression.df){
 # https://github.com/GfellerLab/EPIC

  if(!require(EPIC)){
    devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
  }


  # out is a list containing the various mRNA and cell fractions in each samples as well as some data.frame of the goodness of fit.
  out <- EPIC::EPIC(bulk = expression.df)
  out
}



#' Run Immunedeconv
#'
#' @param expression.df gene × sample. TPM-normalized, not log-transformed
#' @param method c("quantiseq", "timer", "cibersort", "cibersort_abs", "mcp_counter", "xcell", "epic")
#' @param indication When use timer, should be one of immunedeconv::timer_available_cancers
#'
#' @return
#' @export
#'
#' @examples
runImmunedeconv <- function(expression.df, method = c("quantiseq", "timer", "cibersort", "cibersort_abs", "mcp_counter", "xcell", "epic"),
                            indication = ""){
  # https://github.com/icbi-lab/immunedeconv
  # https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html
  if(!require(immunedeconv)){
    remotes::install_github("icbi-lab/immunedeconv")
  }

  method <- match.arg(method)

  if(method=="timer"){
    res = immunedeconv::deconvolute(expression.df, method, indications = indication)
  }else{
    res = immunedeconv::deconvolute(expression.df, method)
  }

  res = t(res)
  colnames(res) = res[c(1),]
  res = res[-c(1),]

  res = data.frame(res, check.names = F, stringsAsFactors = F)

  res

}







#' Perform correlation analysis between two data.frame
#'
#' @param x.axis.df Row is sample, column variable will be column names in new correlation heatmap.
#' @param y.axis.df Row is sample, column variable will be row names in new correlation heatmap.
#' @param cor.method Defaut spearman
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
CorrelationHeatmapTwoDf <- function(x.axis.df, y.axis.df, cor.method = "spearman"){

  # check samples
  sample.names = intersect(rownames(x.axis.df), rownames(y.axis.df))
  if(length(sample.names)!=nrow(x.axis.df) & length(sample.names)!=nrow(y.axis.df) ){
    stop("Sample name is not the same")
  }

  # the same order
  x.axis.df = x.axis.df[sample.names,]
  y.axis.df = y.axis.df[sample.names,]


  # matrix to store rho
  cor.res = data.frame(matrix(nrow = length(colnames(y.axis.df)),
                              ncol = length(colnames(x.axis.df)) )
  )
  rownames(cor.res) = colnames(y.axis.df)
  colnames(cor.res) = colnames(x.axis.df)

  # matrix to store p value
  cor.p.res = cor.res

  # start from x colnames
  for(x.axis.label in colnames(x.axis.df) ) {

    x.axis.label.value = as.numeric(as.vector(x.axis.df[,x.axis.label]) )

    # start from y colnames
    for(y.axis.label in colnames(y.axis.df) )  {

      y.axis.label.value = as.numeric(as.vector(y.axis.df[,y.axis.label]) )

      cor_res = cor.test(x.axis.label.value, y.axis.label.value, method = cor.method,)

      cor.p.res[y.axis.label,x.axis.label] = cor_res$p.value
      cor.res[y.axis.label,x.axis.label]   = cor_res$estimate
    }

  }

  cor.p.res$Rnames = rownames(cor.p.res)
  cor.p.res = reshape2::melt(cor.p.res)

  cor.res$Rnames = rownames(cor.res)
  cor.res = reshape2::melt(cor.res)


  cor.merged = dplyr::full_join(cor.res, cor.p.res, by=c("Rnames"="Rnames", "variable"="variable"), suffix = c(".rho", ".p"))
  cor.merged$logp = -1 * log10(cor.merged$value.p)


  library(ggplot2)
  library(ggpubr)

  cor.merged = cor.merged %>% mutate(
    Sig=case_when(
      `value.p` > 0.05 ~ "",
      `value.p` <= 0.05 &  `value.p` > 0.01 ~ "*",
      `value.p` <= 0.01 &  `value.p` > 0.001 ~ "**",
      `value.p` <= 0.001 ~ "*"
    )
  )


  ggplot(cor.merged, aes(variable, Rnames, fill= value.rho, label = Rnames)) +
    geom_tile() +
    scale_fill_gradient2(low = "red",  mid = "white", high = "blue",  midpoint = 0, name = "Rho") + rotate_x_text(45) + ylab("") + xlab("") + geom_text(aes(label = Sig), col='black', cex=6)


}


#' Perform correlation analysis between two data.frame
#'
#' @param x.axis.df Row is sample, column variable will be point names.
#' @param y.axis.df Row is sample, column variable will be names in new correlation heatmap.
#' @param cor.method Defaut spearman
#' @param spec.select.list list(group1 = c(1,2,3), grou2 = c(4,5)) numbers are col indexin x.axis.df
#' @param rho.cutoff Default 0.3
#' @param p.color c("#D95F02", "#CCCCCC99", "#1B9E77")
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' library(vegan)
#' data(varespec)
#' data(varechem)
#'
#' loonR::CorrelationHeatmapTwoDfMantelTest(varespec[,c(1,2)], varechem)
#'
#' l = list(Community_Diversity= 1:5, Community_structure =6:14)
#' loonR::CorrelationHeatmapTwoDfMantelTest(varespec, varechem, spec.select.list = l)
#'
CorrelationHeatmapTwoDfMantelTest <- function(x.axis.df, y.axis.df, cor.method = "spearman", spec.select.list = NULL, rho.cutoff = 0.3, p.color = c("#D95F02", "#CCCCCC99", "#1B9E77") ){

  # wechat post: 像Science一样绘制相关性热图，看这一篇就够了
  # wechat post: R语言 | 终于实现了Mantel检验的相关性热图
  # check samples
  sample.names = intersect(rownames(x.axis.df), rownames(y.axis.df))
  if(length(sample.names)!=nrow(x.axis.df) & length(sample.names)!=nrow(y.axis.df) ){
    stop("Sample name is not the same")
  }

  if(is.numeric(cor.method)){
    stop("cor.method should be character while rho.cutoff is numeric")
  }

  # the same order
  if(ncol(x.axis.df)==1){ # in case of only one column
    clnm = colnames(x.axis.df)
    x.axis.df = data.frame(x.axis.df[sample.names,], row.names = sample.names)
    colnames(x.axis.df) = clnm
    rm(clnm)
  }else{
    x.axis.df = x.axis.df[sample.names,]
  }


  y.axis.df = y.axis.df[sample.names,]

  if(!require("vegan")){
    BiocManager::install("vegan")
    library("vegan")
  }
  if(!require("linkET")){
    dvetools:install_github("Hy4m/linkET")
    library("linkET")
  }
  library(ggplot2)

  if(is.null(spec.select.list)){

    # single x col is a spec
    # spec.select.list = lapply(1:ncol(x.axis.df), function(x){
    #   x
    # })
    # names(spec.select.list) = colnames(x.axis.df)
    #
    # library(dplyr)
    #
    # mantel <- linkET::mantel_test(x.axis.df, y.axis.df,
    #                               spec_select  = spec.select.list) %>%
    #   mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
    #                   labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    #          pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
    #                   labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

    # start from x colnames
    library(foreach)
    cor.res = foreach(x.axis.label = colnames(x.axis.df), .combine = rbind ) %do% {

      x.axis.label.value = as.numeric(as.vector(x.axis.df[,x.axis.label]) )

      # start from y colnames
      t = foreach(y.axis.label = colnames(y.axis.df), .combine = rbind ) %do% {

        y.axis.label.value = as.numeric(as.vector(y.axis.df[,y.axis.label]) )

        cor_res = cor.test(x.axis.label.value, y.axis.label.value, method = cor.method)

        c(x.axis.label, y.axis.label, cor_res$p.value, cor_res$estimate)
      }
      t
    }
    mantel = data.frame(cor.res)
    colnames(mantel) = c("x", "y", "p", "r")
    mantel$p = as.numeric(mantel$p)
    mantel$r = as.numeric(mantel$r)

    mantel = mantel %>%  mutate(rd = cut(r, breaks = c(-Inf, -1 * rho.cutoff, rho.cutoff, Inf),
                      labels = c(paste0("< -",rho.cutoff), paste0("-", rho.cutoff, " - ", rho.cutoff), paste0(">= ", rho.cutoff) ) ),
             pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                      labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


    test.method = cor.method

    # # # # # # # # # # # # # # # # # # # # # # # #
    plot = qcorrplot(correlate(y.axis.df, method = cor.method), type = "lower", diag = FALSE) +
      geom_square() +
      geom_couple(aes(colour = rd,
                      size = pd),
                  data = mantel, curvature = 0.1) +
      #geom_diag_label(mapping = aes(y = .y + 0.05),#定义对角线上的文字
      #                hjust = 0.15) +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +#定义方块颜色
      scale_size_manual(values = c(2, 1, 0.5)) +
      scale_colour_manual(values = p.color) +
      guides(size = guide_legend(title = paste0(test.method, "'s p"),#定义图例
                                 override.aes = list(colour = "grey35"),
                                 order = 2),
             colour = guide_legend(title = paste0(test.method, "'s r"),
                                   override.aes = list(size = 3),
                                   order = 1),
             fill = guide_colorbar(title = paste0(cor.method,"'s r"), order = 3))

  }else{
    # # mantel test could perform between multiple x (spe) and single y
    # multiple x col is a spec
    library(dplyr)

    mantel <- linkET::mantel_test(x.axis.df, y.axis.df,
                                  spec_select  = spec.select.list) %>%
      mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                      labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
             pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                      labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
    test.method = "Mantel"

    # # # # # # # # # # # # # # # # # # # # # # # #
    plot = qcorrplot(correlate(y.axis.df, method = cor.method), type = "lower", diag = FALSE) +
      geom_square() +
      geom_couple(aes(colour = pd,
                      size = rd),
                  data = mantel, curvature = 0.1) +
      #geom_diag_label(mapping = aes(y = .y + 0.05),#定义对角线上的文字
      #                hjust = 0.15) +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +#定义方块颜色
      scale_size_manual(values = c(0.5, 1, 2)) +
      scale_colour_manual(values = p.color) +
      guides(size = guide_legend(title = paste0(test.method, "'s r"),#定义图例
                                 override.aes = list(colour = "grey35"),
                                 order = 2),
             colour = guide_legend(title = paste0(test.method, "'s p"),
                                   override.aes = list(size = 3),
                                   order = 1),
             fill = guide_colorbar(title = paste0(cor.method,"'s r"), order = 3))
    # + theme(axis.text.y = element_blank())


  }


  res = list(plot = plot, mantel = mantel)
  res
}







#' Immune related gene
#'
#' @return A tibble object
#' @export
#'
#' @examples loonR::return_immune_genes()
return_immune_genes <- function(type=c("default", "checkpoint") ){

  if(length(type)>1){
     type = "default"
  } else{
    type = match.arg(type)
  }


  immune.related.genes = tibble::tribble(
    ~HGNC.Symbol,        ~Super.Category, ~Entrez.ID,    ~Symbol, ~Immune.Checkpoint,
    "ADORA2A",             "Receptor",       135L,  "ADORA2A",       "Inhibitory",
    "ARG1",                "Other",       383L,     "ARG1",       "Inhibitory",
    "BTLA",             "Receptor",    151888L,     "BTLA",       "Inhibitory",
    "BTN3A1",         "Co-inhibitor",     11119L,   "BTN3A1",      "Stimulatory",
    "BTN3A2",         "Co-inhibitor",     11118L,   "BTN3A2",      "Stimulatory",
    "CCL5",               "Ligand",      6352L,     "CCL5",      "Stimulatory",
    "CD27",             "Receptor",       939L,     "CD27",      "Stimulatory",
    "CD274",         "Co-inhibitor",     29126L,    "CD274",       "Inhibitory",
    "CD276",         "Co-inhibitor",     80381L,    "CD276",       "Inhibitory",
    "CD28",        "Co-stimulator",       940L,     "CD28",      "Stimulatory",
    "CD40",             "Receptor",       958L,     "CD40",      "Stimulatory",
    "CD40LG",               "Ligand",       959L,   "CD40LG",      "Stimulatory",
    "CD70",               "Ligand",       970L,     "CD70",      "Stimulatory",
    "CD80",        "Co-stimulator",       941L,     "CD80",      "Stimulatory",
    "CTLA4",             "Receptor",      1493L,    "CTLA4",       "Inhibitory",
    "CX3CL1",               "Ligand",      6376L,   "CX3CL1",      "Stimulatory",
    "CXCL10",               "Ligand",      3627L,   "CXCL10",      "Stimulatory",
    "CXCL9",               "Ligand",      4283L,    "CXCL9",      "Stimulatory",
    "EDNRB",             "Receptor",      1910L,    "EDNRB",       "Inhibitory",
    "ENTPD1",                "Other",       953L,   "ENTPD1",      "Stimulatory",
    "GZMA",                "Other",      3001L,     "GZMA",      "Stimulatory",
    "HAVCR2",             "Receptor",     84868L,   "HAVCR2",       "Inhibitory",
    "HLA-A", "Antigen presentation",      3105L,    "HLA-A",                 "NA",
    "HLA-B", "Antigen presentation",      3106L,    "HLA-B",                 "NA",
    "HLA-C", "Antigen presentation",      3107L,    "HLA-C",                 "NA",
    "HLA-DPA1", "Antigen presentation",      3113L, "HLA-DPA1",                 "NA",
    "HLA-DPB1", "Antigen presentation",      3115L, "HLA-DPB1",                 "NA",
    "HLA-DQA1", "Antigen presentation",      3117L, "HLA-DQA1",                 "NA",
    "HLA-DQA2", "Antigen presentation",      3118L, "HLA-DQA2",                 "NA",
    "HLA-DQB1", "Antigen presentation",      3119L, "HLA-DQB1",                 "NA",
    "HLA-DQB2", "Antigen presentation",      3120L, "HLA-DQB2",                 "NA",
    "HLA-DRA", "Antigen presentation",      3122L,  "HLA-DRA",                 "NA",
    "HLA-DRB1", "Antigen presentation",      3123L, "HLA-DRB1",                 "NA",
    "HLA-DRB3", "Antigen presentation",      3125L, "HLA-DRB3",                 "NA",
    "HLA-DRB4", "Antigen presentation",      3126L, "HLA-DRB4",                 "NA",
    "HLA-DRB5", "Antigen presentation",      3127L, "HLA-DRB5",                 "NA",
    "HMGB1",                "Other",      3146L,    "HMGB1",      "Stimulatory",
    "ICAM1",        "Cell adhesion",      3383L,    "ICAM1",      "Stimulatory",
    "ICOS",             "Receptor",     29851L,     "ICOS",      "Stimulatory",
    "ICOSLG",        "Co-stimulator",     23308L,   "ICOSLG",      "Stimulatory",
    "IDO1",                "Other",      3620L,     "IDO1",       "Inhibitory",
    "IFNA1",               "Ligand",      3439L,    "IFNA1",      "Stimulatory",
    "IFNA2",               "Ligand",      3440L,    "IFNA2",      "Stimulatory",
    "IFNG",               "Ligand",      3458L,     "IFNG",      "Stimulatory",
    "IL10",               "Ligand",      3586L,     "IL10",       "Inhibitory",
    "IL12A",               "Ligand",      3592L,    "IL12A",      "Stimulatory",
    "IL13",               "Ligand",      3596L,     "IL13",       "Inhibitory",
    "IL1A",               "Ligand",      3552L,     "IL1A",      "Stimulatory",
    "IL1B",               "Ligand",      3553L,     "IL1B",      "Stimulatory",
    "IL2",               "Ligand",      3558L,      "IL2",      "Stimulatory",
    "IL2RA",             "Receptor",      3559L,    "IL2RA",      "Stimulatory",
    "IL4",               "Ligand",      3565L,      "IL4",       "Inhibitory",
    "ITGB2",        "Cell adhesion",      3689L,    "ITGB2",      "Stimulatory",
    "KIR2DL1",             "Receptor",      3802L,  "KIR2DL1",       "Inhibitory",
    "KIR2DL2",             "Receptor",      3803L,  "KIR2DL2",       "Inhibitory",
    "KIR2DL3",             "Receptor",      3804L,  "KIR2DL3",       "Inhibitory",
    "LAG3",             "Receptor",      3902L,     "LAG3",       "Inhibitory",
    "MICA", "Antigen presentation", 100507436L,     "MICA",                 "NA",
    "MICB", "Antigen presentation",      4277L,     "MICB",                 "NA",
    "PDCD1",             "Receptor",      5133L,    "PDCD1",       "Inhibitory",
    "PDCD1LG2",         "Co-inhibitor",     80380L, "PDCD1LG2",                 "NA",
    "PRF1",                "Other",      5551L,     "PRF1",      "Stimulatory",
    "SELP",        "Cell adhesion",      6403L,     "SELP",      "Stimulatory",
    "SLAMF7",         "Co-inhibitor",     57823L,   "SLAMF7",       "Inhibitory",
    "TGFB1",               "Ligand",      7040L,    "TGFB1",       "Inhibitory",
    "TIGIT",             "Receptor",    201633L,    "TIGIT",       "Inhibitory",
    "TLR4",             "Receptor",      7099L,     "TLR4",      "Stimulatory",
    "TNF",               "Ligand",      7124L,      "TNF",      "Stimulatory",
    "TNFRSF14",             "Receptor",      8764L, "TNFRSF14",      "Stimulatory",
    "TNFRSF18",             "Receptor",      8784L, "TNFRSF18",      "Stimulatory",
    "TNFRSF4",             "Receptor",      7293L,  "TNFRSF4",      "Stimulatory",
    "TNFRSF9",             "Receptor",      3604L,  "TNFRSF9",      "Stimulatory",
    "TNFSF4",               "Ligand",      7292L,   "TNFSF4",      "Stimulatory",
    "TNFSF9",               "Ligand",      8744L,   "TNFSF9",      "Stimulatory",
    "VEGFA",               "Ligand",      7422L,    "VEGFA",       "Inhibitory",
    "VEGFB",               "Ligand",      7423L,    "VEGFB",       "Inhibitory",
    "C10orf54",         "Co-inhibitor",     64115L, "C10orf54",       "Inhibitory",
    "VTCN1",         "Co-inhibitor",     79679L,    "VTCN1",       "Inhibitory"
  )




  message("Collected from: ")
  if(type=="default"){
    message("https://www.cell.com/immunity/fulltext/S1074-7613(18)30121-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761318301213%3Fshowall%3Dtrue#supplementaryMaterial")
    res = immune.related.genes
  }else if(type=="checkpoint"){
    message("https://www.sciencedirect.com/science/article/pii/S1936523319303389#:~:text=These%20include%20the%20coinhibitory%20molecules,and%20VISTA%20and%20CD277%2C%20the")
    message("Table 1. A List of All Potential Checkpoints on the Cancer/APC Side if the Immunological Synapse.")

    message("https://academic.oup.com/bib/article/22/3/bbaa176/5894466#supplementary-data")
    message("In total, 79 ICGs were identified from a review of the literature [1, 25, 26], and most of the ICGs were ligands, receptors or important molecules in immune checkpoint pathways (Supplementary Table S1). W")

  }


  res

}
