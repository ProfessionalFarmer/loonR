

#' Run estiamte
#'
#' @param expression.df   the row name in an input data must be gene symbols
#' @param platform c("affymetrix", "agilent", "illumina")
#'
#' @return
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
#' @return
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


