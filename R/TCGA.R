
#' Get sample barcode and TCGA project map information
#'
#'
#' @param project A vector.
#'  TCGA-ACC, TCGA-COAD.
#'  See getGDCprojects()$project_id
#' @return A data.frame
#'  X1 proj samples
#'	1	1	TCGA-ACC	TCGA-OR-A5JZ
#'	2	1	TCGA-ACC	TCGA-OR-A5JK
#' @export
#'
#' @examples
get.sample.project.infomap <- function(project=NA){

  library(TCGAbiolinks)

  # If not specify project, load default
  if (project==NA){
    project <- grep("TCGA", sort(getGDCprojects()$project_id),
                    value = TRUE)
  }


  df <- plyr::adply(project, .margins = 1, .fun = function(proj) {
    samples <- TCGAbiolinks:::getSubmitterID(proj)
    return(data.frame(proj, samples))
  })

  return(df)

}
