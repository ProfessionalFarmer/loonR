

#' Export figure to PPT
#'
#' @param obj ggplot object
#' @param file default is "~/test.pptx"
#' @param append 
#'
#' @return NA
#' @export
#'
#' @examples
#' 
export2ppt <- function(obj,file="~/test.pptx", append=TRUE){
  library(export)
  graph2ppt(obj, file=file, append=append)
}