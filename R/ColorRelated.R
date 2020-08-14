#' Generate a most distinctive color palatte
#'
#' @param n Number of colors
#'
#' @return Color vector
#' @export
#'
#' @examples
get.mostDistint.color.palette <- function(n=20){
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  set.seed(123)
  col = sample(col_vector, n)
  return(col)
}



#' Get colors by ggsci
#'
#' @param palette default: nrc. Palette ref: https://nanx.me/ggsci/articles/ggsci.html
#' @param alpha 1
#' @param n number of colors, default n = 10
#'
#' @return
#' @export
#'
#' @examples
get.ggsci.color <- function(palette="nrc", alpha=1, n = 10){

  library(ggsci)
  f <- parse(text=paste("pal_", palette, sep = "")  )
  myPalette <- eval(f)(alpha = alpha)(n)
  myPalette

}






