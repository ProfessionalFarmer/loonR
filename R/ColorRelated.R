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
#' @param n number of colors, default n = 7
#'
#' @return
#' @export
#'
#' @examples get.ggsci.color("nrc", n = 2, alpha = 0.7)
get.ggsci.color <- function(palette="nrc", n = 7, alpha=1){

  library(ggsci)
  f <- parse(text=paste("pal_", palette, sep = "")  )
  myPalette <- eval(f)(alpha = alpha)(n)
  myPalette

}


#' Use ggpubr function to get all kinds of color palette
#' RColorBrewer, ggsci, hue or grey/gray color palettes
#'
#' @param palette Default npg. The color palette to be used for coloring or filling by groups.
#' Allowed values include "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...;
#' or custom color palette e.g. c("blue", "red");
#' and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#' @param alpha default 1
#' @param n number of colors, default n = 7
#'
#' @return
#' @export
#'
#' @examples get.palette.color("nrc", n = 2, alpha = 0.7)
get.palette.color <- function(palette="nrc", n = 7, alpha=1){

  if (palette == "jama_classic"){
    colors = head(c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9",
                    "#CAC27E", "#818282"), n)
  }else{
     colors = ggpubr::get_palette(palette = palette, k = n)
  }

  # set alpha
  scales::alpha(colors, alpha)

}






