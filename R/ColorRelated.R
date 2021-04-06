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
#' Allowed values include "grey" for grey color palettes;
#' Rbrewer palettes e.g. "RdBu", "Blues", ...;
#' Custom color palette e.g. c("blue", "red");
#' Scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#'     https://github.com/nanxstats/ggsci
#' Wes Anderson color palettes: "BottleRocket1"  "BottleRocket2"  "Rushmore1" "Rushmore" "Royal1" "Royal2" "Zissou1" "Darjeeling1" "Darjeeling2" "Chevalier1" "FantasticFox1" "Moonrise1" "Moonrise2" "Moonrise3" "Cavalcanti1" "GrandBudapest1" "GrandBudapest2" "IsleofDogs1" "IsleofDogs2"
#'     https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
#' Color blind palette: cbPalette, cbbPalette
#'
#' @param alpha default 1
#' @param n number of colors, default n = 7
#'
#' @return
#' @export
#'
#' @examples get.palette.color("nrc", n = 2, alpha = 0.7)
#' scales::show_col(colorBlindGrey8)
get.palette.color <- function(palette="nrc", n = 7, alpha=1){


  my_palettes <- list(
    # From Gfplot
    `jama_classic` = c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"),
    # The palette with grey: color blind
    `cbPalette` = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    # The palette with black: color blind
    `cbbPalette` = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  )


  if(palette %in% names(my_palettes) ){
    colors = head(my_palettes[[c(palette)]],n)

  }else if (palette %in% c(names(wesanderson::wes_palettes)) ) {
    colors = head(wesanderson::wes_palettes[[c(palette)]],n, type = c("discrete"))

  }else if (palette == "Most"){
    colors = loonR::get.mostDistint.color.palette(n)

  }else{
     colors = ggpubr::get_palette(palette = palette, k = n)
  }

  # set alpha
  scales::alpha(colors, alpha)
  # scales::show_col(colorBlindGrey8)




}






