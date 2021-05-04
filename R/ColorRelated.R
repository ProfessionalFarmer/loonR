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
#'
#' RColorBrewer, ggsci, hue or grey/gray color palettes
#'
#' @param palette Default npg
#' @param alpha default 1
#' @param n number of colors, default n = 7
#' @param install Default FALSE. TURE if Install related package
#' @param show.color Default FALSE
#'
#' @return
#' @export
#'
#' @examples get.palette.color("npg", n = 2, alpha = 0.7)
#' scales::show_col(colorBlindGrey8)
#'
#'
#' Default npg. The color palette to be used for coloring or filling by groups.
#' Allowed values include "grey" for grey color palettes;
#'
#' Rbrewer palettes e.g. "RdBu", "Blues", ...;
#' Custom color palette e.g. c("blue", "red");
#'
#' Scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#'     https://github.com/nanxstats/ggsci
#'
#' Wes Anderson color palettes: https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
#'     "BottleRocket1"  "BottleRocket2"  "Rushmore1" "Rushmore" "Royal1" "Royal2" "Zissou1" "Darjeeling1" "Darjeeling2" "Chevalier1" "FantasticFox1" "Moonrise1" "Moonrise2" "Moonrise3" "Cavalcanti1" "GrandBudapest1" "GrandBudapest2" "IsleofDogs1" "IsleofDogs2"
#'
#' LaCroix Color Palettes for R: https://github.com/johannesbjork/LaCroixColoR
#'     PassionFruit, Mango, Pure, Lime, Lemon, Orange, Berry, CranRaspberry, Pamplemousse, PeachPear, Coconut, Apricot, Tangerine, KeyLime, PommeBaya, CeriseLimon, PinaFraise, KiwiSandia, MelonPomelo, MurePepino, paired
#'
#' ggthemes::canva_palettes 150 Color Palettes from Canva
#' https://jrnold.github.io/ggthemes/reference/canva_palettes.html
#' Fresh and bright, Subdued and proffesional, Dark and earthy, Crisp and dramatic, Cool blues, Outdoorsy and Natural, Watery blue-greens, Primary colors with a vibrant twist, Refreshing and pretty, Playful greens and blues, Fresh and energetic, Surf and turf, Autumn in vermont, Icy blues and grays, Birds and berries, Day and night, Stylish and retro, Shades of citrus, Sunset to dusk, Bright and tropical, Warm naturals, Bold berries, Summer sunflower, Modern and crisp, Timeless and nautical, Neutral and versatile, Cheerful brights, Garden fresh, Summer barbeque, Berry blues, Lemonade stand, Serene and spa like, Fun and tropical, Spicy neutrals, Pastels, Bold and cultured, Sunny citrus, Crisp complementary colors, Warm and rustic, Neon night, Jewel tones, Polished and inviting, Fresh greens, Wintery reds, Summer fiesta, Chocolaty browns, Naturally elegant, Cozy and warm, Violet sunset, Strawberries and cream, Grecian holiday, Bold and basic, Vineyard neutrals, Modern and urban, Misty greens, Sunkissed village, Sun and sky, Aqua blues, Urban oasis, Candy coated brights, Muted and antique, Classy and timeless, Cosmopolitan, Cheerful and friendly, Nightlife, Coastal, Maritime brights, Vintage charm, Understated and versatile, Artic sunrise, Mediterranean afternoon, Hazy grays, City sights, Retro and relaxing, Green fields, Distintive and unexpected, Sleek and modern, Orange accent, Beyond black and white, Shabby chic neutrals, Warm and cool, Industrial and in control, Autumn oranges and complemtentary neutrals, Pool party, Classic metallics, Subtle and versatile, Professional and traditional, Light and natural, Shadowy and dramatic, Golden afternoon, Dark and handsome, Technology meets nature, Cheerful blues and pink, Exotic and high impact, Back to school, Bright and painterly, Urban living, 1950s kitchen, Smoky purples, Trendy and metropolitan, Fun and professional, Art history inspired, Muted tones, Modern and clean, Neon tones and sharp contrast, Muted and minimal, Warm and bold, Clean and highlighted, Warm tones, Sharp and modern, Cool vs warm, Pretty pastels, Bold and punchy, Tints and tones, Splash of color, Elegant and sophisticated, Summer inspired, Professional and modern, Bold blacks and vibrand highlights, Clean gradient and fresh blues, Cheerful and sleek, Luxurious and modern, Unique and striking, Unexpected color combinations, Retro inspired, Antique and clean, Striking and energetic, Fresh and lively, Clean and crisp, Colorful without clashing, Cool and calm, Modern and muted, Earthy and fresh, High saturation and high energy, Warm and wonderful, Vintage charm 2, Cool jewel tones, Stormy hues, Clean and collegiate, Simple and fresh, Tropical tones, Bold feature colors, Antique tones, Neon and bold, Simple but bold, Corporate and sleek, Modern and minimal, Fun and cheerful, Sunny and calm, Pop art
#'
#' ggthemes color: https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/
#'
#'
get.palette.color <- function(palette="npg", n = 7, alpha=1, install=FALSE, show.color=FALSE){

  if(install){
    # wesanderson
    install.packages("wesanderson")
    # ggsci
    BiocManager::install("ggsci")
    # LaCroix Color Palettes for R
    devtools::install_github("johannesbjork/LaCroixColoR")
    # ggthemes
    devtools::install_github(c("hadley/ggplot2", "jrnold/ggthemes"))
    # ggpubr
    BiocManager::install("ggpubr")
    return()
  }


  my_palettes <- list(
    # From Gfplot
    `jama_classic` = c("#164870", "#10B4F3", "#FAA935", "#2D292A", "#87AAB9", "#CAC27E", "#818282"),
    # The palette with grey: color blind
    `cbPalette` = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    # The palette with black: color blind
    `cbbPalette` = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
    # selected from public
    `sp1` = c("#00AFBB", "#E7B800", "#0392cf", "#7570B3", "#FC4E07", "#BB3099", "#ADC252", "#be9b7b", "#75A3BA", "#bbbbbb", "#4F7175"  )
  )


  if(palette %in% names(my_palettes) ){
    colors = head(my_palettes[[c(palette)]],n)

  }else if (palette %in% c(names(wesanderson::wes_palettes)) ) {
    colors = head(wesanderson::wes_palettes[[c(palette)]],type = c("discrete"), n)

  }else if (palette == "Most"){
    colors = loonR::get.mostDistint.color.palette(n)

  }else if (palette %in% names(LaCroixColoR::lacroix_palettes)){
    colors = head(LaCroixColoR::lacroix_palette(palette, type = "discrete"), n)

  }else if (palette %in% names(ggthemes::canva_palettes)){
    colors = head(ggthemes::canva_palettes[[palette]], n)

  }else{# include ggsci
     colors = ggpubr::get_palette(palette = palette, k = n)
  }



  # set alpha
  col.pal <- scales::alpha(colors, alpha)
  # scales::show_col(colorBlindGrey8)
  if(show.color){
    scales::show_col(col.pal)
  }else{
    col.pal
  }

}






