#' pairwiseAdonis
#'
#' @param data Row. is sample, column is feature
#' @param group
#'
#' @return
#' @export
#'
#' @examples
#' data(iris)
#' loonR::pairwiseAdonisTest(iris[,1:4], iris$Species)
pairwiseAdonisTest <- function(data, group){
  # Adonis，多元方差分析，亦可称为非参数多元方差分析。其原理是利用距离矩阵（比如基于Bray-Curtis距离、Euclidean距离）对总方差进行分解，分析不同分组因素对样品差异的解释度，并使用置换检验对其统计学意义进行显著性分析。
  if(!require(pairwiseAdonis))  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  # https://github.com/pmartinezarbizu/pairwiseAdonis
  # https://www.jianshu.com/p/dfa689f7cafd

  pairwise.adonis(data, group)

}
