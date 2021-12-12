#' Iterate all the combination of housekeeping genes
#'
#' @param raw.ct.df row is sample
#' @param housekeeping.ct.df row is sample
#' @param group Second shuld be TRUE
#' @param alternative Default "less". One of c("less", "greater", "both")
#' @param p.cutoff Default 0.05
#' @param difference.cutoff Default 0. Should be greater than 0. The direction of difference depends on alternative.
#' @param sep Default ", "
#' @param logTen log10(2^(-ΔCT))
#' @param expTwo 2^(-ΔCT)
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
checkHousekeepingCombination <- function(raw.ct.df, housekeeping.ct.df, group, alternative = "less", p.cutoff = 0.05, difference.cutoff = 0, sep = ", ", logTen=FALSE, expTwo=FALSE, scale=TRUE){

  if(difference.cutoff<0){
    stop("Pls set difference.cutoff greater or equal to 0. The direction depends on alternative")
  }


  library(foreach)
  library(dplyr)

  normalized.ct.res = list()
  differential.res = list()

  alternative <- match.arg(alternative, choices = c("less", "greater", "both"), several.ok = FALSE)


  iteration.res <- foreach(comb.size=1:ncol(housekeeping.ct.df), .combine = rbind) %do% {

    comb.df <- loonR::generateCombinations(colnames(housekeeping.ct.df), comb.size, vector = F)

    all.comb.res <- foreach(row.ind = 1:nrow(comb.df), .combine = rbind)%do% {

      hk.combi = unlist(comb.df[row.ind,])

      if(length(hk.combi)!=1){
        hk.ave.value = rowMeans( housekeeping.ct.df[, hk.combi] )
      }else{
        hk.ave.value = housekeeping.ct.df[,unclass(hk.combi) ]
      }

      hk.ave.value = as.numeric(hk.ave.value)

      normalized.df <- sweep(raw.ct.df, 1, hk.ave.value, FUN = "-")

      # if log10(2^(-ΔCT)) or 2^(-ΔCT)
      if(logTen){
        normalized.df = log10( 2^(normalized.df*(-1) )  )
      }
      if(expTwo){
        normalized.df = 2^(normalized.df*(-1) )
      }

      diff.res <- loonR::ttest_differential(t(normalized.df), group, cores = 1, alternative = alternative)

      # store normalized value
      normalized.ct.res[[paste0( comb.df[row.ind,], collapse = sep)]] = normalized.df
      differential.res[[paste0( comb.df[row.ind,], collapse = sep)]] = diff.res

      if(alternative=="less"){
        candidates.res <- diff.res %>% filter(P < p.cutoff & Difference < difference.cutoff)
      }else if (alternative=="greater"){
        candidates.res <- diff.res %>% filter(P < p.cutoff & Difference > difference.cutoff)
      }else if (alternative=="both"){
        candidates.res <- diff.res %>% filter(P < p.cutoff & abs(Difference) > difference.cutoff)
      }

      # cat("Combination size ", comb.size, ": ", hk.combi)
      # cat("# of candidates: ", nrow(candidates.res) , "\n" )

      if(nrow(candidates.res)==0){
        lg.res = list()
        lg.res$AUC = 0.5
        lg.res$eliminationCandidates=NULL
      }else{
        lg.res <- loonR::build.logistic.model(normalized.df[,candidates.res$Name], group, scale = scale)
      }

      # backward eliminatio
      # lg.res <- loonR::build.logistic.model(normalized.df[,lg.res$eliminationCandidates], group)

      # recursive feature selection
      #rfe.res <- loonR::feature.selection.RFE(normalized.df[,candidates.res$Name], group)

      # Lasso selection
      # lasso.res <- loonR::lasso.select.feature(normalized.df[,candidates.res$Name], group)
      # lg.res <- loonR::build.logistic.model(normalized.df[,lasso.res$candidates], group)

      auc = lg.res$AUC

      #cat("AUC is ", auc, "\n\n\n")

      comb.res = c(comb.size, auc, nrow(candidates.res),
            paste0( comb.df[row.ind,], collapse = sep),
            length(lg.res$eliminationCandidates)  )

      comb.res = c(comb.res, paste0( candidates.res$Name, collapse = sep) ) # candidates
      comb.res = c(comb.res, paste0( lg.res$eliminationCandidates, collapse = sep) ) # elimination

      names(comb.res) <- c("Combination Size", "AUC", "# of Candidates", "Housekeeping", "# of elimination candi", "Candidate miR", "Elimination miR")

      comb.res
    }

    all.comb.res
  }
  iteration.res = data.frame(iteration.res, check.names = F)
  iteration.res = iteration.res[order(iteration.res$AUC, decreasing = T),]

  res = list(
    Result.table = iteration.res,
    Normalized.CT = normalized.ct.res,
    Differential.Res = differential.res
  )

  warning("Plsease note the hypothesis is : ", alternative, "\n")
  warning("Plsease note the p.cutoff is : ", p.cutoff, "\n")
  warning("Plsease note the difference.cutoff is : ", difference.cutoff, "\n---------\n")
  res

}


