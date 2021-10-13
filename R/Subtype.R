#' Calculate Jaccard similarity coefficient
#'
#' @param df row is sample, column is study subtype
#' @param returnDf reture a n*n matrix
#' @param concateStudy
#'
#' @return
#' @export
#'
#' @examples
consensusNodeJaccard <- function(df, returnDf = T, concateStudy=F){

  study.Com <- loonR::generateCombinations(colnames(df), size=2, repeats = T)

  library(foreach)

  res <- foreach(i=1:nrow(study.Com), .combine = rbind) %do% {
      study1.name = study.Com[i,1]
      study2.name = study.Com[i,2]

      study1.vector = as.vector(unlist(df[study1.name]))
      study2.vector = as.vector(unlist(df[study2.name]))

      foreach(subtype1 = unique(study1.vector), .combine = rbind) %do%{

        foreach(subtype2 = unique(study2.vector), .combine = rbind) %do%{

          ind = study1.vector == 0 & study2.vector == 0

          j.coef = jaccard::jaccard(study1.vector[ind]==subtype1, study2.vector[ind]==subtype2)

          if(concateStudy){
            c( paste0(study1.name,"-",subtype1), paste0(study2.name, "-", subtype2), j.coef)
          }else{
            c( subtype1, subtype2, j.coef )
          }

        }

      }

  }
  colnames(res) <- c("Subtype1","Subtype2","Jaccard.index")
  res <- data.frame(res, check.names = F)


  if(returnDf){
    res = reshape2::dcast(res, Subtype1~Subtype2, value.var = 'Jaccard.index')
    rownames(res) <- res$Subtype1
    res <- data.frame( res[,-c(1)], check.names = F )
    # reorder
    res <- res[colnames(res),]
    loonR::fillSymmetricNATable(res)

  }else{
    res
  }

}


#' Fill NA in symmetric table
#'
#' @param symmetric.df
#'
#' @return
#' @export
#'
#' @examples
fillSymmetricNATable <- function(symmetric.df){

    if(sum(colnames(symmetric.df)!=rownames(symmetric.df))!=0){
      stop("Row name is not the same with colname")
    }

    for(r in 1:nrow(symmetric.df)){
      for(c in 1:ncol(symmetric.df)){
         value = symmetric.df[r,c]
         rev_value = symmetric.df[c,r]

         if(is.na(value) & !is.na(rev_value)){
           symmetric.df[r,c] = rev_value
         }
         if(!is.na(value) & is.na(rev_value)){
           symmetric.df[c,r] = value
         }

      }

    }
  symmetric.df

}




