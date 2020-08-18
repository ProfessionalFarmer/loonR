
#' Title
#'
#' @param df data.frame. Col is sample, row is gene
#' @param label
#' @param k Folds
#'
#' @return
#'
#' @examples
getOneRoundCVRes <- function(df, label, k){

  set.seed(666+i)
  require(caret)
  flds <- createFolds(label, k = k, list = FALSE, returnTrain = FALSE)

  # Every fold will be used in one round.

  res <- foreach::foreach(i= 1:k, .combine = rbind) %do% {

    piece.ind <- which(flds==i)

    lg.df <- data.frame(label= factor( label[-piece.ind], levels = c(FALSE, TRUE), labels = c(0,1) ),
                        df[-piece.ind, ])

    #colnames(lg.df) <- gsub("\\.", "-", colnames(lg.df) )

    set.seed(666)

    # The type="response" option tells R to output probabilities of the form P(Y = 1|X), as opposed to other information such as the logit.
    glm.fit <- glm(label ~ ., data = lg.df, family = binomial(logit))

    pre.res <- predict(glm.fit, df[piece.ind, ])

    data.frame(Name  = names(pre.res),
               Logit = pre.res,
               Label = label[piece.ind],
               stringsAsFactors = FALSE)
  }

  res

}


#' Title
#'
#' @param df data.frame. Col is sample, row is gene
#' @param label
#' @param k Folds
#' @param n Repeat times
#'
#' @return
#' @export
#'
#' @examples
cross.validation <- function(df = '', label = '', k = 5, n = 100){

  # df = up.mirs.exp
  # label = mir.sample.group
  # k = 5
  # n = 100

  if(df == '' | label == ''){
    stop("provide a datafram or lables")
  }

  library(doParallel)
  library(foreach)
  cl = parallel::makeCluster(40)
  doParallel::registerDoParallel(cl)

  res <- foreach::foreach(i= 1:n, .combine = rbind, .packages="foreach", .export=c("getOneRoundCVRes")) %dopar% {
  #res <- foreach::foreach(i= 1:n, .combine = rbind) %do% {

    getOneRoundCVRes(df, label, k)

  }

  cv.res.mean <- cv.res %>% group_by(Name, Label) %>% summarize(Mean = mean(Logit))

  stopCluster(cl)

  cv.res.mean
}















