
# This script mainly refered the following link
# https://www.machinelearningplus.com/machine-learning/feature-selection/#2variableimportancefrommachinelearningalgorithms


#' Feature selection by Boruta
#'
#' @param feature.df Row is sample, column is feature
#' @param group
#' @param seed Default 111
#' @param scale
#' @param withTentative Default FALSE, not include tentative variable
#'
#' @return
#' @export
#'
#' @examples loonR::feature.selection.boruta(miR.df, group)
feature.selection.boruta <- function(feature.df, group, seed=111, scale=TRUE, withTentative = FALSE){

  library(Boruta)

  trainData = feature.df
  if(scale){
    trainData <- loonR::scaleDF(trainData, byColumn = T)
  }

  trainData <- data.frame(Class=group,
                        trainData,
                        check.names = F)

  trainData$Class <- factor(trainData$Class,
                          levels = unique(trainData$Class),
                          labels = 0:(length(unique(trainData$Class))-1)  )


  # Perform Boruta search
  set.seed(seed)
  boruta_output <- Boruta(Class ~ ., data=trainData, doTrace=0)
  ## Do a tentative rough fix
  #boruta_output <- TentativeRoughFix(boruta_output)
  colnames(boruta_output$ImpHistory) <- stringr::str_remove_all(colnames(boruta_output$ImpHistory), "\`")

  boruta_signif <- getSelectedAttributes(boruta_output, withTentative = withTentative)

  res = list(boruta.res = boruta_output)
  res$candidates <- stringr::str_remove_all(boruta_signif,"`")


  # Variable Importance Scores
  imps <- attStats(boruta_output)
  row.names(imps) <- stringr::str_remove_all(row.names(imps),"`")
  imps = imps[order(-imps$meanImp), ] # descending sort
  res$importance <- imps


  imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
  imps2 = imps2[order(-imps2$meanImp), ] # descending sort
  res$importance.clean = imps2


  plot.obj <- plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

  res$plot <- plot.obj

  res

}




#' Variable Importance from a specific Machine Learning Algorithms
#'
#'
#' @param feature.df Row is sample, column is feature
#' @param group
#' @param seed
#' @param method glm for logistic, rf for random forest, RRF for Regularized Random Forest, rpart for rpart model. https://www.rdocumentation.org/packages/caret/versions/4.47/topics/train
#' @param scale
#'
#' @return
#' @export
#'
#' @examples loonR::feature.selection.4givenMLAlgorithm(miR.df, group, method = 'glm')
feature.selection.4givenMLAlgorithm <- function(feature.df, group, seed=111, scale=TRUE, method = 'glm'){

  library(caret)


  trainData = feature.df
  if(scale){
    trainData <- loonR::scaleDF(trainData, byColumn = T)
  }

  trainData <- data.frame(Class=group,
                          trainData,
                          check.names = F)

  trainData$Class <- factor(trainData$Class,
                            levels = unique(trainData$Class),
                            labels = 0:(length(unique(trainData$Class))-1)  )


  # Perform
  set.seed(seed)

  Mod <- train(Class ~ ., data=trainData, method=method)
  imps <- varImp(Mod)
  rownames(imps$importance) <- stringr::str_remove_all(rownames(imps$importance), "\\\\")
  rownames(imps$importance) <- stringr::str_remove_all(rownames(imps$importance), "`")


  plot.obj <- plot(imps, main='Variable Importance')


  res <- list(importance=data.frame(imps$importance))
  res$model = Mod
  res$method = imps$model
  res$plot = plot.obj

  res
}



#' Feature selection by lasso
#'
#' @param data.matrix Row is sample
#' @param label
#' @param folds
#' @param seed
#' @param family Default binomial. Should be one of “gaussian”, “binomial”, “poisson”, “multinomial”, “cox”, “mgaussian”
#' @param type.measure Default auc. Can be class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.
#' @param s Defalut is lambda.min. User can specify
#' @param scale Default TRUE
#'
#' @return
#' @export
#'
#' @examples
lasso.select.feature <- function(data.matrix, label, folds = 5, seed = 666,
                                 family = "binomial", type.measure = "auc" ,
                                 s = NULL, scale=TRUE){

  library(foreach)
  library(dplyr)
  library(glmnet)


  if(scale){
    data.matrix = scale(data.matrix, center = TRUE, scale = TRUE)
    #data.matrix = data.frame(data.matrix, check.names = T, stringsAsFactors = F)
  }


  set.seed(seed)
  require(caret)
  cvfit = cv.glmnet(as.matrix(data.matrix), label, nfolds = folds,
                    family = family, type.measure = type.measure)

  if (is.null(s)) {
    s = cvfit$lambda.min
  }else{
    s = s
  }

  plot.obj <- plot(cvfit)


  feature.coef = coef(cvfit, s = s)
  feature.coef = data.frame(name = feature.coef@Dimnames[[1]][feature.coef@i + 1], coefficient = feature.coef@x)

  feature.coef = feature.coef[-c(1), ] # remove Intercept
  feature.coef$auc = apply(data.matrix[,feature.coef$name], 2, function(x){
    suppressMessages(roc <- pROC::roc(label, x)  )
    roc$auc
  })
  row.names(feature.coef) <- feature.coef$name


  res = list(model = cvfit,
             coefficient = feature.coef,
             candidates = row.names(feature.coef),
             plot = plot.obj
             )
  res

}


#' Perform multiple round lasso to select stable feature
#'
#' @param data.matrix Row is sample
#' @param label
#' @param folds
#' @param seed
#' @param family Default binomial
#' @param n 100
#' @param cores 50
#' @param type.measure Default auc. Can be class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.
#' @param scale Default TRUE
#'
#' @return
#' @export
#'
#' @examples
lasso.cv.select.feature <- function(data.matrix, label, folds = 5, seed = 666, n = 100,
                                    family = "binomial", type.measure = "auc" ,
                                    cores = 50, scale=TRUE){
  library(foreach)
  library(dplyr)
  library(parallel, doParallel)
  library(glmnet)

  if(scale){
    data.matrix = scale(data.matrix, center = TRUE, scale = TRUE)
    data.matrix = data.frame(data.matrix, check.names = F, stringsAsFactors = F)
  }

  set.seed(seed)
  require(caret)
  doParallel::registerDoParallel(cores=cores)
  parallel::mcaffinity(c(1:cores)) # limit cores to use

  res.raw <- foreach::foreach(i=1:n, .combine = rbind, .packages = c("dplyr")) %dopar% {

    set.seed(seed+i)
    flds <- createFolds(label, k = folds, list = FALSE, returnTrain = FALSE)
    ind = !(flds == 1)

    cvfit = cv.glmnet(as.matrix(data.matrix[ind,]), label[ind], nfolds = folds-1, # cross validatin fold - 1
                      family = family, type.measure = type.measure)


    feature.coef = coef(cvfit, s=cvfit$lambda.min)
    feature.coef = data.frame(name = feature.coef@Dimnames[[1]][feature.coef@i + 1], coefficient = feature.coef@x)

    feature.coef = feature.coef[-c(1), ] # remove Intercept
    feature.coef$auc = apply(data.frame( data.matrix[ind,feature.coef$name]) , 2,function(x){
      suppressMessages(roc <- pROC::roc(label[ind], x)  )
      roc$auc
    })

    feature.coef
  }

  # identify candidates
  candidate.occurence <- data.frame( unlist( table(res.raw$name) ), stringsAsFactors = FALSE )
  colnames(candidate.occurence) <- c("Name","Freq")

  candidate.auc <- aggregate(res.raw%>%select(coefficient, auc), by = list(res.raw$name), FUN = mean)
  names(candidate.auc) <- c("Name", "Coefficient", "AUC")

  candidate.raw.res <- full_join(candidate.occurence, candidate.auc, by="Name")
  candidate.raw.res$Freq <- round(candidate.raw.res$Freq/n, 2)

  res = list(Raw = res.raw, Summarized = candidate.raw.res)
  res


}


#' Recursive feature elimnation
#'
#'
#' @param feature.df Row is sample, column is feature
#' @param group
#' @param functions  Default: lrFuncs. lrFuncs, rfFuncs http://topepo.github.io/caret/available-models.html  There are a number of pre-defined sets of functions for several models, including: linear regression (in the object lmFuncs), random forests (rfFuncs), naive Bayes (nbFuncs), bagged trees (treebagFuncs) and functions that can be used with caret’s train function (caretFuncs).
#' @param seed Default 111
#' @param scale Deafult TRUE
#' @param sizes Default c(1:5), The sizes determines the number of most important features the rfe should iterate.
#' @param repeats For repeated k-fold cross-validation only: the number of complete sets of folds to compute
#' @param method The external resampling method: boot, repeatedcv, cv, LOOCV or LGOCV (for repeated training/test splits)
#' @param number Either the number of folds or number of resampling iterations
#' @param cores cores for parallel
#'
#' @return
#' @export
#'
#' @examples loonR::feature.selection.RFE(miR.df, group, functions="lrFuncs")
#'
#' Recursive feature elimnation (rfe) offers a rigorous way to determine the important variables before you even feed them into a ML algo.
feature.selection.RFE <- function(feature.df, group, functions = "lrFuncs",
        seed = 111, scale = TRUE, sizes = c(1:10),
        repeats = 5, number = 5, method = "cv", cores = 50){

  subsets <- sizes

  trainData = feature.df
  if(scale){
    trainData <- loonR::scaleDF(trainData, byColumn = T)
  }

  trainData <- data.frame(Class=group,
                          trainData,
                          check.names = F)

  trainData$Class <- factor(trainData$Class,
                            levels = unique(trainData$Class),
                            labels = 0:(length(unique(trainData$Class))-1)  )


  library(caret)

  library(doParallel)
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  # Perform
  set.seed(seed)
  # #构建rfe函数的控制参数
  ctrl <- rfeControl(functions = get(functions),
                     method = method,
                     number = number,
                     repeats = repeats,
                     verbose = FALSE)

  profile <- rfe(x=data.frame(trainData[, -c(1)], check.names = T),
                 y=trainData$Class,
                 sizes = subsets,
                 rfeControl = ctrl)

  stopCluster(cl)


  plot.obj <- plot(profile, type = c('g','o'))


  res <- list(result = profile)
  res$candidates <- profile$optVariables
  res$plot <- plot.obj

  res

}












