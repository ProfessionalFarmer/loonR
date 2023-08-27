
#' Get one round CV result
#'
#' @param df data.frame. Row is sample, col is gene
#' @param label TRUE/FALSE vector
#' @param k Folds
#' @param times Default: 1. the number of pieces was use to predict. E.g. times = 5, will run 5 times, each time use a different piece to predict.
#' @param type "response" or "link"
#'
#' @return
#'
#' @examples
getOneRoundCVRes <- function(df, label, k, seed = 1, times = 1, type = "response"){

  set.seed(666+seed)
  require(caret)
  flds <- createFolds(label, k = k, list = FALSE, returnTrain = FALSE)

  # Every fold will be used in one round.
  res <- foreach::foreach(i= 1:sample(k,times), .combine = rbind) %do% {
  # only predict one fold
  #res <- foreach::foreach(i= 1:1, .combine = rbind) %do% {

    piece.ind <- which(flds==i)

    lg.df <- data.frame(label= factor( label[-piece.ind], levels = c(FALSE, TRUE), labels = c(0,1) ),
                        df[-piece.ind, ], check.names = FALSE)

    #colnames(lg.df) <- gsub("\\.", "-", colnames(lg.df) )

    set.seed(666)

    # The type="response" option tells R to output probabilities of the form P(Y = 1|X), as opposed to other information such as the logit.
    suppressWarnings( glm.fit <- glm(label ~ ., data = lg.df, family = binomial(logit)) )

    pre.res <- predict(glm.fit, df[piece.ind, ], type = type)

    data.frame(Name  = names(pre.res),
               Logit = pre.res,
               Label = label[piece.ind],
               stringsAsFactors = FALSE)
  }

  # sometimes logit will be very verylarge e+15. We dont use it
  if(mean(abs(res$Logit)) > 1e+3){
    res$Logit = 0
  }
  res

}


#' k fold n times cross validation
#'
#' @param df data.frame. Row is sample, col is gene
#' @param label True label
#' @param k Folds
#' @param n Repeat times
#' @param scale TRUE
#' @param type "response" or "link"
#' @param seed Default 1
#'
#' @return
#' @export
#'
#' @examples
#' cv.res <- loonR::cross.validation(miR.df[, cf.candidates],
#'             group == "Cancer",
#'             k = 10, n = 100)
cross.validation <- function(df = NULL, label = NULL, k = 5, n = 100, scale=TRUE, type = "response", seed = 1){

  # df = up.mirs.exp
  # label = mir.sample.group
  # k = 5
  # n = 100
  # type <- match.arg(type)

  if(is.null(df) | is.null(label) ){
    stop("provide a datafram or lables")
  }

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}
  df = data.frame(df,check.names = F)


  library(doParallel)
  library(foreach)
  require(dplyr)

  #cl = parallel::makeCluster(40)
  #doParallel::registerDoParallel(cl)

  cv.res <- foreach::foreach(i= 1:n, .combine = rbind, .packages="foreach", .export=c("getOneRoundCVRes")) %do% {
  #res <- foreach::foreach(i= 1:n, .combine = rbind) %do% {

    loonR:::getOneRoundCVRes(df, label, k, seed=seed + i -1, type = type)

  }

  cv.res.mean <- cv.res %>% dplyr::group_by(Name, Label) %>% dplyr::summarize(Mean = mean(Logit))

  #stopCluster(cl)
  cv.res.mean <- cv.res.mean[match(row.names(df), cv.res.mean$Name), ]

  p = loonR::roc_with_ci(cv.res.mean$Label, cv.res.mean$Mean, ci=FALSE, title = "Cross validation")

  cv.res.return = list()
  cv.res.return$data = cv.res.mean
  cv.res.return$Label = cv.res.mean$Label
  cv.res.return$Mean = cv.res.mean$Mean
  cv.res.return$Name = cv.res.mean$Name
  cv.res.return$ROC = p
  cv.res.return$Raw = cv.res

  cv.res.return

}


#' Build a logistic model
#'
#' @param df Column is gene/miRNA, Row is sample
#' @param group Two levels. Second unique variable is defined as experiment group
#' @param seed
#' @param scale
#' @param direction backward  c("both", "backward", "forward")
#' @param rms If TRUE, use rms instead of glm to build the model. Useful for validation and calibration function in rms package.
#'
#' @return  A list. list(model=glm.fit,
#'      StepwiseModel=elimination,
#'      eliminationCandidates=stringr::str_remove_all(names(unclass(coef(elimination))[-c(1)]),'`')
#' )
#' @export
#'
#' @examples
#'
#' data("LIRI")
#' lg.res <- loonR::build.logistic.model(LIRI[,3:5],LIRI$status)
#' lg.res
#'
build.logistic.model <- function(df, group, seed = 666, scale=TRUE, direction = "backward", rms = FALSE){

  cat("Pls note: Second unique variable is defined as experiment group\n")

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  lg.df <- data.frame(
    label = factor(group == unique(group)[2],
      levels = c(FALSE, TRUE), labels = c(0, 1)
    ),
    df, check.names = FALSE
  )

  set.seed(seed)
  if(!rms){

    # The type="response" option tells R to output probabilities of the form P(Y = 1|X), as opposed to other information such as the logit.
    suppressWarnings( glm.fit <- glm(label ~ ., data = lg.df, family = binomial(logit)) )
    elimination = step(glm.fit, direction = direction, trace = 0)
    res <- list(model=glm.fit,
         StepwiseModel=elimination,
         eliminationCandidates=stringr::str_remove_all(names(unclass(coef(elimination))[-c(1)]),'`')
    )


    if(length(res$eliminationCandidates)!=0){
      elimination.df <- lg.df[ , c("label", res$eliminationCandidates) ]
      elimination.df$risk.score = predict(res$StepwiseModel, elimination.df, type = "link")
      res$elimination.result = elimination.df
      res$elimination.ROC = loonR::roc_with_ci(elimination.df$label, elimination.df$risk.score, ci = FALSE)
      res$elimination.AUC = loonR::get.AUC(elimination.df$risk.score, elimination.df$label, raw = F)
      res$elimination.Waterfull = loonR::plot_waterfall(elimination.df$risk.score, elimination.df$label , yticks.labl = NA )
    }


  }else{

    library(rms)
    dd <- rms::datadist(lg.df)
    options(datadist="dd")

    glm.fit <- rms::lrm(label ~ .,lg.df, x = TRUE, y = TRUE)
    res <- list(model=glm.fit)
  }

  lg.df$risk.score <- predict(glm.fit, lg.df)
  res$YoudenIndex <- loonR::get.YoudenIndex(lg.df$risk.score, lg.df$label)

  p <- loonR::roc_with_ci(lg.df$label, lg.df$risk.score, ci = FALSE)

  res$data = lg.df
  res$ROC = p
  res$AUC = loonR::get.AUC(lg.df$risk.score, lg.df$label, raw = F)
  res$Waterfull = loonR::plot_waterfall(lg.df$risk.score, lg.df$label, yticks.labl = NA  )



  res

}

#' Build cox regression model
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param seed Default 666
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
#' # psm (parametric survival model) uses a survival model based on functions and their parameters.
#' # cph (Cox Proportional Hazards Model and Extensions) is using the cox model (and the Anderson-Gill model) which is based on the hazard functions.
#'
#' data(LIRI)
#' res = loonR::build.coxregression.model(LIRI[,3:5],LIRI$status, LIRI$time)
build.coxregression.model <- function(d.frame, status, time, seed=666, scale = TRUE){

  library("survival")
  library("survminer")

  if(scale){
    d.frame = scale(d.frame, center = TRUE, scale = TRUE)
    d.frame = data.frame(d.frame, stringsAsFactors = FALSE, check.names = F)
  }


  covariates <- colnames(d.frame)

  df <- data.frame(d.frame,
                   Time=time,
                   Status =status,
                   check.names = F )

  formula <- as.formula( paste0("Surv(time, status) ~ `",
                                 paste0(covariates, sep='', collapse = '` + `'),
                                 '`',
                                 sep='', collapse = ""
                                 )
               )


  set.seed(seed)

  cox.fit <- coxph( formula, data = df )
  res = list(model=cox.fit, data = df)
  res

}


#' Build a random forests model
#'
#' @param df Row is sample
#' @param group
#' @param seed
#' @param scale Default TRUE
#'
#' @return
#' @export
#'
#' @examples
build.randomforest.model <- function(df, group, seed = 666, scale=TRUE){

  cat("Pls note: Second unique variable is defined as experiment group\n")

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  lg.df <- data.frame(
    label = factor(group == unique(group)[2],
                   levels = c(FALSE, TRUE), labels = c(0, 1)
    ),
    df, check.names = FALSE
  )

  set.seed(seed)

  rf.fit <- randomForest::randomForest(lg.df[,-c(1)], lg.df$label, importance=TRUE, proximity=TRUE)

  colnames(rf.fit$votes) = paste0("class",colnames(rf.fit$votes))
  predicted.prob <- as.vector(rf.fit$votes[,2])
  roc = loonR::roc_with_ci(lg.df$label, predicted.prob, ci=FALSE)
  auc = loonR::get.AUC(predicted.prob, lg.df$label)


  list(model=rf.fit,
       predicted.label=rf.fit$predicted,
       predicted.prob=predicted.prob,
       importance=rf.fit$importance,
       confusion=rf.fit$confusion,
       votes=data.frame(rf.fit$votes, check.names = FALSE),
       ROC = roc,
       AUC = auc
  )


}


#' Build elastic net regression
#'
#' @param df
#' @param group
#' @param seed Default 666
#' @param scale Default TRUE
#' @param nfolds Default 10. number of folds
#' @param summaryFunction Default twoClassSummary for computes sensitivity, specificity and the area under the ROC curve. Please refer https://www.rdocumentation.org/packages/caret/versions/6.0-90/topics/defaultSummary
#'
#' @return
#' @export
#'
#' @examples
#' data("LIRI")
#' reg.res <- loonR::build.elastic.net.regression(LIRI[,3:5],LIRI$status)
#'
build.elastic.net.regression <- function(df, group, seed = 666, scale=TRUE, nfolds = 10, summaryFunction = 'twoClassSummary'){
  # https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r

  cat("Pls note: Second unique variable is defined as experiment group\n")

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  lg.df <- data.frame(
    label = factor(group == unique(group)[2],
                   levels = c(FALSE,TRUE), labels = c("Control","Exp")
    ),
    df, check.names = FALSE
  )

  library(glmnet)
  library(caret)
  set.seed(seed)

  # Set training control
  train_cont <- trainControl(method = "repeatedcv",
                             number = nfolds,
                             repeats = 5,
                             search = "random",
                             classProbs = TRUE,
                             summaryFunction = get(summaryFunction),
                             verboseIter = FALSE)


  # Train the model
  elastic_reg <- train(label ~ .,
                       data = lg.df,
                       method = "glmnet", family = "binomial",
                       tuneLength = 10,
                       trControl = train_cont)

  prob = predict(elastic_reg, df, type = "prob")$Exp
  data = data.frame(
    Prob = prob,
    Group = group
  )
  youden.index = loonR::get.YoudenIndex(data$Prob, data$Group)


  ################ new model by glm
  modelByglmnet <- loonR::build.lassoOrRidge.regression(
    df, group, scale = F, alpha = elastic_reg$bestTune$alpha, lambda = elastic_reg$bestTune$lambda
  )


  res = list(
    model = elastic_reg,
    modelByglmnet = modelByglmnet,
    BestTuningParameter = elastic_reg$bestTune,
    data = data,
    youden.index = youden.index)

  res

}


#' Build lasso or ridge regression
#'
#' @param df
#' @param group
#' @param seed Default 666
#' @param scale Default TRUE
#' @param alpha 1 for lasso, 0 for ridge
#' @param nfolds Default 10, number of folds
#' @param type.measure Default auc. Refer https://www.rdocumentation.org/packages/glmnet/versions/4.1-3/topics/cv.glmnet
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
#' data("LIRI")
#' reg.res <- loonR::build.lassoOrRidge.regression(LIRI[,3:5],LIRI$status)
#'
build.lassoOrRidge.regression <- function(df, group, seed = 666, scale=TRUE, alpha=1, lambda = NA, nfolds = 10, type.measure = "auc"){

  # https://bookdown.org/tpinto_home/Regularisation/ridge-regression.html#ride.prac
  cat("Pls note: Second unique variable is defined as experiment group\n")

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  df = as.matrix(df)

  label = as.numeric(group == unique(group)[2])

  if(alpha < 0 | alpha > 1){
    stop("Alpha should be 0 or 1")
  }else if(alpha==1){
    print("Perform lasso")
  }else if(alpha==0){
    print("Perform ridge")
  }else{
    print("Perform elastic")
  }

  library(glmnet)
  set.seed(seed)

  # find lambda
  if(is.na(lambda)){
  # First we need to find the amount of penalty, λ by cross-validation. We will search for the
  # λ that give the minimum MSE
  #Penalty type (alpha=1 is lasso and alpha=0 is the ridge)
  cv.lambda.lassoOrRidge <- cv.glmnet(x=df, y=label, family="binomial",
                               alpha = alpha, type.measure = type.measure, nfolds = nfolds)
  #MSE for several lambdas
  plot(cv.lambda.lassoOrRidge)

  #We can also see the impact of different λ s in the estimated coefficients. When λ is very high, all the coefficients are shrunk exactly to zero.
  #Lasso path
  plot(cv.lambda.lassoOrRidge$glmnet.fit,
       "lambda", label=FALSE)

  l.lasso.min <- cv.lambda.lassoOrRidge$lambda.min

  }else{
    l.lasso.min = lambda
  }

  lasso.model <- glmnet(x=df, y=label,
                        alpha  = alpha, family="binomial",
                        lambda = l.lasso.min)

  prob = predict(lasso.model, df, type = "response")[,1]
  data = data.frame(
    Prob = prob,
    Group = group
  )
  youden.index = loonR::get.YoudenIndex(data$Prob, data$Group)

  # elimination
  candidates = as.matrix(lasso.model$beta)
  candidates = rownames(candidates)[candidates!=0]

  # ROC
  p.roc = loonR::roc_with_ci(data$Prob, label = data$Group, ci = F)
  auc = loonR::get.AUC(data$Prob, label = data$Group)

  res=list(model = lasso.model,
           alpha = alpha,
           lambda = l.lasso.min,
           data = data,
           `beta(coef)`= data.frame( as.matrix( coef(lasso.model) ) ),
           candidates = candidates,
           youden.index = youden.index,
           ROC = p.roc,
           AUC = auc)
  res

}


#' Build Support Vector Machines
#'
#' @param df row is sample
#' @param group
#' @param seed
#' @param scale Default TRUE
#'
#' @return
#' @export
#'
#' @examples
build.SVM <- function(df, group, seed = 666, scale=TRUE){

  cat("Pls note: Second unique variable is defined as experiment group\n")

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  svm.df <- data.frame(
    label = factor(group == unique(group)[2],
                   levels = c(FALSE, TRUE), labels = c(0, 1)
    ),
    df, check.names = FALSE
  )

  if(!require(e1071)){
    BiocManager::install("e1071")
    require(e1071)
  }

  set.seed(seed)
  svm.fit <- svm(label ~ ., data = svm.df, probability = TRUE)

  pred <- predict(svm.fit, svm.df)
  tab <- table(Predicted = pred, Actual = svm.df$label)
  error.rate = 1-sum(diag(tab))/sum(tab)


  predicted.prob <- predict(svm.fit, svm.df, probability=TRUE )
  predicted.prob = attr(predicted.prob,"probabilities")[,2]
  roc = loonR::roc_with_ci(svm.df$label, predicted.prob, ci=FALSE)
  auc = loonR::get.AUC(predicted.prob, svm.df$label)


  list(model = svm.fit,
       predicted.label = pred,
       predicted.prob = predicted.prob,
       confusion = tab,
       error.rate = error.rate,
       ROC = roc,
       AUC = auc
  )


}



#' Get confusion matrix
#'
#' @param groups TRUE/FALSE label
#' @param rs Predicted score
#' @param cancer Cancer Name
#' @param best.cutoff If not set, use youden index instead
#'
#' @return A data.frame object
#' @export
#'
#' @examples
#' label = c(1,1,1,1,1,2,2,2,2,2) == 1
#' risk.probability = runif(10, min=0, max=100)
#' confusion_matrix(label, risk.probability, cancer = "ESCC")
confusion_matrix <- function(groups, risk.pro, cancer="Cancer", best.cutoff = NA){

  library(caret)
  if( anyNA(best.cutoff) ){

    best.cutoff <- coords(roc(groups, risk.pro), "best", transpose = TRUE, input="threshold", best.method="youden")
    if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
    }else{
      best.cutoff <- best.cutoff[c("threshold")]
    }

  }

  # remove NA. some risk score may have NA value
  ind <- which(!is.na(risk.pro))

  results <- caret::confusionMatrix( factor(risk.pro[ind] > best.cutoff), factor(groups[ind]) , positive = "TRUE" )
  results.tl <- as.table(results)


  result.reformat <- data.frame(matrix(nrow = 6, ncol = 3))
  rownames(result.reformat) <- c("Normal",cancer,"Totals","Correct","Sensitivity","Specificity")
  colnames(result.reformat)  <- c("Normal",cancer,NA)

  result.reformat[1:2,1:2] <- results.tl
  result.reformat[2,3] <- c("Totals")

  result.reformat[c("Totals"),1:2] <- colSums(results.tl)
  result.reformat[c("Totals"),3]   <- sum(results.tl)

  result.reformat[c("Correct"), 1:2] <- diag(results.tl)
  result.reformat[c("Correct"), 3] <-  sum( diag(results.tl) )


  value = round(as.data.frame( as.matrix(results, what = "classes"))[c("Sensitivity"),][1], 3)
  result.reformat[c("Sensitivity"), 2] <- as.character(value)

  value = round(as.data.frame( as.matrix(results, what = "classes"))[c("Specificity"),][1], 3)
  result.reformat[c("Specificity"), 1] <- as.character(value)
  result.reformat[c("Specificity"), 3] <- as.character( round(sum( diag(results.tl) )/sum(results.tl),  3)  )

  result.reformat
  # as.matrix(results, what = "overall")
  # as.matrix(results, what = "classes")

}




#' Confusion matrix visulization
#'
#' @param groups True label
#' @param rs Predicted score
#' @param best.cutoff If not set, use youden index instead
#'
#' @return
#' @export
#'
#' @examples
#' draw_confusion_matrix(label,risk.probability)
draw_confusion_matrix <- function(groups, risk.pro, best.cutoff = NA) {

  cat("Pls note: Second unique variable is defined as experiment group\n")

  library(caret)
  library(pROC)
  if( anyNA(best.cutoff) ){

    best.cutoff <- coords(roc(groups, risk.pro), "best", transpose = TRUE, input="threshold", best.method="youden")
    if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
    }else{
      best.cutoff <- best.cutoff[c("threshold")]
    }

  }

  # remove NA. some risk score may have NA value
  ind <- which(!is.na(risk.pro))


  # https://stackoverflow.com/questions/23891140/r-how-to-visualize-confusion-matrix-using-the-caret-package/42940553
  cm <- caret::confusionMatrix( factor(risk.pro[ind] > best.cutoff), factor(groups[ind]==unique(groups)[2]) , positive = "TRUE" )


  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)

  # create the matrix
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, unique(groups)[1], cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, unique(groups)[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, unique(groups)[1], cex=1.2, srt=90)
  text(140, 335, unique(groups)[2], cex=1.2, srt=90)

  # add in the cm results
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')

  # add in the specifics
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)

  # add in the accuracy information
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)

}




#' Leave one out cross validation
#'
#' @param df Row is sample, column is gene/minRA/variat
#' @param label
#' @param seed Default 999
#' @param scale TRUE
#'
#' @return
#' @export
#'
#' @examples
loo.cv <- function(df, label, seed=999, scale=TRUE){

  if(scale){df = scale(df, center = TRUE, scale = TRUE)}

  lg.df = data.frame(Label = label, df,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

  n.samples = nrow(lg.df)
  library(foreach)
  library(dplyr)
  loo.res <- foreach::foreach(i=1:n.samples, .combine = rbind) %do%{
    loo.sample.data = lg.df[i,]

    set.seed(seed)
    train.sample.data = data.frame( lg.df[-i,],
                                    stringsAsFactors = FALSE,
                                    check.names = FALSE)

    suppressWarnings( glm.fit <- glm(Label ~ ., data = train.sample.data, family = binomial(logit)) )

    predicted.score = predict(glm.fit, loo.sample.data)
    sample.label = lg.df[i,"Label"]
    res <- c(predicted.score, as.character(sample.label))
    names(res) = c("Score", "Label")
    res
  }

  row.names(loo.res) <- row.names(lg.df)
  loo.res = data.frame( loo.res, check.names = FALSE, stringsAsFactors = FALSE)
  loo.res$Score = as.numeric(loo.res$Score)
  loo.res
}

#' Leave one out cross validation for cox regression
#'
#' @param df Row is sample, column is gene/minRA/variat
#' @param seed Default 999
#' @param status
#' @param time
#' @param label Default NA. Recommend to provide such information
#' @param scale TRUE
#'
#' @return
#' @export
#'
#' @examples
loo.cv.cox <- function(df, status, time,  seed=999, label=NA, scale =TRUE){

  if(scale){
    df = scale(df, center = TRUE, scale = TRUE)
    df = data.frame(df, check.names = T, stringsAsFactors = F)
  }

  library("survival")
  library("survminer")

  covariates <- colnames(df)

  n.samples = nrow(df)

  library(foreach)
  library(dplyr)
  loo.res <- foreach::foreach(i=1:n.samples, .combine = rbind) %do%{
    loo.sample.data = df[i, ]

    set.seed(seed)
    cox.fit = loonR::build.coxregression.model(df[-c(i),], status[-c(i)], time[-c(i)] )

    predicted.score = predict(cox.fit, loo.sample.data)

    sample.label = label[i]
    res <- c(predicted.score, as.character(sample.label))
    names(res) = c("Score", "Label")
    res
  }

  row.names(loo.res) <- row.names(df)
  loo.res = data.frame( loo.res, check.names = FALSE, stringsAsFactors = FALSE)
  loo.res$Score = as.numeric(loo.res$Score)
  loo.res
}



#' Get model performace, sensitivity, sepcificity and others
#'
#' @param pred Predicted score
#' @param labels T/F label
#' @param best.cutoff
#' @param digit 2 for 0.01, 3 for 0.001
#' @param boot.n Default 2000
#' @param specify.sen User can input one or more specific sensitivity
#' @param specify.spe User can input one or more specific specificity
#' @param OR Default F. If calculate OR or not
#'
#' @return
#' @export
#'
#' @examples
#' risk.probability = c(1:10)
#' labels = rep(c(1,0,1),c(2,6,2))
#' loonR::get_performance(risk.probability, labels)
get_performance <- function(pred, labels, best.cutoff =NA, digit = 2, boot.n = 2000, specify.sen = NULL, specify.spe = NULL, OR = F){  #  x="best", input = "threshold"

  labels = factor(labels)

  if(length( setdiff(unique(labels), c(FALSE,TRUE)) ) != 0 ){
    warning("Pls note second label will be used as TRUE label when get performance")
    labels = labels == unique(labels)[2]
  }


  # 这个更简单
  # reportROC::reportROC(gold = groups=="T", predictor = tmp.mir.exp, plot = F)

  library(pROC)
  library(caret)

  #pred <- df_plot[tmp.ind,]$RS
  #lables <- as.factor(tmp.label)
  # pred=kumamoto_rs_cv
  # labels= kumamoto$label
  # x="best"
  # input="threshold"

  input="threshold"

  string.format <- paste0("%.",digit,"f (%.",digit,"f-%.",digit,"f)")


  if( anyNA(best.cutoff) ){

    best.cutoff <- coords(roc(labels, pred), "best", transpose = TRUE, input="threshold", best.method="youden")
    if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
    }else{
      best.cutoff <- best.cutoff[c("threshold")]
    }

  }

  # remove NA. some risk score may have NA value
  ind <- which(!is.na(pred))

  # 20220427  confusion matrix
  confusion.matrix = loonR::confusion_matrix(labels,pred)
  true.cancer = as.integer(confusion.matrix[2,2])
  true.healthy = as.integer(confusion.matrix[1,1])

  # calculate OR
  results <- caret::confusionMatrix( factor(pred[ind] > best.cutoff), factor(labels[ind]) , positive = "TRUE" )
  results.tl <- as.table(results)

  or <- vcd::oddsratio(results.tl,log=FALSE)
  or.ci <- confint(vcd::oddsratio(results.tl,log=FALSE), level = 0.95)
  or <- sprintf(string.format, exp(or$coefficients),or.ci[1],or.ci[2]) # Odds ratio

  # count sample
  t.c <- sum(labels==TRUE)
  n.c <- length(labels) - t.c

  # get AUC CI
  set.seed(100)
  roc.obj <- pROC::roc(labels, pred, ci=TRUE, plot=FALSE)
  auc <- pROC::ci(roc.obj, boot.n = boot.n)[c(2, 1, 3)]
  auc <- round(auc, digit)
  auc <- sprintf(string.format, auc[1], auc[2], auc[3])

  # get performance
  set.seed(100)
  rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
            "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "precision", "recall")
  others <- pROC::ci.coords(roc.obj, x = best.cutoff, boot.n = boot.n, input = input, ret = rets, best.policy = "random", transpose = TRUE)
  # to be continue

  # ppv Precision  https://www.rdocumentation.org/packages/pROC/versions/1.16.2/topics/coords
  # sensitivity recall
  #f1.score <- (2 * others$precision[1,c("50%")] * others$recall[1,c("50%")]) / (others$precision[1,c("50%")] + others$recall[1,c("50%")])
  #f1.score <- round(f1.score,2)

  res <- c(
    paste0(n.c, " [", true.healthy,"]"),
    paste0(t.c, " [", true.cancer,"]"),
    auc[1], # AUC
    sprintf(string.format, others$accuracy[1,c("50%")],  others$accuracy[1,c("2.5%")],  others$accuracy[1,c("97.5%")]), # accuracy
    sprintf(string.format, others$precision[1,c("50%")], others$precision[1,c("2.5%")], others$precision[1,c("97.5%")]),# precision
    sprintf(string.format, others$recall[1,c("50%")],    others$recall[1,c("2.5%")],    others$recall[1,c("97.5%")]),#recall
    sprintf(string.format, others$specificity[1,c("50%")],others$specificity[1,c("2.5%")],others$specificity[1,c("97.5%")]),#specificity
    sprintf(string.format, others$npv[1,c("50%")],others$npv[1,c("2.5%")],others$npv[1,c("97.5%")]),#npv
    or  # Odds ratio
  )

  #names(res) <- c("Ncount", "Tcount", "AUC (CI)",  "Accuracy", "Precision", "Recall", "Specificity", "NPV","Odds Ratio")
  #recall = sensititivity, precision=PPV
  names(res) <- c("Ncount [TN]", "Tcount [TP]", "AUC",  "Accuracy",  "PPV", "Sensitivity", "Specificity", "NPV","Odds Ratio")

  # decide whether to show OR
  if(!OR){
    res = res[-c(length(res))] # the last one is OR ratio
  }



  # 20220310 update. Get specific sensitivity and specificity
  if(! is.null(specify.sen) ){
    specify.sen = as.numeric(specify.sen)

    library(foreach)
    updated.20220310 = foreach(x=specify.sen, .combine = c) %do%{
      single.point.res = ci.coords(roc.obj, x, input = "sensitivity", ret = "specificity", boot.n = boot.n)
      sprintf(string.format, single.point.res$specificity[1,c("50%")],  single.point.res$specificity[1,c("2.5%")],  single.point.res$specificity[1,c("97.5%")])
    }
    names(updated.20220310) = paste("Sp @ ", specify.sen, " Sn", sep = "")

    res = c(res, updated.20220310)
  }
  if(! is.null(specify.spe) ){
    specify.spe = as.numeric(specify.spe)

    library(foreach)
    updated.20220310 = foreach(x=specify.spe, .combine = c) %do%{
      single.point.res = ci.coords(roc.obj, x, input = "specificity", ret = "sensitivity", boot.n = boot.n)
      sprintf(string.format, single.point.res$sensitivity[1,c("50%")],  single.point.res$sensitivity[1,c("2.5%")],  single.point.res$sensitivity[1,c("97.5%")])
    }
    names(updated.20220310) = paste("Sn @ ", specify.spe, " Sp", sep = "")

    res = c(res, updated.20220310)
  }


  # value (conf) -> value  conf
  rname <- names(res)
  value <- sapply(strsplit(res," \\("), function(x) paste(x[1]) )


  tmp <- sapply(strsplit(res," \\("), function(x) paste(x[2]) )
  conf <- sapply(strsplit(tmp,"\\)"), function(x) paste(x[1]) )
  # N and T don't have CI
  conf[1:2] <- ''


  df <- data.frame(Name = rname, value = value, confidence = conf)


  df$Formated <- paste0(df$value, " ", "(", df$confidence, ")"  )
  df$Formated <- stringr::str_remove_all(df$Formated, " \\(\\)")



  df

}


#' Get individual_candidates_performance
#'
#' @param scores Df or list, if data.frame, column is feature.
#' @param labels Vector or list
#' @param best.cutoff
#' @param digit 2 for 0.01, 3 for 0.001
#' @param specify.sen User can input one or more specific sensitivity
#' @param specify.spe User can input one or more specific specificity
#' @param OR Default F. If calculate OR or not
#'
#' @return
#' @export
#'
#' @examples
get_individual_candidates_performance <- function(scores, labels, best.cutoff =NA, digit = 2, specify.sen = NULL, specify.spe = NULL, OR=F){

  oldw <- getOption("warn")
  options(warn = -1)

  set.seed(100)

  require(pROC)

  if(inherits(scores, "list") ) {
    performance.list <- lapply(1:length(scores), function(i){
      index <- !is.na(scores[[i]])
      loonR::get_performance(scores[[i]][index], labels[[i]][index], best.cutoff = best.cutoff, digit = digit, specify.sen = specify.sen, specify.spe = specify.spe)

    })
    names(performance.list) <- names(scores)

  }else{
    performance.list <- apply(scores, 2, function(x) loonR::get_performance(x, labels, best.cutoff = best.cutoff, digit = digit, specify.sen = specify.sen, specify.spe = specify.spe, OR = OR) )
    # https://stackoverflow.com/questions/57608056/how-to-change-legend-description-on-ggroc-function
  }

  performance.formated.list <- lapply(performance.list, function(x){
    dplyr::pull(x,4)
  })

  performance.formated.df <- do.call(cbind,performance.formated.list)
  row.names(performance.formated.df) <- row.names(performance.list[[1]])

  data.frame(performance.formated.df, check.names = F)

}




#' Return a table including all metrics across different threshold
#'
#' @param score
#' @param label
#'
#' @return
#' @export
#'
#' @examples
#' score = c(10:20)
#' label = c(0,0,1,0,1,0,1,1,0,1,1)
#' loonR::get_all_performance_table(score, label)
get_all_performance_table <- function(score = NULL, label = NULL){

  if(is.null(score) | is.null(label) ){
    stop("Pls input score and label")
  }
  library(pROC)
  roc.obj <- roc(response = label, predictor = score)
  coordinates <- coords(roc.obj, x = "all", input = "threshold",
                        ret = c("threshold", "sensitivity", "specificity", "accuracy", "npv", "ppv"))
  coordinates

}




#' Variate logistic analysis
#' Row: sample, Column: gene expression
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param label True label
#' @param scale
#'
#' @return
#' @export
#'
#' @examples multivariate_or(data.frame, label)
multivariate_or <- function(d.frame, label, scale=TRUE){

  if(scale){df = scale(d.frame, center = TRUE, scale = TRUE)}

  res <- glm(Event ~ . , data = data.frame(d.frame, Event=label, check.names = FALSE), family=binomial(logit))
  #exp( coef(res) )
  #summary(res)

  res <- data.frame(
    exp( cbind(coef(res), confint(res) )  ) ,
    Estimate = as.vector(summary(res)$coefficients[,1] ),
    Pr = as.vector(summary(res)$coefficients[,4] )
  )

  row.names(res) <- stringr::str_remove_all(row.names(res), "`")

  res <- data.frame(Vairate=row.names(res), round(res,3) )
  res[,2:ncol(res)] <- data.frame(lapply(res[,2:ncol(res)],as.numeric))

  colnames(res) <- c("Variate","OR", "2.5%", "97.5%", "Estimate" , "Pr")
  res

}



#' Variate logistic analysis
#' Row: sample, Column: gene expression
#' score E.g.: Gene or miRNA expression, or risk score
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param label True Sample label
#' @param scale
#'
#' @return c(OR, 2.5% CI, 97.5% CI)
#' @export
#'
#' @examples
#' data("LIRI")
#' loonR::univariate_or(LIRI[,3:6], LIRI$status)
univariate_or <- function(d.frame, label, scale=TRUE){

  library(foreach)

  if(scale){d.frame = scale(d.frame, center = TRUE, scale = TRUE)}

  all.res <- foreach(i=1:ncol(d.frame), .combine = rbind) %do%{
    #   for(i in 1:ncol(d.frame) ){
    res <- glm(Event ~ . , data = data.frame(Score = d.frame[,i], Event=label), family=binomial(logit))
    #exp( coef(res) )
    #summary(res)

    # column name
    c.name <- colnames(d.frame)[i]

    # 用t= as.numeric(questionr::odds.ratio(res)[2,])比较简单
    # t= as.numeric(questionr::odds.ratio(res)[2,])

    res <- c( exp( cbind(coef(res), confint(res) )  )[2,] , #  c("OR", "2.5 %", "97.5 %")
              summary(res)$coefficients[2,1], # "Estimate"
              summary(res)$coefficients[2,4]  # "Pr"
    )
    # names(res) <- c("OR", "2.5 %", "97.5 %")
    c( Name=c.name, round(res,3))

  }

  all.res <- data.frame(all.res, stringsAsFactors = F)
  row.names(all.res) <- all.res$Name
  colnames(all.res) <- c("Variate", "OR", "2.5%", "97.5%","Estimate" , "Pr")

  all.res[,2:ncol(all.res)] <- data.frame(lapply(all.res[,2:ncol(all.res)],as.numeric))
  all.res

}



#' Univariate Cox analysis
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param scale
#' @param max.time NA
#'
#' @return
#' @export
#'
#' @examples
univariate_cox <- function(d.frame, status, time, scale=TRUE, max.time = NULL){

  library(foreach)
  library("survival")
  library("survminer")

  if(scale){d.frame = scale(d.frame, center = TRUE, scale = TRUE)}


  if(!is.null(max.time)){
    if(!is.numeric(max.time)){stop("Pls input days or months for max.time")}
    # There not event in the maximum point.
    status[time >= max.time & status] = 0
    time[time >= max.time] = max.time
  }

  all.res <- foreach(i=1:ncol(d.frame), .combine = rbind) %do%{


    res <- coxph( Surv(time, status) ~ Score ,
                  data = data.frame(Score = d.frame[,i],
                                    Time=time, Status =status, check.names = F ) )
    #exp( coef(res) )
    #summary(res)

    # column name
    c.name <- colnames(d.frame)[i]
    i.coef <- coef(res)
    i.hazzardRation <- exp( coef(res) )
    i.hr.upper = summary(res)$conf.int[4]
    i.hr.lower = summary(res)$conf.int[3]
    i.pvalue = summary(res)$waldtest[3]
    i.SE <- sqrt( diag( res$var ))

    res <- c( i.coef, # coef
              i.hazzardRation, # HR
              i.hr.lower, # lower
              i.hr.upper, # upper
              i.pvalue, i.SE
    )
    c( Name=c.name, formatC(res, format = "e", digits = 2) )
    #c( Name=c.name, res )

  }

  all.res <- data.frame(all.res, stringsAsFactors = F)
  row.names(all.res) <- all.res$Name
  colnames(all.res) <- c("Variate", "coefficient" , "HR", "lower 95%", "upper 95%", "Pr", "SE")

  all.res[,2:ncol(all.res)] <- data.frame(lapply(all.res[,2:ncol(all.res)],as.numeric))
  all.res

}

#' Univariate Cox analysis
#'
#' Copy http://www.sthda.com/english/wiki/cox-proportional-hazards-model
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param scale
#' @param max.time NA
#'
#' @return
#' @export
#'
#' @examples
univariate_cox_sthda <- function(d.frame, status, time, scale=TRUE, max.time = NULL){

  if(scale){d.frame = scale(d.frame, center = TRUE, scale = TRUE)}

  library(survival)
  # Clone frome http://www.sthda.com/english/wiki/cox-proportional-hazards-model
  covariates <- colnames(d.frame)


  if(!is.null(max.time)){
    if(!is.numeric(max.time)){stop("Pls input days or months for max.time")}
    # There not event in the maximum point.
    status[time >= max.time & status] = 0
    time[time >= max.time] = max.time
  }

  df <- data.frame(d.frame,
                   Time=time,
                   Status =status,
                   check.names = F )

  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(Time, Status)~`', x,'`', sep="")))

  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df)})
  # Extract data
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  as.data.frame(res)

}

#' Title
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
multivariate_cox <- function(d.frame, status, time, scale=TRUE){

  library("survival")
  library("survminer")

  if(scale){
    d.frame = scale(d.frame, center = TRUE, scale = TRUE)
    d.frame = data.frame(d.frame, stringsAsFactors = FALSE, check.names = F)
  }

  covariates <- colnames(d.frame)

  df <- data.frame(d.frame,
                   Time=time,
                   Status =status,
                   check.names = F )

  formula <- as.formula( paste0("Surv(time, status) ~ `",
                                paste0(covariates, sep='', collapse = '` + `'),
                                '`',
                                sep='', collapse = ""
  )
  )


  res <- coxph( formula , data = df )

  res <- data.frame(
    exp( cbind(coef(res), confint(res) )  ) ,
    Estimate = as.vector(summary(res)$coefficients[,1] ),
    Pr = as.vector(summary(res)$coefficients[,5] )
  )

  res <- data.frame(Vairate=row.names(res), round(res,3) )
  res[,2:ncol(res)] <- data.frame(lapply(res[,2:ncol(res)],as.numeric))

  colnames(res) <- c("Variate","HR", "2.5%", "97.5%", "Estimate" , "Pr")
  res

}


#' Plot miR's correlation
#'
#' @param df Row: miR expression, Column: Sample
#' @param title miRs' correlation
#' @param cl.lim c(-0.5,1) correlation lims
#'
#' @return Plot
#' @export
#'
#' @examples plot_miRCorrelation(data[candi,])
#'
plot_miRCorrelation <- function(df, title="miRs' correlation", cl.lim = c(-0.5,1)){
  cor.res <- cor(t(df))
  cor.res <- round(cor.res, 3)#保留两位小数

  library(corrplot)#先加载包
  corrplot(cor.res, type = "upper",
           order = "hclust", tl.col = "black", tl.srt = 90, mar=c(0,0,2,0),
           cl.lim = c(-0.5,1), addgrid.col=FALSE, title = title ) +
  cowplot::theme_cowplot(font_family = "Arial")

}






#' plot ROC with CI. Stratified bootstrap 2000 times
#'
#' @param label True label
#' @param rs Predicted score
#' @param font Font
#' @param palette
#' @param legend.pos
#' @param title
#' @param fontsize
#' @param panel panel的第一列为factor，event，class等
#' @param SE.SP Whether to show SE and SP instead of CI of AUC
#' @param ci Default draw CI interval
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data("LIRI")
#' loonR::roc_with_ci(LIRI$status, LIRI$ANLN,
#'    panel = data.frame(Evebt=LIRI$status, LIRI[,c(3,4)])  )
#'
#'
#' roc_with_ci(labels, rs,
#'             font = "Arial",
#'             palette = "jama_classic",
#'             title = "HGD vs Healthy",
#'             panel = data.frame(Evebt=labels, data)
#' )
roc_with_ci <- function(label, rs, font = "Arial", palette = "jama", legend.pos = c(0.4, 0.2), title = NULL, fontsize = 16, panel = NULL, SE.SP = FALSE, ci = TRUE) {

  library(pROC)
  library(ggplot2)
  obj = roc(label, rs, ci=TRUE, plot=FALSE)

  # panel的第一列为factor，event，class等

  library(doParallel)
  registerDoParallel(40)
  set.seed(100)


  if(ci){
    library(doParallel)
    registerDoParallel(40)
    set.seed(100)

    ciobj <- ci.se(obj, specificities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)
    dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                         se.lower = ciobj[, 1],
                         se.upper = ciobj[, 3])

    ciobj <- ci.sp(obj, sensitivities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)
    dat.ci$y <- as.numeric(rownames(ciobj))
    dat.ci$sp.lower <- ciobj[, 1]
    dat.ci$sp.upper <- ciobj[, 3]
  }else{

  }

  aucs <- pROC::ci(obj)[c(2, 1, 3)]
  others <- pROC::coords(obj, "best", ret = c("sensitivity", "specificity"), best.policy = "omit")

  if(!SE.SP){
    annot <- sprintf("AUC %.2f\n(%.2f-%.2f)", aucs[1], aucs[2], aucs[3])
  }else{
    annot <- sprintf("AUC %.2f\nSensitivity %.2f\nSpecificity %.2f", aucs[1],others[1,1],others[1,2])
  }

  if(ci){
    p <- ggroc(obj, colour = loonR::get.palette.color(palette, n=length(annot)) , size=0.93, legacy.axes = TRUE ) +
      labs(x = "1 - Specificity", y = "Sensitivity") +
      scale_color_manual(labels = annot) + annotate("text", x = 0.55, y = 0.1, label =annot, size = fontsize/3) +
      theme(legend.position = legend.pos, legend.title = element_blank()) +
      theme_classic() + cowplot::theme_cowplot(font_family = font) +
      #geom_abline( slope = 1,  intercept = 1, linetype = "dashed", alpha = 0.7) +
      geom_abline(linetype = "dashed", alpha = 0.3) +
      coord_equal() +
      geom_ribbon(
        data = dat.ci,
        aes(x = 1-x, xmin = 1-sp.upper, xmax = 1-sp.lower, y=y, ymin = se.lower, ymax = se.upper), # note 1 -
        fill = loonR::get.palette.color(palette, length(annot)),
        alpha = 0.1
      ) + ggtitle(title)  + theme(plot.title = element_text(hjust = 0.5)) # capture.output(obj$ci)
  }else{
    p <- pROC::ggroc(obj,
               colour = loonR::get.palette.color(palette, n=length(annot)) ,
               size=0.93,
               legacy.axes = TRUE ) +
      labs(x = "1 - Specificity", y = "Sensitivity") +
      scale_color_manual(labels = annot) + annotate("text", x = 0.55, y = 0.1, label =annot, size = fontsize/3) +
      theme(legend.position = legend.pos, legend.title = element_blank()) +
      theme_classic() + cowplot::theme_cowplot(font_family = font) +
      #geom_abline( slope = 1,  intercept = 1, linetype = "dashed", alpha = 0.7) +
      geom_abline(linetype = "dashed", alpha = 0.3) +
      coord_equal() +
      ggtitle(title)  + theme(plot.title = element_text(hjust = 0.5)) # capture.output(obj$ci)
  }

  ### if else
  if(is.null(panel)){
    p
  }else{

    panel = as.data.frame(panel)

    ## calculate panel member's sensitivity and specificity
    library(foreach)
    library(doParallel)
    registerDoParallel(20)

    p.mem.res <-
      foreach(p.mem = colnames(panel[,-c(1)]), .combine = rbind) %dopar% {
        cat(p.mem,"\n")
        set.seed(100)

        perf = loonR::get_performance(panel[,p.mem], panel[,c(1)] )

        tmp.res <- c(
          stringr::str_split( perf[c("Specificity"),c("confidence")],"-")[[1]][1],
          perf[c("Specificity"),c("value")],
          stringr::str_split( perf[c("Specificity"),c("confidence")],"-")[[1]][2],
          stringr::str_split( perf[c("Sensitivity"),c("confidence")],"-")[[1]][1],
          perf[c("Sensitivity"),c("value")],
          stringr::str_split( perf[c("Sensitivity"),c("confidence")],"-")[[1]][2]
        )

        as.numeric(tmp.res)

      }
    p.mem.res <- as.data.frame(p.mem.res)

    colnames(p.mem.res) <- c("specificity.low", "specificity", "specificity.high", "sensitivity.low", "sensitivity", "sensitivity.high")
    row.names(p.mem.res) <- colnames(panel[,-c(1)])


    # 1-specificity
    p + geom_segment(data=p.mem.res, color="#009090",
                     aes(x=1-specificity.high,
                         xend=1-specificity.low,
                         #ymin=(sensitivity)/2,
                         #ymax=(sensitivity)/2,
                         y=(sensitivity),
                         yend=(sensitivity))
    )  + # sensitivity
      geom_segment(data=p.mem.res,  color="#009090",
                   aes(x=1-specificity,
                       xend=1-specificity,
                       #ymin=sensitivity.low,
                       #ymax=sensitivity.high,
                       y=sensitivity.low,
                       yend=sensitivity.high)
      ) +
      geom_point(data=p.mem.res, mapping=aes(x=(1-specificity),
                                             y=(sensitivity)  #, ymin=0,ymax=0
                                             ),
                 size=3.5, shape=16, fill="#009090",color="#009090")

  }
  ## if else

}



#' Plot multiple ROCs in one figure
#'
#' @param scores A list or a data.frame. If list, labels shoule also be a list
#' @param labels A list or vector
#' @param font Arial
#' @param palette jama
#' @param legend.pos
#' @param title
#' @param panel
#' @param color specify the color manually
#' @param ci
#' @param SE.SP Whether to show SE and SP instead of CI of AUC
#'
#' @return
#' @export
#'
#' @examples multi_roc_with_ci(rss, labels, font = "Arial", palette = "jama")
multi_roc_with_ci <- function(scores, labels, font = "Arial", palette = "jama", legend.pos = c(0.4, 0.2), title = NULL, panel = NULL, color = NULL, ci =TRUE, SE.SP = FALSE) {
  oldw <- getOption("warn")
  options(warn = -1)


  # panel的第一列为factor，event，class等
  set.seed(100)

  require(pROC)
  require(ggplot2)

  if(inherits(scores, "list") ) {
    roclist <- lapply(1:length(scores), function(i){
      index <- !is.na(scores[[i]])
      pROC::roc(labels[[i]][index], scores[[i]][index])
    })
    names(roclist) <- names(scores)

  }else{
    roclist <- apply(scores, 2, function(x) pROC::roc(labels,x) )
    # https://stackoverflow.com/questions/57608056/how-to-change-legend-description-on-ggroc-function
  }


  dat.ci <- data.frame( x=NA, se.lower=NA, se.upper=NA, group = NA, y =NA, sp.lower = NA, sp.upper = NA )

  if (ci){
      for(group_name in names(roclist)){

        library(doParallel)
        registerDoParallel(40)
        set.seed(100)

        se.ciobj <- ci.se(roclist[[group_name]], specificities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)

        library(doParallel)
        registerDoParallel(40)
        set.seed(100)

        sp.ciobj <- ci.sp(roclist[[group_name]], sensitivities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)

        dat.ci.tmp <- data.frame(x = as.numeric(rownames(se.ciobj)),
                                 se.lower = se.ciobj[, 1],
                                 se.upper = se.ciobj[, 3],
                                 group = group_name,
                                 y = as.numeric(rownames(sp.ciobj)),
                                 sp.lower = sp.ciobj[, 1],
                                 sp.upper = se.ciobj[, 3]
        )

        dat.ci <- rbind(dat.ci, dat.ci.tmp)
        rm(dat.ci.tmp)
      }
      # remove first line
      dat.ci <- dat.ci[-c(1),]
  }else {
  }



  annot <- c()
  aucs <- c()
  for(group_name in names(roclist)){
    auc <- pROC::ci(roclist[[group_name]])[c(2, 1, 3)]
    others <- pROC::coords(roclist[[group_name]], "best", ret = c("sensitivity", "specificity"), best.policy = "omit")

    cat(sprintf("%s %.2f (%.2f-%.2f)\n", group_name, auc[1], auc[2], auc[3]) )
    #annot <- c(annot, sprintf("%s %.2f (%.2f-%.2f)", group_name, auc[1], auc[2], auc[3])  )
    #annot <- c(annot, sprintf("AUC %.2f\nSensitivity %.2f\nSpecificity %.2f", aucs[1],others[1,1],others[1,2]) )
    if(!SE.SP){
      annot <- c(annot, sprintf("%.2f (%.2f-%.2f)", auc[1], auc[2], auc[3])  ) # 常用，AUC CI
    }else{
      annot <- c(annot, sprintf("%.2f (SE %.2f, SP %.2f)", auc[1], others[1,1],others[1,2] )  ) # SE, SP
    }
    #



    aucs <- c(aucs, auc[1])
  }

  # Without CI information
  # annot <- paste0(stringr::str_pad(names(roclist), max(sapply(names(roclist), nchar))+1, "right"), "\t", sprintf("%.2f", aucs))
  # With Ci information
  annot <- paste0(stringr::str_pad(names(roclist), max(sapply(names(roclist), nchar))+1, "right"), "\t", annot)



  if(is.null(color)){
    colors <- loonR::get.palette.color(palette, length(annot))
  }else{
    colors <- color
  }

  names(colors) <- names(roclist)


  p <- pROC::ggroc(roclist, legacy.axes = TRUE, size=0.93) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    scale_color_manual(labels = annot, values = colors ) +
    theme_classic() + cowplot::theme_cowplot(font_family = font) +
    #geom_abline( slope = 1,  intercept = 1, linetype = "dashed", alpha = 0.7) +
    geom_abline(linetype = "dashed", alpha = 0.3) +
    coord_equal()

  if(ci){

    p <- p  +
      geom_ribbon(
        data = dat.ci, inherit.aes = FALSE, show.legend = FALSE,
        aes(x = 1-x, xmin = 1-sp.upper, xmax = 1-sp.lower, y =y, ymin = se.lower, ymax = se.upper, group=group, fill=as.factor(group)), # note 1 -
        alpha = 0.1
      ) + scale_fill_manual(values=colors)

  }

  # if so many ROCs, put legend in the right
  if(length(roclist)>5){
    p <- p + theme(legend.title = element_blank()) +
      ggtitle(title)  # capture.output(obj$ci)
  }else{
    p <- p + theme(legend.position = legend.pos, legend.title = element_blank()) +
      ggtitle(title)  # capture.output(obj$ci)
  }


  ### if else
  if(is.null(panel)){
    p
  }else{
    ## calculate panel member's sensitivity and specificity
    library(foreach)
    library(doParallel)
    registerDoParallel(20)

    p.mem.res <-
      foreach(p.mem = colnames(panel[,-c(1)]), .combine = rbind) %dopar% {
        cat(p.mem,"\n")
        set.seed(100)

        perf = loonR::get_performance(panel[,p.mem], panel[,c(1)] )

        tmp.res <- c(
          stringr::str_split( perf[c("Specificity"),c("confidence")],"-")[[1]][1],
          perf[c("Specificity"),c("value")],
          stringr::str_split( perf[c("Specificity"),c("confidence")],"-")[[1]][2],
          stringr::str_split( perf[c("Sensitivity"),c("confidence")],"-")[[1]][1],
          perf[c("Sensitivity"),c("value")],
          stringr::str_split( perf[c("Sensitivity"),c("confidence")],"-")[[1]][2]
        )

        as.numeric(tmp.res)
      }
    p.mem.res <- as.data.frame(p.mem.res)

    #cat(dim(p.mem.res))


    colnames(p.mem.res) <- c("specificity.low", "specificity", "specificity.high", "sensitivity.low", "sensitivity", "sensitivity.high")
    row.names(p.mem.res) <- colnames(panel[,-c(1)])
    options(warn = oldw)

    # 1-specificity
    p + geom_segment(data=p.mem.res, color="#009090",
                     aes(x=1-specificity.high,
                         xend=1-specificity.low,
                         #ymin=(sensitivity)/2,
                         #ymax=(sensitivity)/2,
                         y=(sensitivity),
                         yend=(sensitivity))
    )  + # sensitivity
      geom_segment(data=p.mem.res,  color="#009090",
                   aes(x=1-specificity,
                       xend=1-specificity,
                       #ymin=sensitivity.low,
                       #ymax=sensitivity.high,
                       y=sensitivity.low,
                       yend=sensitivity.high)
      ) +
      geom_point(data=p.mem.res, mapping=aes(x=(1-specificity),
                                             y=(sensitivity)  #, ymin=0,ymax=0
      ),
      size=3.5, shape=16, fill="#009090",color="#009090")


  }
  ## if else

}

library(pROC)
data(aSAH)
roc_obj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)


#' Get all point if you want to draw yourself
#'
#' @param pred
#' @param label
#' @param ci if CI
#' @param boot.n 2000
#'
#' @return
#' @export
#'
#' @examples
get.all.point.ROC <- function(pred, label, ci=FALSE, boot.n = 2000){

  res = list()

  library(pROC)
  roc_obj <- roc(label, pred, quiet=TRUE)

  all.points = coords(roc = roc_obj, transpose = T, ret = c("all") )
  all.points = data.frame(t(all.points), check.names = F)
  res$all.points = all.points

  if(ci){

    all.points.CI.raw = ci.coords(roc = roc_obj, x = all.points$threshold, input = "threshold", transpose = T, ret = c("all") , boot.n = boot.n)
    tmp = lapply(names(all.points.CI.raw), function(x){

        tmp.df = all.points.CI.raw[[x]]
        colnames(tmp.df) = paste(x, " (", colnames(tmp.df), ")", sep ="" )
        colnames(tmp.df) = stringr::str_remove_all(colnames(tmp.df), " \\(50%\\)")

        tmp.df

    })
    names(tmp) = names(all.points.CI.raw)
    all.points.CI.raw = tmp
    rm(tmp)

    res$all.points.CI.raw = all.points.CI.raw

    # merge all the tables
    res$all.points.CI.merged =  do.call("cbind",all.points.CI.raw)
  }
  res

}



#' Convert logti to probability
#'
#' @param logit
#'
#' @return
#' @export
#'
#' @examples
logit2prob <- function(logit){
  odds <- exp(logit)
  odds[ is.infinite(odds) ] = 500
  prob <- odds / (1 + odds)
  return(prob)
}

#' Convert Probability to Logit
#'
#' Defined simply as \code{log(x / (1 - x))}.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector.
#'
#' @export
prob2logit <- function(x) {
  out <- log(x / (1 - x))
  return(out)
}


#' Waterfall plot
#'
#' @param risk.score
#' @param label
#' @param xlab Risk probability
#' @param palette
#' @param title
#' @param yticks.labl seq(0,1,by = 0.25) c(0,0.25,0.5,0.75,1)
#' @param sample If you want to show sample label, provide a vector
#' @param rotate.x Default 90. Rotate x axis text.
#'
#' @return
#' @export
#'
#' @examples loonR::plot_waterfall(average.riskscore$Mean, average.riskscore$Label, xlab = "Risk probability")
plot_waterfall <- function(risk.score, label, xlab = "Risk probability", palette = "jco", title = "", yticks.labl = c(0,0.25,0.5,0.75,1), sample = NA, rotate.x = 90 ){

    if( ( max(risk.score)>1|min(risk.score)<(-1) ) & xlab == "Risk probability")  {
      xlab = "Risk score"
      if( !is.na(yticks.labl) & sum(yticks.labl == c(0,0.25,0.5,0.75,1)) == length(yticks.labl) ){
        yticks.labl=NA
      }
    }

    library(ggpubr)
    risk.score = risk.score
    idx = order(risk.score)
    risk.score <- risk.score[idx]
    label = label[idx]
    if(anyNA(sample)){
      tmp.df <- data.frame(Risk=risk.score, Class = label, ID=1:length(risk.score) )
    }else{
      sample = sample[idx]
      tmp.df <- data.frame(Risk=risk.score, Class = label, ID=sample )
    }

    colnames(tmp.df)[1] <- xlab



    p <- ggbarplot(tmp.df, x = "ID", y = xlab, xlab = "",
                   color = "Class", fill = "Class",
                   palette = palette, legend = "right", title = title)
    if(anyNA(sample)){
      p <- p + rremove("x.text") + rremove("x.axis") + rremove("x.ticks")
    }else{
      p <- p + rotate_x_text(rotate.x)
    }

    if(!anyNA(yticks.labl)){
      p <- p + scale_y_continuous(labels = yticks.labl)
    }

   p
}



#' plot predictive probability distribution
#'
#' @param group
#' @param pred Must within 0-1
#' @param palette
#' @param bins
#'
#' @return Contain multiple plots
#' @export
#'
#' @examples
multiplePlotPredictedPro <- function(group, pred, palette="aaas", bins = 10){
  # https://darrendahly.github.io/post/homr/
  df <- data.frame(pred=pred, Class=group)

  waterfall_plot = loonR::plot_waterfall(df$pred-0.5, df$Class)

  df$ClassFactor = factor(group==unique(group)[2], levels = c(FALSE,TRUE))
  df$ClassFactorNumric = as.numeric(df$ClassFactor)-1

  require(ggpubr)

  logistic_curve_fit <-
  ggplot(df, aes(x = pred, y = as.numeric(factor(group==unique(group)[2], levels = c(FALSE,TRUE))) - 1)) +
    geom_jitter(height = 0.1, size =1, alpha = 0.5) +
    geom_smooth(method = "glm", se = FALSE,
                method.args = list(family = "binomial")) +
    theme_minimal() +
    scale_y_continuous(breaks = c(0, 1), labels = unique(group)) +
    ylab("") +
    xlab("Predicted probability")

  probability_distribution <-
    gghistogram(df, x = "pred", fill = "Class", palette = "aaas", rug = T, bins = bins) +
    ylab("Count") + xlab("Predicted probability")



  probability_proportion_distribution <- ggplot(df, aes(x = pred, fill = Class)) +
    geom_histogram(position = "fill", bins = bins) +
    theme_pubr() +
    xlab("Predicted probability") +
    ylab("Proportion") + scale_fill_manual(values = loonR::get.palette.color(palette))


  boxplot <- ggboxplot(df, x = "Class", y = "pred", fill = "Class", palette = palette, add = "jitter") +
    ylab("Predicted probability") + xlab("")


  res = list(waterfall_plot = waterfall_plot,
             logistic_curve_fit = logistic_curve_fit,
             probability_distribution = probability_distribution,
             probability_proportion_distribution = probability_proportion_distribution,
             boxplot = boxplot)
  res
}




#' Get best lamda by perform multiple round lasso-sv
#' https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
#' https://stackoverflow.com/questions/62170965/looping-cv-glmnet-and-get-the-best-coefficients
#'
#' @param d.matrix Row is sample
#' @param group
#' @param family Default binomial
#' @param type.measure class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.
#' @param nfolds Default 5
#' @param nreps Default 1000
#' @param scale Default TRUE
#'
#' @return
#' @export
#'
#' @examples
lasso_best_lamda <- function(d.matrix, group, family = "binomial", type.measure = "auc", nfolds = 5, nreps = 100, scale=TRUE, seed = 66 ){

  if(scale){
    d.matrix = scale(d.matrix, center = TRUE, scale = TRUE)
    #data.matrix = data.frame(data.matrix, check.names = T, stringsAsFactors = F)
  }



  library(glmnet)
  X = as.matrix(d.matrix)
  Y = group



  library(foreach)
  library(parallel)

  doParallel::registerDoParallel(cores=10)
  parallel::mcaffinity(c(1:10))


  res = foreach(i = 1:nreps, .combine = rbind, .export = c("glmnet") )%dopar%{
    set.seed(666+i)
    fit = cv.glmnet(x=X, y=Y, nfolds=nfolds, family = family, type.measure = type.measure)
    data.frame(`MSE_mean`=fit$cvm, lambda= fit$lambda, se = fit$cvsd)
  }


  library(dplyr)

  summarized_res = res %>%
    group_by(lambda) %>%
    summarise(MSE=mean(`MSE_mean`),se=mean(se)) %>%
    arrange(desc(lambda))



  idx = which.min(summarized_res$MSE)
  lambda.min = summarized_res$lambda[idx]


  index_1se = with(summarized_res, which(MSE < MSE[idx]+se[idx])[1])
  lambda_1se = summarized_res$lambda[index_1se]


  library(ggplot2)
  p <- ggplot(res,aes(x=log10(lambda),y=MSE_mean)) + stat_summary(fun=mean,size=2,geom="point") +
    geom_vline(xintercept=log10(c(lambda.min,lambda_1se))) + theme_bw()


  fit = glmnet(x=X, y=Y, family = family, type.measure = type.measure, lambda = lambda.min)

  feature.coef = coef(fit, s=lambda.min)
  feature.coef = data.frame(name = feature.coef@Dimnames[[1]][feature.coef@i + 1], coefficient = feature.coef@x)

  library(dplyr)
  PredictedValue = predict(fit,d.matrix)[,1]
  d.matrix = data.frame(d.matrix, check.names = F)
  d.matrix$PredictedValue = PredictedValue


  res = list(lambda.min = lambda.min, lambda_1se = lambda_1se, res = res,
             summarized_res = summarized_res, plot = p,
             fit = fit, feature.coef = feature.coef,
             data = d.matrix
             )
  res
}



#' Perform multiple rounds differential analysis to select stable feature
#'
#' @param log.df Row is gene, col is samples
#' @param label Normal first
#' @param folds 5
#' @param seed 666
#' @param n 1000
#' @param cores 50
#' @param AUC.cut.off 0.8
#' @param FoldC.cut.off 1
#'
#' @return
#' @export
#'
#' @examples
limma.differential.cv.selection <- function(log.df, label,
                  folds = 5, seed = 666, n = 1000, cores = 50, AUC.cut.off = 0.8, FoldC.cut.off = 1){



  discovery.design <- data.frame(Group = label,
                                 Sample = colnames(log.df),
                                 stringsAsFactors = F)


  library(foreach)
  library(parallel)
  library(doParallel)
  registerDoParallel(cores=cores)
  parallel::mcaffinity(c(1:cores)) # limit cores to use


  # cross validation
  cv.raw.res <- foreach::foreach(i=1:n, .combine = rbind, .packages = c("dplyr")) %dopar% {
    set.seed(i+seed)
    train.design <- discovery.design %>% group_by(Group) %>% sample_frac((folds-1)/folds)
    validation.design <- anti_join(discovery.design, train.design, by = 'Sample')

    train.df <- log.df[, train.design$Sample ]
    train.diff.res <- loonR::limma_differential(train.df, train.design$Group, pre.filter = 1)

    train.diff.res
  }

  cv.discover.raw.res <- cv.raw.res %>% filter(logFC > FoldC.cut.off & AUC > AUC.cut.off & adj.P.Val < 0.05)

  # identify candidates frequency
  candidate.occurence <- data.frame( unlist( table(cv.discover.raw.res$REF) ), stringsAsFactors = FALSE )
  colnames(candidate.occurence) <- c("Name","Freq")

  # candidate mean in all rounds
  candidate.auc <- aggregate(cv.raw.res %>% select(logFC, AUC, AveExpr), by = list(cv.raw.res$REF), FUN = mean)

  names(candidate.auc) <- c("Name", "Mean (LogFC)", "Mean AUC","Mean (log2 CPM)")
  summarize.raw <-candidate.auc
  candidate.auc <- candidate.auc %>% filter(Name %in% candidate.occurence$Name)


  candidate.raw.res <- full_join(candidate.occurence, candidate.auc, by="Name")
  candidate.raw.res$Freq <- round(candidate.raw.res$Freq/n, 2)

  tmp.plot <- cv.raw.res %>% filter(REF %in% candidate.raw.res$Name) %>% select(REF, AUC)
  candidate.auc.in.cv.boxplot <- ggboxplot(tmp.plot, x = "REF", y = "AUC", xlab = "") + rotate_x_text(45)


  res = list(Raw = cv.raw.res,
             Candidate.Raw = cv.discover.raw.res,
             Summarized.Candidate = candidate.raw.res,
             Summarized.Raw = summarize.raw,
             Plot = candidate.auc.in.cv.boxplot)
  res


}



#' Get Youden index value
#'
#' @param pred
#' @param label
#'
#' @return
#' @export
#'
#' @examples get.YoudenIndex(risk.scores, labels)
get.YoudenIndex <- function(pred, label){


  library(pROC)
  library(caret)

  #pred <- df_plot[tmp.ind,]$RS

  input="threshold"

  best.cutoff <- coords(roc(label, pred), "best", transpose = TRUE, input="threshold", best.method="youden")
  if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
  }else{
      best.cutoff <- best.cutoff[c("threshold")]
  }

  best.cutoff
}


#' Obtain AUC value
#'
#' @param pred Predicted score or probability
#' @param label label/class/group
#' @param raw Default TRUE, return raw auc result with description
#' @param direction <: (controls < t <= cases) ######  >: (controls > t >= cases)
#' @return
#' @export
#'
#' @examples
#' data("LIRI")
#' loonR::get.AUC(LIRI$ANLN, LIRI$status)
#'
get.AUC <- function(pred, label, raw=TRUE, direction = "auto"){

  roc_obj <- pROC::roc(label, pred, quiet=TRUE, direction = direction)
  if(raw){
    pROC::auc(roc_obj)
  }else{
    auc = stringr::str_remove( pROC::auc(roc_obj), "Area under the curve: " )
    round( as.numeric(auc), digits = 3)
  }


}



#' Split sample by group into traning and validation set
#'
#' @param sample Vector
#' @param group Vector
#' @param seed Default 666
#' @param fraction Default 0.5
#'
#' @return list(Train=train.design, Validation=validation.design)
#' @export
#'
#' @examples
splitSampleByGroup <- function(sample=NULL, group=NULL, seed = 666, fraction = 0.5){

  if(is.null(sample) | is.null(group)){
    stop("Sample or group should not be NA")
  }else{
    df = data.frame(Sample= sample, Group = group)
  }

  set.seed(seed)

  train.design <- df %>% group_by(Group) %>% sample_frac(fraction)
  validation.design <- anti_join(df, train.design, by = 'Sample')

  list(Train=train.design,
       Validation=validation.design)

}



#' split data frame by group
#'
#' @param data.frame row is sample
#' @param group
#' @param seed 666
#' @param fraction 0.5
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' res <- splitDataByGroup(data.frame = LIRI, group = LIRI$status)
splitDataByGroup <- function(data.frame=NULL, group=NULL, seed = 666, fraction = 0.5){
  if(is.null(data.frame) | is.null(group)){
    stop("data.frame or group should not be NA")
  }
  set.seed(seed)
  train_idx = caret::createDataPartition(y = group,p = fraction,list = FALSE)

  train.df <- data.frame[train_idx,]
  train.group <- group[train_idx]

  validation.df <- data.frame[-train_idx,]
  validation.group <- group[-train_idx]

  list(Train.df=data.frame(train.df, check.names = F),
       Validation.df=data.frame(validation.df, check.names = F),
       Train.group=train.group,
       Validation.group=validation.group)

}



#' confidence interval
#'
#' @param vector
#' @param interval
#'
#' @return
#' @export
#'
#' @examples
#' https://stackoverflow.com/questions/48612153/how-to-calculate-confidence-intervals-for-a-vector
#'
#' vector <- c(12, 17, 24, 35, 23, 34, 56)
#' confidence_interval(vector, 0.90)
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}







#' Forest plot
#'
#' @param tabletext data.frame. Column names will be used as header.
#' @param estimate.data Estimate, upper CI, lower CI. Corresponded with tabletext
#' @param appendHeader If provide a vector, top header will be added
#' @param specify.summary Specify the row index. Row will be bold character.
#' @param clipping Range to cut. Will use arrow at two cutting points
#' @param graph.pos Column index, or "left" or "right". Column index to put in right section
#' @param xlab xlab
#' @param xlog If in log-format. If TRUE, vectical line in 0, otherwise a vectical line at 0
#' @param xticks Default c(0.1, 0.5, 1, 2, 3, 4)
#' @param x.font.size x-axis label font size.
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' or.res <- loonR::univariate_or(LIRI[,c(1,3,4)],LIRI$status)
#'
#' estimate.data = or.res[,c(2,3,4)]
#' text.data = data.frame(Variate = or.res$Variate, OR = or.res$OR)
#'
#' loonR::forest_plot(text.data, estimate.data, graph.pos = 2, specify.summary = 1)
forest_plot <- function(tabletext, estimate.data, appendHeader = NULL, specify.summary = NULL,
                        clipping = c(0.1, 4), graph.pos = "right", xlab = "", xlog = TRUE,
                        xticks = c( 0.1, 0.5, seq(1,4,1) ), x.font.size = 1  ){

  if(!require(forestplot)){
    BiocManager::install("forestplot")
  }

  summary.label = c(TRUE, rep(FALSE,nrow(tabletext)))
  tabletext = rbind( colnames(tabletext), tabletext)
  estimate.data <- rbind(NA, estimate.data)
  hline = c(1,2,nrow(tabletext)+1)


  if(!is.null(specify.summary)){
    summary.label[(specify.summary+1)] = TRUE
  }


  if(!is.null(appendHeader)){
    summary.label = c(TRUE, summary.label)
    tabletext = rbind( appendHeader, tabletext)
    estimate.data <- rbind(NA, estimate.data)
    hline[2] = 3
  }

  hrzl_lines_width = 1:ncol(tabletext)
  if(graph.pos=="left"){
    hrzl_lines_width = hrzl_lines_width + 1
  }else if( graph.pos!="left" & graph.pos!="right" ){
    # add plot index to the last column
    hrzl_lines_width = c( hrzl_lines_width, ncol(tabletext) + 1 )
    # remove plot ined
    hrzl_lines_width= hrzl_lines_width[ hrzl_lines_width != (graph.pos) ]

  }


  hrzl_lines = list(
    gpar(lty = 1, col = "black", lwd = 1.5),
    gpar(lty = 1, col = "black", columns = hrzl_lines_width, lwd = 1.5 ),
    gpar(lty = 1, col = "black", columns = hrzl_lines_width, lwd = 1.5 )
  )

  names(hrzl_lines) = hline

  forestplot(tabletext, estimate.data,
             hrzl_lines = hrzl_lines,
             new_page = TRUE, vertices = TRUE,
             is.summary = summary.label,
             graph.pos = graph.pos,
             clip = clipping,
             xlog = TRUE, boxsize = 0.15,
             xticks = xticks,
             col = fpColors(box       = "black",
                            line      = "black",
                            summary   = "black",
                            hrz_lines = "black"),
             xlab = xlab,
             txt_gp = fpTxtGp(ticks=gpar(cex = x.font.size))
             )


}



#' Try different mtry and select the best fitted model
#'
#' @param rf.df Row is sample
#' @param group
#' @param ntree Default 500
#' @param seed Default 111
#' @param scale If perform scale
#'
#' @return
#' @export
#'
#' @examples
build.best.random.forest <- function(rf.df, group, ntree = 500, seed=111, scale = TRUE){

  if(scale==TRUE){
    rf.df = loonR::scaleDF(rf.df, byColumn = TRUE)
  }


  library(randomForest)

  min  = 100
  mtry = 0

  for (i in 1:20){

    set.seed(seed)
    rf.classifier <- randomForest(x = rf.df,
                                  y = as.factor(group),
                                  mtry = i,
                                  ntree = ntree)

    err<-mean(rf.classifier$err.rate)
    print(err)
    if(err<min) {
      min = err
      mtry = i }
  }

  print(min)
  print(mtry)

  set.seed(seed)
  rf.classifier <- randomForest(x = rf.df,
                                y = as.factor(group),
                                mtry = mtry,
                                ntree = ntree)
  rf.classifier

}



#' Build parametric survival model
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param seed Default 666
#' @param scale
#' @param time.point Default 36.  Calculate the prob at a specific time point
#'
#' @return
#' @export
#'
#' @examples
#' # psm (parametric survival model) uses a survival model based on functions and their parameters.
#' # cph (Cox Proportional Hazards Model and Extensions) is using the cox model (and the Anderson-Gill model) which is based on the hazard functions.
#'
#' data(LIRI)
#' res = build.psm.regression.model(LIRI[,3:5],LIRI$status, LIRI$time)
build.psm.regression.model <- function(d.frame, status, time, seed=666, scale = TRUE, time.point = 36){

  if(scale){
    d.frame = scale(d.frame, center = TRUE, scale = TRUE)
    d.frame = data.frame(d.frame, stringsAsFactors = FALSE, check.names = F)
  }


  covariates <- colnames(d.frame)

  df <- data.frame(d.frame,
                   Time=time,
                   Status =status,
                   check.names = F )

  formula <- as.formula( paste0("Surv(time, status) ~ `",
                                paste0(covariates, sep='', collapse = '` + `'),
                                '`',
                                sep='', collapse = "")
  )


  library(rms)
  library(survival)
  set.seed(seed)
  # 建立参数性生存分析模型
  sur_model <- psm(formula,
                   data = df,
                   dist = "weibull") # weibull分布

  timepoint.prob = surv(time.point, lp = sur_model$linear.predictors)

  res = list(model = sur_model, data = df, timepoint.prob = timepoint.prob)
  res
}



#' Nomogram plot by rms
#'
#' @param fit rms model
#' @param data The data used to build the model
#' @param fun.list an optional function to transform the linear predictors, and to plot on another axis. If more than one transformation is plotted, put them in a list, e.g. list(function(x) x/2, function(x) 2*x). Any function values equal to NA will be ignored.
#' @param lp 线性预测If fun.list is NA, lp will be TRUE. Set to FALSE to suppress creation of an axis for scoring X beta
#'
#' @return
#' @export
#'
#' @examples
#' # Logistic model
#' data(LIRI)
#' res = loonR::build.logistic.model(LIRI[,3:5], LIRI$status, rms = T, scale = F)
#' f1 = list(Risk = loonR::logit2prob)
#' loonR::nomogram.plot(res$model, res$data, lp = T)
#' loonR::nomogram.plot(res$model, res$data, f1, lp =F)
#'
#' # Survial model
#' res = build.psm.regression.model(LIRI[,3:5],LIRI$status, LIRI$time, scale = F)
#'
#' surv <- Survival(res$model) # This would also work if f was from cph
#' surv_100 <- function(x) surv(100, lp = x)
#' surv_300 <- function(x) surv(300, lp = x)
#'
#' med <- Quantile(res$model)
#' med_f <- function(x) med(lp=x)
#'
#' f.list=list(
#'   `Median Survival Time`= med_f,
#'   `Probability of 100-day Survival`=surv_100,
#'   `Probability of 300-day Survival`=surv_300
#' )
#' nomogram.plot(res$model, res$data, fun.list = f.list, lp =F)
nomogram.plot <- function(fit=NULL, data = NULL, fun.list = NA, lp =F){

  message("If the warining raise: ",
          "Pls run manualy outsite the function:
  ddist <- datadist(res$data)
  options(datadist = \"ddist\")  ")

  if( is.null(fit) | is.null(data) ){
    stop("Please set all the parameter correctly")
  }

  library(rms)
  # 转化为datadist
  ddist <- datadist(data)
  options(datadist = "ddist")

  if(is.na(fun.list)){
    nomo <- rms::nomogram(fit, lp = lp)
  }else{
    nomo <- rms::nomogram(fit, fun = fun.list, lp = lp)
  }


  plot(nomo, col.grid =c("tomato", "darkcyan"))

}



#' Compare two ROC curves
#'
#' @param risk1
#' @param lab1
#' @param risk2
#' @param lab2
#' @param method the method to use, either “delong”, “bootstrap” or “venkatraman”. The first letter is sufficient.
#' @param alternative "two.sided", "less", "greater"
#'
#' @return
#' @export
#'
#' @examples
compareROC <- function(risk1, lab1, risk2, lab2, method = c("delong", "bootstrap", "venkatraman") , alternative = c("two.sided", "less", "greater")){

  method  = match.arg(method)
  alternative = match.arg(alternative)

  cat("\n===============\n\nUsed method: ",method,"\n Used alternative: ",alternative, "\n\n\n")

  pROC::roc.test( pROC::roc(lab1, risk1),
                  pROC::roc(lab2, risk2),
                  method=method,
                  alternative = alternative)



}



#' Generate gene pair value
#'
#' @param value
#' @param names
#' @param sep Pair seperation
#'
#' @return
#' @export
#'
#' @examples
#'
#' set.seed(111)
#' v = rnorm(10)
#' n = paste0("ID", 1:10)
#' names(v) = n
#' loonR::generatePairValue(v,n)
#'
generatePairValue <- function(value=NULL, names = NULL, sep ="_"){

  if(is.null(value)){
      stop("Pls set parameter value")
  }
  if(is.null(names)){
    names = names(value)
  }else{
    if(sum(names(value) == names) != length(names) ){
      stop("Names should be the same")
    }
  }

  combs = loonR::generateCombinations(names, size = 2, vector = T, sep=sep)
  combs = combs[,1]

  value.rank = rank(-1*value)

  library(foreach)

  # based on value
  diff = foreach(comb = combs, .combine = rbind) %do%{
    gene1 = stringr::str_split(comb,sep)[[1]][1]
    gene2 = stringr::str_split(comb,sep)[[1]][2]

    c(gene1, gene2,
      as.numeric(c(
      value[gene1], value[gene2], value.rank[gene1], value.rank[gene2],
      value[gene1] - value[gene2],
      value.rank[gene2] - value.rank[gene1],
      value.rank[gene2] > value.rank[gene1]))
      )
  }
  colnames(diff) = c("G1","G2","Value1","Value2","Rank1","Rank2",
                           "Value.diff",
                           "Rank.diff",
                           "Binary")
  rownames(diff) = combs
  diff = data.frame(diff, stringsAsFactors = F, check.names = F)

  diff[,3:ncol(diff)] = loonR::convertDfToNumeric(diff[,3:ncol(diff)])
  diff = data.frame(diff, stringsAsFactors = F, check.names = F)

  diff
}

#' Generate gene pair value data.frame
#'
#' @param df row is gene, column is sample
#'
#' @return
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' df = matrix(rnorm(100*10, 1, .5), ncol=10)
#' colnames(df) = paste0("ID",1:10)
#' rownames(df) = paste0("Gene",1:100)
#' res = loonR::generateGenePairValueDf(df)
#' head(res$Binary.res)
generateGenePairValueDf <- function(df){

  df = as.matrix(df)

  sample.gene.pair.value.list = list()

  for(sample in colnames(df)){
    sample.value = df[,sample]
    sample.paired.value = loonR::generatePairValue(sample.value)
    sample.gene.pair.value.list[[sample]]=sample.paired.value
  }
  names(sample.gene.pair.value.list) = colnames(df)

  # make sure row names is the same
  rname = rownames(sample.gene.pair.value.list[[1]])

  library(foreach)
  value.diff.res = foreach(sample=colnames(df), .combine = 'cbind') %do%{
     t = sample.gene.pair.value.list[[sample]][rname,"Value.diff"]
     as.vector(t)
  }
  colnames(value.diff.res) = colnames(df)
  rownames(value.diff.res) = rname
  value.diff.res = loonR::convertDfToNumeric(value.diff.res)


  rank.diff.res = foreach(sample=colnames(df), .combine = 'cbind') %do%{
    t = sample.gene.pair.value.list[[sample]][rname,"Rank.diff"]
    as.vector(t)
  }
  colnames(rank.diff.res) = colnames(df)
  rownames(rank.diff.res) = rname
  rank.diff.res = loonR::convertDfToNumeric(rank.diff.res)


  binary.res = foreach(sample=colnames(df), .combine = cbind) %do%{
    t = sample.gene.pair.value.list[[sample]][rname,"Binary"]
    as.vector(t)
  }
  colnames(binary.res) = colnames(df)
  rownames(binary.res) = rname
  binary.res = loonR::convertDfToNumeric(binary.res)

  res = list(
    RawSampleDataList=sample.gene.pair.value.list,
    Binary.res = binary.res,
    Value.diff.res = value.diff.res,
    Rank.diff.res = rank.diff.res
  )

  res

}















