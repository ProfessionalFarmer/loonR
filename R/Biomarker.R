
#' Title
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


#' Title
#'
#' @param df data.frame. Row is sample, col is gene
#' @param label True label
#' @param k Folds
#' @param n Repeat times
#' @param scale TRUE
#' @param type "response" or "link"
#'
#' @return
#' @export
#'
#' @examples
#' cv.res <- loonR::cross.validation(miR.df[, cf.candidates],
#'             group == "Cancer",
#'             k = 10, n = 100)
cross.validation <- function(df = '', label = '', k = 5, n = 100, scale=TRUE, type = "response"){

  # df = up.mirs.exp
  # label = mir.sample.group
  # k = 5
  # n = 100

  if(df == '' | label == ''){
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

    getOneRoundCVRes(df, label, k, seed=i, type = type)

  }

  cv.res.mean <- cv.res %>% group_by(Name, Label) %>% summarize(Mean = mean(Logit))

  #stopCluster(cl)
  cv.res.mean <- cv.res.mean[match(row.names(df), cv.res.mean$Name), ]
  cv.res.mean
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
  }else{
    glm.fit <- rms::lrm(label ~ .,lg.df, x = TRUE, y = TRUE)
    res <- list(model=glm.fit)
  }

  lg.df$risk.score <- predict(glm.fit, lg.df)
  p <- loonR::roc_with_ci(lg.df$label, lg.df$risk.score, ci = FALSE)

  res$data = lg.df
  res$ROC = p

  res

}

#' Build co regression model
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param status
#' @param time OS, DFS, RFS et al.....
#' @param seed Default 666
#'
#' @return
#' @export
#'
#' @examples
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
  cox.fit

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

  list(model=rf.fit,
       predicted.label=rf.fit$predicted,
       importance=rf.fit$importance,
       confusion=rf.fit$confusion,
       votes=data.frame(rf.fit$votes, check.names = FALSE)
  )


}
# todo elastic regression
# https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r





#' Get confusion matrix
#'
#' @param groups True label
#' @param rs Predicted score
#' @param cancer Cancer Name
#' @param best.cutoff If not set, use youden index instead
#'
#' @return A data.frame object
#' @export
#'
#' @examples confusion_matrix(label,risk.probability,cancer = "ESCC")
#'
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



#' Title Get model performace, sensitivity, sepcificity and others
#'
#' @param pred Predicted score
#' @param labels True label
#' @param best.cutoff
#' @param digit 2 for 0.01, 3 for 0.001
#'
#' @return
#' @export
#'
#' @examples get_performance(risk.probability, label)
get_performance <- function(pred, labels, best.cutoff =NA, digit = 2){  #  x="best", input = "threshold"

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
  auc <- pROC::ci(roc.obj,boot.n=2000)[c(2, 1, 3)]
  auc <- round(auc, digit)
  auc <- sprintf(string.format, auc[1], auc[2], auc[3])

  # get performance
  set.seed(100)
  rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv",
            "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "precision", "recall")
  others <- pROC::ci.coords(roc.obj, x = best.cutoff, boot.n = 2000, input = input, ret = rets, best.policy = "random", transpose = TRUE)
  # to be continue

  # ppv Precision  https://www.rdocumentation.org/packages/pROC/versions/1.16.2/topics/coords
  # sensitivity recall
  #f1.score <- (2 * others$precision[1,c("50%")] * others$recall[1,c("50%")]) / (others$precision[1,c("50%")] + others$recall[1,c("50%")])
  #f1.score <- round(f1.score,2)

  res <- c(
    n.c,
    t.c,
    auc[1], # AUC
    sprintf(string.format, others$accuracy[1,c("50%")],  others$accuracy[1,c("2.5%")],  others$accuracy[1,c("97.5%")]), # accuracy
    sprintf(string.format, others$precision[1,c("50%")], others$precision[1,c("2.5%")], others$precision[1,c("97.5%")]),# precision
    sprintf(string.format, others$recall[1,c("50%")],    others$recall[1,c("2.5%")],    others$recall[1,c("97.5%")]),#recall
    sprintf(string.format, others$specificity[1,c("50%")],others$specificity[1,c("2.5%")],others$specificity[1,c("97.5%")]),#specificity
    sprintf(string.format, others$npv[1,c("50%")],others$npv[1,c("2.5%")],others$npv[1,c("97.5%")]),#npv
    or  # Odds ratio
  )
  names(res) <- c("Ncount", "Tcount", "AUC (CI)",  "Accuracy", "Precision", "Recall", "Specificity", "NPV","Odds Ratio")

  # value (conf) -> value  conf
  rname <- names(res)
  value <- sapply(strsplit(res," \\("), function(x) paste(x[1]) )

  tmp <- sapply(strsplit(res," \\("), function(x) paste(x[2]) )
  conf <- sapply(strsplit(tmp,"\\)"), function(x) paste(x[1]) )

  df <- data.frame(Name = rname, value = value, confidence = conf)

  df

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
#' @examples univariate_or(risk.df, label)
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
univariate_cox <- function(d.frame, status, time, scale=TRUE){

  library(foreach)
  library("survival")
  library("survminer")

  if(scale){d.frame = scale(d.frame, center = TRUE, scale = TRUE)}


  all.res <- foreach(i=1:ncol(d.frame), .combine = rbind) %do%{


    res <- coxph( Surv(time, status) ~ Score , data = data.frame(Score = d.frame[,i], Time=time, Status =status, check.names = F ) )
    #exp( coef(res) )
    #summary(res)

    # column name
    c.name <- colnames(d.frame)[i]
    i.coef <- coef(res)
    i.hazzardRation <- exp( coef(res) )
    i.hr.upper = summary(res)$conf.int[4]
    i.hr.lower = summary(res)$conf.int[3]
    i.pvalue = summary(res)$waldtest[3]

    res <- c( i.coef, # coef
              i.hazzardRation, # HR
              i.hr.lower, # lower
              i.hr.upper, # upper
              i.pvalue
    )
    c( Name=c.name, round(res,3))

  }

  all.res <- data.frame(all.res, stringsAsFactors = F)
  row.names(all.res) <- all.res$Name
  colnames(all.res) <- c("Variate", "coefficient" , "HR", "lower 95%", "upper 95%", "Pr")

  all.res[,2:ncol(all.res)] <- data.frame(lapply(all.res[,2:ncol(all.res)],as.numeric))
  all.res

}

#' Title
#'
#' Copy http://www.sthda.com/english/wiki/cox-proportional-hazards-model
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
univariate_cox_sthda <- function(d.frame, status, time, scale=TRUE){

  if(scale){d.frame = scale(d.frame, center = TRUE, scale = TRUE)}


  # Clone frome http://www.sthda.com/english/wiki/cox-proportional-hazards-model
  covariates <- colnames(d.frame)

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
#' @examples roc_with_ci(labels, rs,
#' font = "Arial",
#' palette = "jama_classic",
#' title = "HGD vs Healthy",
#' panel = data.frame(Evebt=labels, data)
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
        obj <- pROC::roc( panel[,c(1)], panel[,p.mem], ci=TRUE, plot=FALSE)

        # get performance
        set.seed(100)
        others <- ci.coords(obj, x="best", input="threshold", best.policy = "random")

        tmp.res <- c(as.numeric( round(as.vector(others[["specificity"]]),2)  )  )
        tmp.res <- c(tmp.res, as.numeric( round(as.vector(others[["sensitivity"]]),2)  )  )
        as.numeric(tmp.res)
      }
    p.mem.res <- as.data.frame(p.mem.res)

    colnames(p.mem.res) <- c("specificity.low", "specificity.median", "specificity.high", "sensitivity.low", "sensitivity.median", "sensitivity.high")
    row.names(p.mem.res) <- colnames(panel[,-c(1)])


    p + geom_segment(data=p.mem.res, color="#009090",
                     aes(x=1-specificity.high,
                         xend=1-specificity.low,
                         ymin=(sensitivity.low+sensitivity.high)/2,
                         ymax=(sensitivity.low+sensitivity.high)/2,
                         y=(sensitivity.low+sensitivity.high)/2,
                         yend=(sensitivity.low+sensitivity.high)/2)
    )  +
      geom_segment(data=p.mem.res,  color="#009090",
                   aes(x=(1-specificity.high+1-specificity.low)/2,
                       xend=(1-specificity.high+1-specificity.low)/2,
                       ymin=sensitivity.low,
                       ymax=sensitivity.high,
                       y=sensitivity.low,
                       yend=sensitivity.high)
      ) +
      geom_point(data=p.mem.res, mapping=aes(x=(1-specificity.high+1-specificity.low)/2,
                                             y=(sensitivity.low+sensitivity.high)/2,
                                             ymin=0,ymax=0),
                 size=3.5, shape=16, fill="#009090",color="#009090")

  }
  ## if else

}



#' Plot multiple ROCs in one figure
#'
#' @param scores A list or a data.frame. If list, labels shoule also be a list
#' @param labels
#' @param font
#' @param palette
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
    roclist <- apply(scores, 2, function(x) roc(labels,x) )
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
        obj <- pROC::roc( panel[,c(1)], panel[,p.mem], ci=TRUE, plot=FALSE)

        # get performance
        set.seed(100)
        others <- ci.coords(obj, x="best", input="threshold", best.policy = "random")

        tmp.res <- c(as.numeric( round(as.vector(others[["specificity"]]),2)  )  )
        tmp.res <- c(tmp.res, as.numeric( round(as.vector(others[["sensitivity"]]),2)  )  )
        as.numeric(tmp.res)
      }
    p.mem.res <- as.data.frame(p.mem.res)

    #cat(dim(p.mem.res))


    colnames(p.mem.res) <- c("specificity.low", "specificity.median", "specificity.high", "sensitivity.low", "sensitivity.median", "sensitivity.high")
    row.names(p.mem.res) <- colnames(panel[,-c(1)])
    options(warn = oldw)

    p + geom_segment(data=p.mem.res, color="#009090",
                     aes(x=1-specificity.high,
                         xend=1-specificity.low,
                         ymin=(sensitivity.low+sensitivity.high)/2,
                         ymax=(sensitivity.low+sensitivity.high)/2,
                         y=(sensitivity.low+sensitivity.high)/2,
                         yend=(sensitivity.low+sensitivity.high)/2)
    )  +
      geom_segment(data=p.mem.res,  color="#009090",
                   aes(x=(1-specificity.high+1-specificity.low)/2,
                       xend=(1-specificity.high+1-specificity.low)/2,
                       ymin=sensitivity.low,
                       ymax=sensitivity.high,
                       y=sensitivity.low,
                       yend=sensitivity.high)
      ) +
      geom_point(data=p.mem.res, mapping=aes(x=(1-specificity.high+1-specificity.low)/2,
                                             y=(sensitivity.low+sensitivity.high)/2,
                                             ymin=0,ymax=0),
                 size=3.5, shape=16, fill="#009090",color="#009090")

  }
  ## if else

}




#' Download dataset from GEO by accession ID and platform
#'
#' @param geo.accession.id GEO Accession ID
#' @param platform Platform
#' @param destdir Default tempdir()
#'
#' @return list(expression, phenotype, probe.annotation)
#' @export
#'
#' @examples
download.geo.dataset <- function(geo.accession.id, platform, destdir = tempdir()){

  if (missing("geo.accession.id") | missing("platform")){
    stop("Please provide both geo.accession.id and platform")
  }


  library(GEOquery)

  # load series and platform data from GEO
  gset <- getGEO(geo.accession.id, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = destdir)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  gpl.annotation <- fData(gset)

  # show(gset)

  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  phenotype <- pData(gset)

  # eliminate samples
  exp.df <- exprs(gset)

  # if log2 transform.    VALUE	Normalized log2 data represent microRNA expression levels
  qx <- as.numeric(quantile(exp.df, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    exp.df[which(exp.df <= 0)] <- NaN
    exp.df <- log2(exp.df)
    print("Perform log2 transformation")
  }else{
    print("Note: here not perform log2 transformation")
  }
  rm(qx, LogC)

  result <- list(expression = exp.df,
                 rawExpression = exprs(gset),
                 phenotype  = phenotype,
                 probe.annotation = gpl.annotation)

  return(result)

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
plot_waterfall <- function(risk.score, label, xlab = "Risk probability", palette = "jco", title = "", yticks.labl = c(0,0.25,0.5,0.75,1), sample = NA, rotate.x = 90){

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
#'
#' @return Contain multiple plots
#' @export
#'
#' @examples
#'
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


#' Title
#'
#' @param pred
#' @param label
#'
#' @return
#' @export
#'
#' @examples
get.AUC <- function(pred, label){

  library(pROC)
  library(caret)

  roc_obj <- roc(label, pred, quiet=TRUE)
  auc(roc_obj)

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
splitSampleByGroup <- function(sample=NA, group=NA, seed = 666, fraction = 0.5){

  if(is.na(sample) | is.na(group)){
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
#'
splitDataByGroup <- function(data.frame=NULL, group=NULL, seed = 666, fraction = 0.5){
  if(is.null(data.frame) | is.null(group)){
    stop("data.frame or group should not be NA")
  }
  train_idx = caret::createDataPartition(y = group,p = fraction,list = FALSE)

  train.df <- data.frame[train_idx,]
  train.group <- group[train_idx]

  validation.df <- data.frame[-train_idx,]
  validation.group <- group[-train_idx]

  list(Train.df=train.df, Validation.df=validation.df,
       Train.group=train.group, Validation.group=validation.group)

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







riskCalibrationPlot<-function(x, ...) {
  # Generics function
  UseMethod('riskCalibrationPlot')
}


#' Calibration plot
#'
#' @param group Must be a TRUE/FALSE factor
#' @param pred predicted probability
#' @param rms.method If TRUE, use rms::val.prob function instead
#' @param title
#' @param show.oberved.ci
#' @param bins Number of bins. Default 20
#' @param color Default npg
#' @param show.group whether to plot histogram by group
#' @param ticks.unit 0.25 seq(0, 1, by = 0.25)
#' @param full.range Default TRUE. loess smoothing between 0-1 or first bin to last bin
#'
#' @return
#' @export
#'
#' @examples
#' data(BreastCancer)
#' BreastCancer = BreastCancer[,-c(1)]
#' BreastCancer = na.omit(BreastCancer)
#' m <- glm(Class ~ ., data = BreastCancer, family = binomial)
#' BreastCancer$pred <- predict(m, type = "response")
#' riskCalibrationPlot(factor(BreastCancer$Class=="malignant", levels=c(FALSE, TRUE)),
#'                    BreastCancer$pred)
#'
#' data(LIRI)
#'
#' d1 <- LIRI[,-c(1,5)]
#' m <- glm(status ~ ., data = d1, family = binomial(logit))
#' d1$pred <- predict(m, type = "response")
#' riskCalibrationPlot(factor(LIRI$status), d1$pred)
riskCalibrationPlot.default <- function(group, pred, rms.method = FALSE, title = "Calibration plot", show.oberved.ci = FALSE,  bins = 10, color="npg", show.group = FALSE, ticks.unit=0.25, full.range=TRUE){

  # Thanks reference: https://darrendahly.github.io/post/homr/

  df <- data.frame(pred=pred, RawClass=group)
  if(!is.factor(group)){
    warning("group Must be a TRUE/FALSE factor")
    group = factor(group==unique(group)[2], levels = c(FALSE,TRUE))
  }
  df$Class=group


  if(rms.method){
    print(rms::val.prob(df$pred, as.numeric(df$Class)))
    return()
  }



  require(gridExtra)
  require(dplyr)
  require(ggpubr)



  df.cal.grouped <- arrange(df, pred) %>%
    mutate(bin = ntile(pred, bins)) %>%
    group_by(bin) %>%
    mutate(n = n(),
           bin_pred = mean(pred),
           bin_prob = mean(as.numeric(Class)-1),
           se = sqrt((bin_prob * (1 - bin_prob)) / n),
           ul = bin_prob + 1.96 * se,
           ll = bin_prob - 1.96 * se
    ) # %>% select(n, bin, bin_pred, bin_prob, se, ul ,ll) %>% unique() %>% arrange(bin)
    #ungroup()

  # 不在用，这里显示的CI是根据标准差计算得到
  # if(show.oberved.ci){
  #   p1 = ggplot(df.cal.grouped, aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  #     geom_pointrange(size = 0.5, color = "black") +
  #     scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = ticks.unit)) +
  #     scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = ticks.unit)) +
  #     geom_abline() + # 45 degree line indicating perfect calibration
  #     theme_pubr() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  # }else{
  #  p1 = ggscatter(df.cal.grouped, x="bin_pred", y ="bin_prob" ) +
  #    geom_abline(intercept = 0, slope = 1)
  # }


  p1 = ggscatter(df.cal.grouped, x="bin_pred", y ="bin_prob" ) +
        geom_abline(intercept = 0, slope = 1)
  # loess fit through estimates，并且显示标准差
  if(full.range){
     # This plot the full range 0-1
    p1 = p1 + geom_smooth( aes(x = pred, y = as.numeric(Class)-1 ), level = 0.95,
                           color = "red", se = show.oberved.ci, method = "loess")
  }else{
    #this only plot between bins: the first bin_pred to the last bin_pred
    p1 = p1 + geom_smooth(aes(x = bin_pred, y = bin_prob ), level = 0.95,
                          color = "red", se = show.oberved.ci, method = "loess")
  }


  p1 = ggpar(p1, xlim = c(0,1), ylim = c(0,1), xticks.by =  ticks.unit,
             xlab = "", ylab="Observed proportion", title = title)


  # The distribution plot
  if(show.group){
    p2 <- gghistogram(df, x="pred", fill = "RawClass", bins = bins, rug = TRUE,
                      palette = color, position = "stack") +
      xlab("Predicted probability") + ylab("Count")
  }else{
    p2 <- gghistogram(df, x="pred", fill = "black", bins = bins, rug = TRUE,
                      palette = "black", color="black") +
      xlab("Predicted probability") +  ylab("Count")
  }
  p2 <- ggpar(p2, legend = "none", xlim = c(0,1), xticks.by = ticks.unit)


  plot_row <- cowplot::plot_grid(p1, p2, nrow = 2, ncol = 1, rel_heights = c(2,1))

  print(ResourceSelection::hoslem.test(as.numeric(df$Class)-1, df$pred, g = 10))

  #### rms plot
  # this function can also plot calbiration curve
  library(rms)
  val.prob(df$pred, as.numeric(df$Class) )
  ####

  plot_row
}




#' Calibration plot by rms.calibrate
#'
#' @param rms.model A model built by rms package
#' @param cox If is a cox model
#'
#' @return
#' @export
#'
#' @examples
#' #' data(LIRI)
#' m=loonR::build.logistic.model(LIRI[,c(3,4)],LIRI$status, rms = T)
#' m=m$model
#' riskCalibrationPlot(m)
riskCalibrationPlot.lrm <- function(rms.model, cox=FALSE){
  plot( rms::calibrate(rms.model,
                       method=c("boot"),
                       B=100,
                       smoother="lowess" ),
        xlab = "Predicted probability"
        )
}




#' Decision curve analysis
#'
#' @param data.frame.list row is sample. List should inlcude names
#' @param label
#' @param rms Default FALSE
#' @param validation.df shold include all variable in data.frame.list
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#'
#' d1 <- LIRI[,c(3)]
#' d2 <- LIRI[,c(3,4)]
#' d3 <- LIRI[,c(3,4,5)]
#' data.frame.list = list(d1=d1,d2=d2,d3=d3)
#' decisionCurveAnalysis(data.frame.list, label=LIRI$status)
decisionCurveAnalysis <- function(data.frame.list=NULL, label = NULL, rms=FALSE, validation.df = NULL, palette="aaas"){

  if(!require(ggDCA)){
    install.packages('ggDCA')
  }

  if(is.null(data.frame.list) | is.null(label) ){
     stop("data.frame.list or label should not be null")
  }


  new.df.list = sapply(data.frame.list, function(x){
    d = data.frame(x, label=label)
    #colnames(d) = c(colnames(x),"label")
  })


  n.df <- length(new.df.list)
  if(rms){
    if(n.df==1){
      m1 <- rms::lrm(label~., new.df.list[[1]])
    }else if(n.df==2){
      m1 <- rms::lrm(label~., new.df.list[[1]])
      m2 <- rms::lrm(label~., new.df.list[[2]])
    }else if(n.df==3){
      m1 <- rms::lrm(label~., new.df.list[[1]])
      m2 <- rms::lrm(label~., new.df.list[[2]])
      m3 <- rms::lrm(label~., new.df.list[[3]])
    }else if(n.df==4){
      m1 <- rms::lrm(label~., new.df.list[[1]])
      m2 <- rms::lrm(label~., new.df.list[[2]])
      m3 <- rms::lrm(label~., new.df.list[[3]])
      m4 <- rms::lrm(label~., new.df.list[[4]])
    }
  }else{
    if(n.df==1){
      m1 <- glm(label ~ ., data = new.df.list[[1]], family = binomial(logit))
    }else if(n.df==2){
      m1 <- glm(label ~ ., data = new.df.list[[1]], family = binomial(logit))
      m2 <- glm(label ~ ., data = new.df.list[[2]], family = binomial(logit))
    }else if(n.df==3){
      m1 <- glm(label ~ ., data = new.df.list[[1]], family = binomial(logit))
      m2 <- glm(label ~ ., data = new.df.list[[2]], family = binomial(logit))
      m3 <- glm(label ~ ., data = new.df.list[[3]], family = binomial(logit))
    }else if(n.df==4){
      m1 <- glm(label ~ ., data = new.df.list[[1]], family = binomial(logit))
      m2 <- glm(label ~ ., data = new.df.list[[2]], family = binomial(logit))
      m3 <- glm(label ~ ., data = new.df.list[[3]], family = binomial(logit))
      m4 <- glm(label ~ ., data = new.df.list[[4]], family = binomial(logit))
    }
  }

  if(n.df==1){
    dca_res <- dca(m1,m2, model.names=names(new.df.list), new.data=validation.df)
  }else if(n.df==2){
    dca_res <- dca(m1,m2, model.names=names(new.df.list), new.data=validation.df)
  }else if(n.df==3){
    dca_res <- dca(m1,m2,m3, model.names=names(new.df.list), new.data=validation.df)
  }else if(n.df==4){
    dca_res <- dca(m1,m2,m3,m4, model.names=names(new.df.list), new.data=validation.df)
  }

  ggplot(dca_res,
         linetype=F, #线型
         lwd = 1.2,
         color = c(loonR::get.ggsci.color(palette,n=length(new.df.list)),'black', 'gray')
         )

}




#' Resampling Validation of a Fitted Model's Indexes of Fit
#'
#' @param rms.model
#' @param B Bootstrap times
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' m=loonR::build.logistic.model(LIRI[,c(3,4)],LIRI$status, rms = T)
#' m=m$model
#' rms::validate(m)
rms_validate <- function(rms.model, B = 100){

   # https://blog.csdn.net/fjsd155/article/details/84669331
   v = rms::validate(rms.model, dxy=TRUE, B = B, method="boot", fastbw = FALSE)

   # Get the Dxy
   Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
   orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]

   # The c-statistic according to Dxy=2(c-0.5)
   bias_corrected_c_index  <- abs(Dxy)/2+0.5
   orig_c_index <- abs(orig_Dxy)/2+0.5

   list(
     validate.res = v,
     original.Dxy=orig_Dxy,
     corrected.Dxy=Dxy,
     original.C_index=orig_c_index,
     corrected.C_index=bias_corrected_c_index
        )


}




#' Draw C index across times
#'
#' @param list.cox.model
#' @param df include time and status column names
#' @param palette aaas
#' @param main Default ""
#'
#' @return
#' @export
#'
#' @examples
time_serials_Cindex <- function(list.cox.model, df, palette="aaas", main=""){

  #pec包的cindex()可计算多个模型在多个时间点的C指数
  pk<- cindex(list.cox.model,
              formula=Surv(time,status==0)~1,
              data=df#, eval.times=seq(3,83,1)
              )
  #绘制时间C-index
  plot(pk,
       col=loonR::get.palette.color(palette, n=length(list.cox.model)),#曲线颜色
       xlab="Month",#xy名字
       ylab="C-index",
       #ylim = c(0.4,1),#xy轴范围
       #xlim = c(3,83),
       legend.x=1,     #图例位置
       legend.y=0.6,
       legend.cex=1,    #图例字号
  );title(main = main)
}



