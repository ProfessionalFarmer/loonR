
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





#' Get confusion matrix
#'
#' @param groups
#' @param rs
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
  if( is.na(best.cutoff) ){

    best.cutoff <- coords(roc(groups, risk.pro), "best", transpose = TRUE, input="threshold", best.method="youden")
    if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
    }else{
      best.cutoff <- best.cutoff[c("threshold")]
    }

  }

  results <- confusionMatrix( factor(risk.pro > best.cutoff), factor(groups) , positive = "TRUE" )
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







#' Title Get model performace, sensitivity, sepcificity and others
#'
#' @param pred
#' @param labels
#'
#' @return
#' @export
#'
#' @examples get_performance(risk.probability, label)
get_performance <- function(pred, labels, best.cutoff =NA){  #  x="best", input = "threshold"

  #pred <- df_plot[tmp.ind,]$RS
  #lables <- as.factor(tmp.label)
  # pred=kumamoto_rs_cv
  # labels= kumamoto$label
  # x="best"
  # input="threshold"

  input="threshold"

  if( is.na(best.cutoff) ){

    best.cutoff <- coords(roc(labels, pred), "best", transpose = TRUE, input="threshold", best.method="youden")
    if( inherits(best.cutoff, "matrix")  ){ # 当youden index有重复的时候，取第一个
      best.cutoff <- best.cutoff[c("threshold"),c(1)]
    }else{
      best.cutoff <- best.cutoff[c("threshold")]
    }

  }
  # calculate OR
  results <- confusionMatrix( factor(pred > best.cutoff), factor(labels) , positive = "TRUE" )
  results.tl <- as.table(results)

  or <- vcd::oddsratio(results.tl,log=FALSE)
  or.ci <- confint(vcd::oddsratio(results.tl,log=FALSE), level = 0.95)
  or <- sprintf("%.2f (%.2f-%.2f)",exp(or$coefficients),or.ci[1],or.ci[2]) # Odds ratio

  # count sample
  t.c <- sum(labels==TRUE)
  n.c <- length(labels) - t.c

  # get AUC CI
  set.seed(100)
  roc.obj <- pROC::roc(labels, pred, ci=TRUE, plot=FALSE)
  auc <- pROC::ci(roc.obj,boot.n=2000)[c(2, 1, 3)]
  auc <- round(auc,2)
  auc <- sprintf("%s (%.2f-%.2f)", auc[1], auc[2], auc[3])

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
    sprintf("%.2f (%.2f-%.2f)", others$accuracy[1,c("50%")],  others$accuracy[1,c("2.5%")],  others$accuracy[1,c("97.5%")]), # accuracy
    sprintf("%.2f (%.2f-%.2f)", others$precision[1,c("50%")], others$precision[1,c("2.5%")], others$precision[1,c("97.5%")]),# precision
    sprintf("%.2f (%.2f-%.2f)", others$recall[1,c("50%")],    others$recall[1,c("2.5%")],    others$recall[1,c("97.5%")]),#recall
    sprintf("%.2f (%.2f-%.2f)", others$specificity[1,c("50%")],others$specificity[1,c("2.5%")],others$specificity[1,c("97.5%")]),#specificity
    sprintf("%.2f (%.2f-%.2f)", others$npv[1,c("50%")],others$npv[1,c("2.5%")],others$npv[1,c("97.5%")]),#npv
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

#' Title Row: sample, Column: gene expression
#'
#' @param df
#' @param label
#'
#' @return
#' @export
#'
#' @examples
multivariate_or <- function(df, label){

  res <- glm(Event ~ . , data = as.data.frame(tmp.df, Event=label), family=binomial(logit))
  exp( coef(res) )
  summary(res)

  exp( cbind(coef(res), confint(res) )  )
}







#' Plot heatmap with log 2 fold change and risk probability information
#'
#' @param heatmap.df Row: miRNA, Column: Sample
#' @param label A factor with labels TURE/FALSE
#' @param lgfold corresponded to row name sequence
#' @param risk.procorresponded to sample sequence
#' @param group.name Default "Cancer"
#' @param scale Default "TRUE"
#'
#' @return
#' @export
#'
#' @examples heatmap.with.lgfold.riskpro(data.tmp[candi,],label, logfd,  risk.pro)
heatmap.with.lgfold.riskpro <- function(heatmap.df, label, lgfold, risk.pro, scale=TRUE, group.name="Cancer", ylim = c(0, 1) ){

  label = factor(label)

  if(scale){
    heatmap.df = t(scale(t(heatmap.df)))
    heatmap.df[ heatmap.df > 2] <- 2
    heatmap.df[ heatmap.df < -2] <- -2
  }

  library(ComplexHeatmap)
  # 根据risk score排序
  heatmap.df <- cbind(
    heatmap.df[, label==levels(label)[1] ][, order(risk.pro[ label==levels(label)[1] ] ) ],
    heatmap.df[, label==levels(label)[2] ][, order(risk.pro[ label==levels(label)[2] ] ) ]  )

  # 对riskscore排序
  risk.pro <- c(
    risk.pro[label==levels(label)[1] ][order(risk.pro[label==levels(label)[1] ] ) ],
    risk.pro[label==levels(label)[2] ][order(risk.pro[label==levels(label)[2] ] ) ]  )



  # heatmap和barplot一起画
  row_ha = rowAnnotation( Log2FC = anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))

  label = factor(label)
  Tumor = loonR::get.palette.color("jama_classic", n = 2)
  names(Tumor) = levels(label)

  # rename annotation names
  annotation <- data.frame(Tmp = label)
  colnames(annotation) <- group.name

  ann_colors = list(Tmp = Tumor)
  names(ann_colors) = group.name

  ha = HeatmapAnnotation(df = annotation,
                         col = ann_colors,
                         Risk = anno_points(risk.pro, pch = 16, size = unit(1, "mm"),gp = gpar(col = "black"), ylim = ylim, axis_param = list( side = "left", at = ylim, labels = as.character(ylim) )) )


  # 根据risk probility先组内排序
  Heatmap(heatmap.df,
          name = " ", cluster_rows = FALSE, cluster_columns = FALSE,
          show_row_names = TRUE, show_column_names = FALSE, height = unit(5, "cm"),
          top_annotation = ha,
          right_annotation = row_ha)

}



#' Plot miR's correlation
#'
#' @param df Row: miR expression, Column: Sample
#'
#' @return
#' @export
#'
#' @examples plot_miRCorrelation(data[candi,])
#'
plot_miRCorrelation <- function(df){
  cor.res <- cor(t(df))
  cor.res <- round(cor.res, 3)#保留两位小数

  library(corrplot)#先加载包
  corrplot(cor.res, type = "upper",
           order = "hclust", tl.col = "black", tl.srt = 90, mar=c(0,0,2,0),
           cl.lim = c(-0.5,1), addgrid.col=FALSE, title ="miRs' correlation - TCGA" ) +
  cowplot::theme_cowplot(font_family = "Arial")

}





