
#' Title
#'
#' @param df data.frame. Row is sample, col is gene
#' @param label TRUE/FALSE vector
#' @param k Folds
#' @param times Default: 1. the number of pieces was use to predict. E.g. times = 5, will run 5 times, each time use a different piece to predict.
#'
#' @return
#'
#' @examples
getOneRoundCVRes <- function(df, label, k, seed = 1, times = 1){

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

    pre.res <- predict(glm.fit, df[piece.ind, ])

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
  require(dplyr)

  #cl = parallel::makeCluster(40)
  #doParallel::registerDoParallel(cl)

  cv.res <- foreach::foreach(i= 1:n, .combine = rbind, .packages="foreach", .export=c("getOneRoundCVRes")) %do% {
  #res <- foreach::foreach(i= 1:n, .combine = rbind) %do% {

    getOneRoundCVRes(df, label, k, seed=i)

  }

  cv.res.mean <- cv.res %>% group_by(Name, Label) %>% summarize(Mean = mean(Logit))

  #stopCluster(cl)
  cv.res.mean <- cv.res.mean[match(row.names(df), cv.res.mean$Name), ]
  cv.res.mean
}





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


#' Leave one out cross validation
#'
#' @param df
#' @param label
#' @param seed Default 999
#'
#' @return
#' @export
#'
#' @examples
loo.cv <- function(df, label, seed=999){

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





#' Title Get model performace, sensitivity, sepcificity and others
#'
#' @param pred Predicted score
#' @param labels True label
#'
#' @return
#' @export
#'
#' @examples get_performance(risk.probability, label)
get_performance <- function(pred, labels, best.cutoff =NA){  #  x="best", input = "threshold"

  library(pROC)
  library(caret)

  #pred <- df_plot[tmp.ind,]$RS
  #lables <- as.factor(tmp.label)
  # pred=kumamoto_rs_cv
  # labels= kumamoto$label
  # x="best"
  # input="threshold"

  input="threshold"

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

#' Variate logistic analysis
#' Row: sample, Column: gene expression
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param label True label
#'
#' @return
#' @export
#'
#' @examples multivariate_or(data.frame, label)
multivariate_or <- function(d.frame, label){

  res <- glm(Event ~ . , data = data.frame(d.frame, Event=label, check.names = FALSE), family=binomial(logit))
  #exp( coef(res) )
  #summary(res)

  res <- data.frame(
    exp( cbind(coef(res), confint(res) )  ) ,
    Estimate = as.vector(summary(res)$coefficients[,1] ),
    Pr = as.vector(summary(res)$coefficients[,4] )
  )

  res <- data.frame(Vairate=row.names(res), round(res,3) )
  res[,2:ncol(res)] <- data.frame(lapply(res[,2:ncol(res)],as.numeric))

  colnames(res) <- c("Variate","OR", "2.5%", "97.5%", "Estimate" , "Pr")
  res

}



#' Variate logistic analysis
#' Row: sample, Column: gene expression
#' score E.g.: Gene or miRNA expression, or risk score
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param label True Sample label
#'
#' @return c(OR, 2.5% CI, 97.5% CI)
#' @export
#'
#' @examples univariate_or(risk.df, label)
univariate_or <- function(d.frame, label){

  library(foreach)

  all.res <- foreach(i=1:ncol(d.frame), .combine = rbind) %do%{
    #   for(i in 1:ncol(d.frame) ){
    res <- glm(Event ~ . , data = data.frame(Score = d.frame[,i], Event=label), family=binomial(logit))
    #exp( coef(res) )
    #summary(res)

    # column name
    c.name <- colnames(d.frame)[i]


    res <- c( exp( cbind(coef(res), confint(res) )  )[2,] , #  c("OR", "2.5 %", "97.5 %")
              summary(res)$coefficients[2,1], # "Estimate"
              summary(res)$coefficients[2,1]  # "Pr"
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


#' Plot heatmap with log 2 fold change and risk probability information
#'
#' @param heatmap.df Row: miRNA, Column: Sample
#' @param label A factor with labels TURE/FALSE
#' @param risk.pro
#' @param lgfold corresponded to row name sequence
#' @param group.name Default "Cancer"
#' @param scale Default "TRUE"
#' @param ylim = c(0, 1), risk score range.
#' @param show.lgfold = TRUE。 Whether to show right panel.
#' @param show.risk.pro
#' @param bar.name
#' @param height
#' @param show_column_names Default False
#' @param cluster_rows
#' @param z.score.cutoff Default 2
#'
#' @return A heatmap plot by complex heatmap
#' @export
#'
#' @examples heatmap.with.lgfold.riskpro(data.tmp[candi,],label, logfd,  risk.pro)
heatmap.with.lgfold.riskpro <- function(heatmap.df, label, risk.pro, lgfold=NA, scale=TRUE, group.name="Cancer", bar.name = "Log2FC", ylim = c(0, 1),
                                        show.lgfold = TRUE, show.risk.pro = TRUE, height = 5, show_column_names = FALSE, cluster_rows = FALSE, z.score.cutoff = 2 ){

  if (anyNA(lgfold)){
    show.lgfold = FALSE
    lgfold = replicate(nrow(heatmap.df),1)
  }

  label = factor(label, levels = unique(label))

  if(scale){
    heatmap.df = t(scale(t(heatmap.df)))
    heatmap.df[ heatmap.df > z.score.cutoff] <- z.score.cutoff
    heatmap.df[ heatmap.df < -z.score.cutoff] <- -z.score.cutoff
  }

  library(ComplexHeatmap)
  if(show.risk.pro){
  # 根据risk score排序
    heatmap.df <- cbind(
      heatmap.df[, label==levels(label)[1] ][, order(risk.pro[ label==levels(label)[1] ] ) ],
      heatmap.df[, label==levels(label)[2] ][, order(risk.pro[ label==levels(label)[2] ] ) ]  )

    # 对riskscore排序
    risk.pro <- c(
      risk.pro[label==levels(label)[1] ][order(risk.pro[label==levels(label)[1] ] ) ],
      risk.pro[label==levels(label)[2] ][order(risk.pro[label==levels(label)[2] ] ) ]  )
    }


  # heatmap和barplot一起画

  #row_ha = rowAnnotation( assign(eval(bar.name), anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))   )
  row_ha = rowAnnotation( Log2FC = anno_barplot(lgfold, gp = gpar(fill = "black",col="black")))

  label = factor(label)
  Tumor = loonR::get.palette.color("jama_classic", n = 2)
  names(Tumor) = levels(label)

  # rename annotation names
  annotation <- data.frame(Tmp = label[ c(which(label==levels(label)[1]),
                                          which(label==levels(label)[2]) )
                                        ]
                           )
  colnames(annotation) <- group.name

  ann_colors = list(Tmp = Tumor)
  names(ann_colors) = group.name


  if(show.risk.pro){
  ha = HeatmapAnnotation(df = annotation,
                         col = ann_colors,
                         Risk = anno_points(risk.pro, pch = 16, size = unit(1, "mm"),
                                            gp = gpar(col = "black"),
                                            ylim = ylim,
                                            axis_param = list( side = "left", at = ylim, labels = as.character(ylim) )
                                            )
                         )
  }else{
    ha = HeatmapAnnotation(df = annotation,
                           col = ann_colors    )
  }

  #
  if(show.lgfold){

    Heatmap(heatmap.df, col = c("#0c3e74","#77a8cd","white","#d86652","#7e0821"),
            name = " ", cluster_rows = cluster_rows, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = show_column_names, height = unit(height, "cm"),
            top_annotation = ha,
            right_annotation = row_ha  )

  }else{

    Heatmap(heatmap.df, col = c("#0c3e74","#77a8cd","white","#d86652","#7e0821"),
            name = " ", cluster_rows = cluster_rows, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = show_column_names, height = unit(height, "cm"),
            top_annotation = ha  )

  }


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
  if (LogC) { exp.df[which(exp.df <= 0)] <- NaN
  exp.df <- log2(exp.df) }
  rm(qx, LogC)

  result <- list(expression = exp.df,
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
#'
#' @return
#' @export
#'
#' @examples loonR::plot_waterfall(average.riskscore$Mean, average.riskscore$Label, xlab = "Risk probability")
plot_waterfall <- function(risk.score, label, xlab = "Risk probability", palette = "jco", title = "", yticks.labl = c(0,0.25,0.5,0.75,1)){

    library(ggpubr)
    risk.score = risk.score
    idx = order(risk.score)
    risk.score <- risk.score[idx]
    label = label[idx]


    tmp.df <- data.frame(Risk=risk.score, Class = label, ID=1:length(risk.score) )
    colnames(tmp.df)[1] <- xlab

    p <- ggbarplot(tmp.df, x = "ID", y = xlab, xlab = "", color = "Class", fill = "Class", palette = palette, legend = "right", title = title) +
      rremove("x.axis") + rremove("x.text") + rremove("x.ticks")

    if(!anyNA(yticks.labl)){
      p <- p + scale_y_continuous(labels = yticks.labl)
    }

   p
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
#'
#' @return
#' @export
#'
#' @examples
lasso_best_lamda <- function(d.matrix, group, family = "binomial", type.measure = "auc", nfolds = 5, nreps = 1000 ){

  library(glmnet)
  X = as.matrix(d.matrix)
  Y = group

  library(foreach)
  library(parallel)
  registerDoParallel(cores=10)
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


  res = list(lambda.min = lambda.min, lambda_1se = lambda_1se, res = res,
             summarized_res = summarized_res, plot = p,
             fit = fit, feature.coef = feature.coef
             )
  res
}


#' Feature selection by lasso
#'
#' @param data.matrix Row is sample
#' @param label
#' @param folds
#' @param seed
#' @param family Default binomial
#' @param type.measure class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.
#' @param s Defalut is lambda.min. User can specify
#'
#' @return
#' @export
#'
#' @examples
lasso.select.feature <- function(data.matrix, label, folds = 5, seed = 666,
                                    family = "binomial", type.measure = "auc" , s = NA){
  library(foreach)
  library(dplyr)


  set.seed(seed)
  require(caret)
  cvfit = cv.glmnet(as.matrix(data.matrix), label, nfolds = folds,
                      family = family, type.measure = type.measure)

  if (is.null(s)) {
    s = lambda.min
  }else{
    s = cvfit$lambda.min
  }


  feature.coef = coef(cvfit, s = s)
  feature.coef = data.frame(name = feature.coef@Dimnames[[1]][feature.coef@i + 1], coefficient = feature.coef@x)

  feature.coef = feature.coef[-c(1), ] # remove Intercept
  feature.coef$auc = apply(data.matrix[,feature.coef$name], 2,function(x){
    suppressMessages(roc <- pROC::roc(label, x)  )
    roc$auc
  })

  feature.coef


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
#' @param type.measure class, auc, deviance, mae. “deviance” uses actual deviance. “mae” uses mean absolute error. “class” gives misclassification error. “auc” (for two-class logistic regression ONLY) gives area under the ROC curve.
#'
#' @return
#' @export
#'
#' @examples
lasso.cv.select.feature <- function(data.matrix, label, folds = 5, seed = 666, n = 100,
                                 family = "binomial", type.measure = "class" , cores = 50){
  library(foreach)
  library(dplyr)
  library(parallel, doParallel)
  library(glmnet)

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


