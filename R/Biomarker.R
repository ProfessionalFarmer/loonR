
#' Title
#'
#' @param df data.frame. Col is sample, row is gene
#' @param label
#' @param k Folds
#'
#' @return
#'
#' @examples
getOneRoundCVRes <- function(df, label, k, seed = 1){

  set.seed(666+seed)
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

  cl = parallel::makeCluster(40)
  doParallel::registerDoParallel(cl)

  cv.res <- foreach::foreach(i= 1:n, .combine = rbind, .packages="foreach", .export=c("getOneRoundCVRes")) %dopar% {
  #res <- foreach::foreach(i= 1:n, .combine = rbind) %do% {

    getOneRoundCVRes(df, label, k, seed=i)

  }

  cv.res.mean <- cv.res %>% group_by(Name, Label) %>% summarize(Mean = mean(Logit))

  stopCluster(cl)

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
#' @param pred Predicted score
#' @param labels True label
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

#' Variate logistic analysis
#' Row: sample, Column: gene expression
#'
#' @param d.frame Data.frame --- Row: sample, Column: gene expression
#' @param label True label
#'
#' @return
#' @export
#'
#' @examples
multivariate_or <- function(d.frame, label){

  res <- glm(Event ~ . , data = data.frame(d.frame, Event=label), family=binomial(logit))
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
#'
#' @return
#' @export
#'
#' @examples heatmap.with.lgfold.riskpro(data.tmp[candi,],label, logfd,  risk.pro)
heatmap.with.lgfold.riskpro <- function(heatmap.df, label, risk.pro, lgfold=NA, scale=TRUE, group.name="Cancer", ylim = c(0, 1), show.lgfold = TRUE, show.risk.pro = TRUE, height = 5 ){

  if (is.na(lgfold)){
    show.lgfold = FALSE
    lgfold = replicate(nrow(heatmap.df),1)
  }

  label = factor(label)

  if(scale){
    heatmap.df = t(scale(t(heatmap.df)))
    heatmap.df[ heatmap.df > 2] <- 2
    heatmap.df[ heatmap.df < -2] <- -2
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

    Heatmap(heatmap.df,
            name = " ", cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = FALSE, height = unit(height, "cm"),
            top_annotation = ha,
            right_annotation = row_ha  )

  }else{

    Heatmap(heatmap.df,
            name = " ", cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = TRUE, show_column_names = FALSE, height = unit(height, "cm"),
            top_annotation = ha  )

  }


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
roc_with_ci <- function(label, rs, font = "Arial", palette = "jama", legend.pos = c(0.4, 0.2), title = NULL, fontsize = 16, panel = NULL) {

  library(pROC)
  obj = roc(label, rs, ci=TRUE, plot=FALSE)

  # panel的第一列为factor，event，class等

  library(doParallel)
  registerDoParallel(40)
  set.seed(100)

  ciobj <- ci.se(obj, specificities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       se.lower = ciobj[, 1],
                       se.upper = ciobj[, 3])

  library(doParallel)
  registerDoParallel(40)
  set.seed(100)

  ciobj <- ci.sp(obj, sensitivities = seq(0, 1, l = 100), boot.n = 2000, parallel = TRUE)
  dat.ci$y <- as.numeric(rownames(ciobj))
  dat.ci$sp.lower <- ciobj[, 1]
  dat.ci$sp.upper <- ciobj[, 3]

  aucs <- pROC::ci(obj)[c(2, 1, 3)]
  others <- pROC::coords(obj, "best", ret = c("sensitivity", "specificity"), best.policy = "omit")

  annot <- sprintf("AUC %.2f\n(%.2f-%.2f)", aucs[1], aucs[2], aucs[3])
  #annot <- sprintf("AUC %.2f\nSensitivity %.2f\nSpecificity %.2f", aucs[1],others[1,1],others[2,1])


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
#'
#' @return
#' @export
#'
#' @examples multi_roc_with_ci(rss, labels, font = "Arial", palette = "jama")
multi_roc_with_ci <- function(scores, labels, font = "Arial", palette = "jama", legend.pos = c(0.4, 0.2), title = NULL, panel = NULL, color = NULL, ci =TRUE) {

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

    cat(sprintf("%s %.2f (%.2f-%.2f)\n", group_name, auc[1], auc[2], auc[3]) )
    #annot <- c(annot, sprintf("%s %.2f (%.2f-%.2f)", group_name, auc[1], auc[2], auc[3])  )
    annot <- c(annot, sprintf("%.2f (%.2f-%.2f)", auc[1], auc[2], auc[3])  )

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




#' Download RAN expression by TCGAbiolinks
#' @param project E.g. TCGA-LIHC
#' @param remove.Raw If to remove raw data
#' @return list(clinical = clinical, expression = expression.df, group = group)
#' Expression was log2 transformed.
#' Group is a data.frame including label, short name and sample barcode
#' Clinical is a list
#' @export
#'
#' @examples
tcgabiolinks.get.RNA.expression.log2tpm <- function(project, remove.Raw = FALSE){

  fpkm2tpm <- function(fpkm){
    tpm <- exp(log(fpkm) - log(sum(fpkm,na.rm=T)) + log(1e6))
    tpm[which(is.na(tpm))] <- 0
    return(tpm)
  }


  query <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Transcriptome Profiling",
                                  workflow.type = "HTSeq - FPKM", # HTSeq - FPKM-UQ or HTSeq - Counts
                                  data.type = "Gene Expression Quantification",
                                  experimental.strategy = "RNA-Seq",
                                  legacy = FALSE
  )
  TCGAbiolinks::GDCdownload(query)
  project.data <- TCGAbiolinks::GDCprepare(query = query)


  # Expression
  expression.df <- SummarizedExperiment::assay(project.data)
  rm(project.data)
  # convert to tpm
  expression.df <- apply(expression.df, 2, fpkm2tpm)

  normal.sample <- TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                         typesample = "NT")
  tumor.sample <- TCGAquery_SampleTypes(barcode = colnames(expression.df),
                                        typesample = "TP")

  group <- data.frame(
    Label = rep( c("Normal","Tumor"), c(length(normal.sample), length(tumor.sample))  ),
    Short.Name = substr(c(normal.sample,tumor.sample),1,12),
    Barcode = c(normal.sample,tumor.sample),
    stringsAsFactors = FALSE
  )

  # log2 transformation
  expression.df <- log2(expression.df[,c(normal.sample,tumor.sample)])



  ## clinical
  query <- TCGAbiolinks::GDCquery(project = "TCGA-LIHC",
                                  data.category = "Clinical",
                                  data.type = "Clinical Supplement",
                                  data.format = "BCR Biotab")
  TCGAbiolinks::GDCdownload(query)
  clinical <- TCGAbiolinks::GDCprepare(query)


  #clinical <- lihc.clinical$clinical_patient_lihc


  if(remove.Raw){
    file.remove( paste("GDCdata/",project,sep="",collapse = "")  )
  }

  result <- list(clinical = clinical,
                 expression = expression.df,
                 group = group
  )
  return(result)
}



