#' Decision curve analysis. Similar with plot.ntbft
#'
#' @param data.frame.list row is sample. List should inlcude names
#' @param label
#' @param rms Default FALSE
#' @param validation.df should include all variable in data.frame.list
#' @param palette
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#'
#' d1 <- data.frame( LIRI[,c(3)])
#' d2 <- LIRI[,c(3,4)]
#' d3 <- LIRI[,c(3,4,5)]
#' data.frame.list = list(d1=d1,d2=d2,d3=d3)
#' decisionCurveAnalysis(data.frame.list, label=LIRI$status)
#'
decisionCurveAnalysis <- function(data.frame.list=NULL, label = NULL, rms=FALSE, validation.df = NULL, palette="aaas"){
  # https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis
  # https://atm.amegroups.com/article/view/20389/pdf

  # 除了ggDCA dcurves 也可以
  # http://www.danieldsjoberg.com/dcurves/articles/dca.html

  if(!require(ggDCA)){
    install.packages('ggDCA')
  }

  if(is.null(data.frame.list) | is.null(label) ){
    stop("data.frame.list or label should not be null")
  }


  new.df.list = sapply(data.frame.list, function(x){
    d = data.frame(x, label=label)
    colnames(d) = c(colnames(x),"label")
    d
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
  # DCA cox model analysis
  # devtools::install_github('yikeshu0611/ggDCA')
  # library(rms)
  # library(ggDCA)
  #
  # bb<-datadist(train.cox.df)
  # options(datadist='bb')
  #
  # f1 <- cph(Surv(time,OS)~ADAMTS12 + CHST11 + DCBLD2 + FN1 + FRMD6 + KRT17 + LOXL2 + MMP14 + NRP2 + PPFIBP1 + TGFBI + VCL,train.cox.df)
  # dca.res <- dca(f1, times = c(12,36,60))
  # ggplot(dca.res)

}




#' Net benifit
#'
#' @param response Should be 0 and 1
#' @param pred
#' @param xstart the starting point of the threshold probability, the default value is 0.01.
#' @param xstop the end point of the threshold probability, the default value is 0.99
#' @param step a numerical value specifying the incremental step of the threshold probability, the default value is 0.01
#' @param type controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”
#' @param model.name
#' @param exterdt specify whether the prediction should be made in external dataset; this is useful for cross validation of the model. By default, the value is NULL, indicating the prediction is made in the same dataset as the training dataset
#'
#' @return
#'
#' @examples
#' data(LIRI)
#'
#' pred <- as.vector( unlist( LIRI[,c(3)] ) )
#' ntbft(LIRI$status, pred)
#'
ntbft <- function(response, pred, xstart=0.01, xstop=0.99, step=0.01, type="treated", model.name = "Model", exterdt = NULL){
  # https://atm.amegroups.com/article/view/20389/pdf
  # For decision curve analysis. not exported

  if(length( setdiff(unique(response), c(0,1))) != 0 ){
    stop("Response or outcome should be 0 and 1")
  }

  pt <- seq(from=xstart,to=xstop,by=step)

  lpt <- length(pt)

  if(type=="treated") coef<-cbind(rep(1,lpt),rep(0,lpt))
  if(type=="untreated") coef<-cbind(rep(0,lpt),rep(1,lpt))
  if(type=="overall") coef<-cbind(rep(1,lpt),rep(1,lpt))
  if(type=="adapt") coef<-cbind(1-pt,pt)

  if(!is.null(exterdt)) response <- exterdt

  event.rate <- mean(response)

  nball  <- event.rate - (1-event.rate) * pt/(1-pt)
  nbnone <- 1 - event.rate - event.rate * (1-pt)/pt

  # pred and response should be of the same length
  N<-length(pred)
  nbt<-rep(NA,lpt)
  nbu<-rep(NA,lpt)
  for(t in 1:lpt){
    tp<-sum(pred>=pt[t] & response==1)
    fp<-sum(pred>=pt[t] & response==0)
    fn<-sum(pred<pt[t] & response==1)
    tn<-sum(pred<pt[t] & response==0)
    nbt[t]<-tp/N-fp/N*(pt[t]/(1-pt[t]))
    nbu[t]<-tn/N-fn/N*((1-pt[t])/pt[t])
  }
  nb<-data.frame(pt)
  names(nb)<-"threshold"
  nb["All"]<-coef[,1]*nball
  nb["None"]<-coef[,2]*nbnone
  nb[model.name]<-coef[,1]*nbt+coef[,2]*nbu
  return(nb)
}


#' Net benefit Bootstrap method to correct overfitting
#'
#' @param response Should be 0 and 1
#' @param pred
#' @param xstart the starting point of the threshold probability, the default value is 0.01.
#' @param xstop the end point of the threshold probability, the default value is 0.99
#' @param step a numerical value specifying the incremental step of the threshold probability, the default value is 0.01
#' @param type controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”
#' @param boots 500
#' @param model.name
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#'
#' pred <- as.vector( unlist( LIRI[,c(3)] ) )
#' ntbft.boot(LIRI$status, pred)
#'
ntbft.boot <- function(outcome, pred, boots = 500, xstart=0.01, xstop=0.99, step=0.01, type="treated", model.name = "Model"){

  diffnet <- function(pred, ii, outcome){
    # ii to allow boot() function to select a sample
    #cat(ii)
    nb <- ntbft(outcome[ii], pred[ii],
                xstart=xstart,xstop=xstop,
                step=step,type=type, model.name = model.name)

    nb0 <- ntbft(outcome, pred,
                 exterdt = outcome, xstart=xstart,
                 xstop=xstop, step=step, type=type, model.name = model.name)
    diff <- as.vector( unlist(nb[model.name] - nb0[model.name] ) )
    cat(".")
    return(diff)
  }

  library(boot)
  set.seed(124)
  rstls<- boot(pred, statistic=diffnet, R=boots, outcome=outcome)

  rstls
}


#' Plotting the Net benefit function. Coule be multiple curves
#'
#' @param nb Object from
#' @param nolines the number of the columns of nb which should be plotted using lines. The default is to plot all columns (except the first one containing the threshold).
#' @param nobands the number of the columns of nb which should be plotted using bands (useful to plot confidence intervals). The default is to plot no bands
#' @param ymin
#' @param ymax
#' @param legpos a vector of two coordinates indicating where the legend should be in the graph
#' @param palette
#'
#' @return
#'
#' @examples
#' p = loonR::decisionCurveAnalysisSimple(label, risk)
#' p$Plot
#' v.netbenefit = p$NetBenefit
#'
#' # similar to getr.netbenefit
#'
#' nb <- data.frame(threshold = v.netbenefit$threshold,
#'                  All = v.netbenefit$All,
#'                  None = v.netbenefit$None,
#'                  `Transcriptomic panel` = v.netbenefit$Model,
#'                  `Risk-stratification model` = r.netbenefit$Model,
#'                  stringsAsFactors = F,
#'                  check.names = F)
#'
#' p = loonR:::plot.ntbft(nb, 2:5)
#' p
#'
plot.ntbft <- function(nb, nolines = 2:dim(nb)[2], nobands = NULL, ymin = -0.1, ymax = max(nb[, c(nolines, nobands)], na.rm = T), legpos = c(0.9, 0.8), palette="aaas"){
  # https://atm.amegroups.com/article/view/20389/pdf
  # For decision curve analysis. not exported

  library(ggplot2)
  library(reshape2)
  library(ggpubr)


  ylow <- nb[, 1]
  yup <- nb[, 1]
  if (!is.null(nobands)) {
    ylow <- nb[, min(nobands)]
    yup <- nb[, max(nobands)]
  }
  nb.melt <- melt(nb[, c(1, nolines)],
                  id.vars = "threshold",
                  value.name = "Netbenefit", variable.name = "Models"
  )


  p <- ggplot(nb.melt) +
    geom_line(aes(
      x = threshold, y = Netbenefit,
      colour = Models, linetype = Models
    )) + scale_color_manual(values=loonR::get.ggsci.color(palette, n =length(unique(nb.melt$Models))  ) )
  p <- p + theme_pubr()

  # Plot confident interval
  if(!is.null(nobands)){
    p <- p +
      geom_ribbon(
        data = nb, aes(
          x = threshold,
          ymin = ylow, ymax = yup
        ),
        linetype = 2, alpha = 0.2
      )
  }
  p <- p + theme_pubr(legend = legpos) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    xlab("Threshold probability") +
    ylab("Net benefit")

  # modify the figure
  p <- p + rremove("legend.title")
  p

}



#' Decision curve analysis using risk probability and label
#'
#' @param label Should be 0 and 1
#' @param pred
#' @param xstart the starting point of the threshold probability, the default value is 0.01.
#' @param xstop the end point of the threshold probability, the default value is 0.99
#' @param step a numerical value specifying the incremental step of the threshold probability, the default value is 0.01
#' @param type controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”
#' @param model.name Default Model
#' @param boots 100 Bootstrape to avoid overfitting
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#'
#' pred <- as.vector( unlist( LIRI[,c(3)] ) )
#' decisionCurveAnalysisSimple(LIRI$status, pred)
#'
decisionCurveAnalysisSimple <- function(label, pred, xstart=0.01, xstop=0.99, step=0.01, type = "treated", model.name = "Model", boots = 0){

  # https://atm.amegroups.com/article/view/20389/pdf
  nb <- ntbft(label, pred, xstart = xstart, xstop = xstop, step = step, type = type, model.name = model.name)
  if(boots!=0){
    nb.boots <- ntbft.boot(label, pred, boots = boots, xstart = xstart, xstop = xstop, step = step, type = type, model.name = model.name)
    nb[model.name] = nb[model.name] - rowMeans(t(nb.boots$t))
  }

  p <- plot.ntbft(nb)
  print(p)

  res <- list(NetBenefit = nb, Plot = p)

}



#' Comparing two decision curves and reporting P value
#'
#' @param outcome
#' @param pred1
#' @param pred2
#' @param xstart the starting point of the threshold probability, the default value is 0.01.
#' @param xstop the end point of the threshold probability, the default value is 0.99
#' @param step a numerical value specifying the incremental step of the threshold probability, the default value is 0.01
#' @param type controls the type of net benefit to be computed. The allowed values correspond to the treated (“treated”), untreated (“untreated”) and overall (“overall”) patients, or to the ADAPT index (“adapt”). The default is the “treated”
#' @param model.name Default Model
#' @param boots 100 Bootstrape to avoid overfitting
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#'
#' pred1 <- as.vector( unlist( LIRI[,c(3)] ) )
#' pred2 <- as.vector( unlist( LIRI[,c(2)] ) )
#' decisionCurve.diff.pvalue(LIRI$status, pred1, pred2)
#'
decisionCurve.diff.pvalue <- function(outcome, pred1, pred2, boots = 100, xstart=0.01, xstop=0.99, step=0.01, type = "treated", model.name = "Model"){

  nbdiff <- function(outcome.infunction, ii, pred1, pred2){

    # ii to allow boot() function to select a sample
    #cat(ii)
    nb1 <- ntbft(outcome.infunction[ii], pred1[ii],
                 xstart=xstart, xstop=xstop,
                 step=step, type=type, model.name = model.name)

    nb2 <- ntbft(outcome.infunction[ii], pred2[ii],
                 xstart=xstart, xstop=xstop,
                 step=step, type=type, model.name = model.name)

    nb.diff <- as.vector( unlist(nb2[model.name] - nb1[model.name] ) )

    cat(".")
    return(nb.diff)
  }

  set.seed(127)
  library(boot)
  boot.diff <- boot(outcome, statistic = nbdiff,
                  R = boots, pred1 = pred1, pred2 = pred2)
  pvalue<-NULL
  for(i in 1:length(boot.diff$t0))
    pvalue<-c(pvalue,mean(abs(boot.diff$t[,i]-
                                boot.diff$t0[i])>abs(boot.diff$t0[i])))
  cat("\n","number of significant differences over threshold probabilities",
      xstart,"-",xstop,"=",sum(pvalue<=0.05),"\n")
  cat("\n","number of non-significant differences over threshold probabilities", xstart,"-",xstop,"=",sum(pvalue>0.05),"\n")

}





#' riskCalibrationPlot
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
riskCalibrationPlot<-function(group, pred, rms.method = FALSE, title = "Calibration plot", show.oberved.ci = FALSE,  bins = 10, color="npg",ticks.unit=0.25, full.range=TRUE, smooth.method = "loess", ...) {
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
#' @param ticks.unit 0.25 seq(0, 1, by = 0.25)
#' @param full.range Default TRUE. loess smoothing between 0-1 or first bin to last bin
#' @param smooth.method Smoothing method (function) to use, accepts either NULL or a character vector, e.g. "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, stats::lm, or stats::loess. "auto" is also accepted for backwards compatibility.
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
#' loonR::riskCalibrationPlot(factor(LIRI$status), d1$pred)
riskCalibrationPlot.default <- function(group, pred, rms.method = FALSE, title = "Calibration plot", show.oberved.ci = FALSE,  bins = 10, color="npg", ticks.unit=0.25, full.range=TRUE, smooth.method = "loess"){

  if(sum(pred > 1)!=0){
    stop("Pls set pred ranges from 0 to 1")
  }
  # Thanks reference: https://darrendahly.github.io/post/homr/
  # 公众号《绘制预测模型的校准曲线》用的riskRegression包（https://cran.r-project.org/web/packages/riskRegression/）也不错

  df <- data.frame(pred=pred, RawClass=group)
  if(!is.factor(group)){
    stop("group Must be a TRUE/FALSE factor")
    group = factor(group==unique(group)[2], levels = c(FALSE,TRUE))
  }
  df$Class=group


  if(rms.method){
    print( rms::val.prob(df$pred, as.numeric(df$Class))  )
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
    geom_abline(intercept = 0, slope = 1, linetype="dotted", color="gray")
  # loess fit through estimates，并且显示标准差
  if(full.range){
    # This plot the full range 0-1
    p1 = p1 + geom_smooth( aes(x = pred, y = as.numeric(Class)-1 ), level = 0.95,
                           color = "red", se = show.oberved.ci, method = smooth.method )
  }else{
    #this only plot between bins: the first bin_pred to the last bin_pred
    p1 = p1 + geom_smooth(aes(x = bin_pred, y = bin_prob ), level = 0.95,
                          color = "red", se = show.oberved.ci, method = smooth.method )
  }


  p1 = ggpar(p1, xlim = c(0,1), ylim = c(0,1), xticks.by =  ticks.unit,
             xlab = "", ylab="Observed proportion", title = title)


  # The distribution plot
    library(ggplot2)
  # p2 <- ggplot(df, aes(x=pred, color=RawClass, fill = RawClass)) +
  #     geom_histogram( position = "stack", color="#e9ecef", bins=bins, breaks = seq(0,1,1/bins) ) +
  #     scale_fill_manual(values=loonR::get.palette.color(color)) +
  #     theme_pubr()

    p2 <- ggdensity(df, x="pred", fill = "RawClass", bins = bins, rug = TRUE,
                      palette = color, position = "stack" )

  p2 =  p2 + xlab("Predicted probability") +  ylab("\nDensity")
  #p2 = p2 + scale_x_continuous(breaks = seq(0, 1, by=1/(bins)))
  p2 <- ggpar(p2, legend = "none", xlim = c(0,1), xticks.by = ticks.unit)


  plot_row <- cowplot::plot_grid(p1, p2, nrow = 2, ncol = 1, rel_heights = c(2,1))

  print(ResourceSelection::hoslem.test(as.numeric(df$Class)-1, df$pred, g = bins))

  #### rms plot
  # this function can also plot calbiration curve
  library(rms)

  rms::val.prob(  df$pred, as.integer(df$Class)-1,
                  logistic.cal= F, statloc = F, legendloc = F,
                  ylab = "Observed proportion"  )
  ####

  # The calibration performance was assessed with a predicted and observed mortalities plot and summarized across the full range of risk scores using the Hosmer-Lemeshow statistic.

  plot_row
}




#' Plot multiple calibration plot by riskRegression
#'
#' @param risk.list Data.frame with column risk or a named list
#' @param label a vector
#' @param palette
#'
#' @return
#' @export
#'
#' @examples
#'
#' data("LIRI")
#' risk.list = list(ANLN  = LIRI$ANLN,
#'                  CENPA = LIRI$CENPA)
#' label = LIRI$status
#'
#' risk.list = LIRI[,c(3,4)]
#'
riskCalibrationPlotMultiple <- function(risk.list = NULL, label = NULL, palette = "aaas"){

  if(!require(riskRegression)){
    devtools::install_github("tagteam/riskRegression")
    require(riskRegression)
  }

  if(is.null(risk.list)|is.null(label)){
    stop("Pls input risk and label list")
  }

  # if is a data frame
  if(is.data.frame(risk.list) | is.matrix(risk.list)){
    t = lapply(colnames(risk.list), function(x){ risk.list[[x]]})
    names(t) = colnames(risk.list)
    risk.list = t
    rm(t)
  }

  fit.res = lapply(names(risk.list), function(x){

    if(length(risk.list[[x]])!=length(label)){
      stop("Risk length should be the same as label length")
    }

    model.formula <- as.formula(paste('Label ~ ', x))

    df = data.frame(Risk  =  risk.list[[x]],
                    Label = label, check.names = F, stringsAsFactors = F)
    colnames(df)[1] = x

    fit <- glm(model.formula, data = df, family = binomial(link=logit) )
    fit
  })

  names(fit.res) = names(risk.list)


  ## merge list to data.frame
  risk.df = data.frame( do.call(cbind, risk.list), check.names = F)
  risk.df$Label = label

  xb <- Score(fit.res,
              formula = Label~1,
              null.model = FALSE,
              conf.int = TRUE,
              plots = c("calibration","ROC"),
              metrics = c("auc","brier"),
              B=1000, M=50,
              data = risk.df)  # risk df column names should be included all the variables and label

  plotCalibration(xb, rug = F, show.frequencies	= T, auc.in.legend = T,
                  col = loonR::get.palette.color(palette)[1:length(fit.res)],
                  )


}




#' Calibration plot
#'
#' @param risk Should be risk probability
#' @param label
#' @param color
#' @param bins Default 10
#'
#' @return
#' @export
#'
#' @examples
#' data("LIRI")
#' risk = LIRI$ANLN
#' risk = loonR::logit2prob(scale(risk))
#' label = LIRI$status
#' loonR::riskCalibrationPlotSingleComplex(risk, label)
#'
riskCalibrationPlotSingleComplex <- function(risk, label, color = "#4682B4", bins = 10){
  # inherite from rms.val
  if(!require(CalibrationCurves)){
    remotes::install_github("BavoDC/CalibrationCurves")
    require(CalibrationCurves)
  }


  val.prob.ci.2(risk, label,
                col.smooth = color,
                col.ideal = "gray",
                smooth ="loess", CL.smooth="fill",
                # col.smooth = scales::alpha(color, 0.5),
                lty.log=9, lwd.log=1.5, lwd.ideal = 0.5,
                g = bins
                )

}



#' Plot clinical impact: detail number of patients
#'
#' @param risk.score
#' @param label
#' @param name Default "Study"
#' @param palette Default "aaas"
#'
#' @return
#' @export
#'
#' @examples
clincial_impact <- function(risk.score, label, name="Study", palette = "aaas"){

    if(!require(rmda)){
      BiocManager::install("rmda")
      require(rmad)
    }
    # https://mdbrown.github.io/rmda/
    df = data.frame(Risk = risk.score, Label =label, stringsAsFactors = F, check.names = F)
    # https://blog.csdn.net/fjsd155/article/details/88951676
    simple<- decision_curve(Label~Risk, data = df, family = binomial(link ='logit'),
                            thresholds= seq(0,1, by = 0.01),
                            confidence.intervals =0.95,study.design = 'case-control',
                            population.prevalence = 0.3)

    plot_decision_curve(list(Study=simple),curve.names= name,
                        cost.benefit.axis =FALSE,col = c('red','blue'),
                        confidence.intervals =FALSE,standardize = FALSE)

    plot_clinical_impact(simple, population.size = 1000,cost.benefit.axis = T,
                         n.cost.benefits= 8,col = loonR::get.palette.color(palette)[1:2],
                         confidence.intervals= T)


}



#' Hosmer-Lemeshow Goodness of Fit (GOF) Test
#'
#' @param pred pred score
#' @param label Label
#' @param bin Default 10
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(123)
# n <- 500
# x <- rnorm(n)
# y <- rbinom(n, 1, plogis(0.1 + 0.5*x))
# HosmerLemesGoodnessOfFitTest(x, y)
HosmerLemesGoodnessOfFitTest <- function(pred=NULL, label=NULL, bin=10){

  # Method 1: generalhoslem::logitgof(repdata$type, fitted(multinom.fit))
  # Method 2: ResourceSelection::hoslem.test(mod$y, fitted(mod), g=10)
  # Method 3: DescTools::HosmerLemeshowTest(fit = fitted(m),obs= y)

  if(!require(ResourceSelection, quietly = T)){
    BiocManager::install("ResourceSelection")
  }

  ResourceSelection::hoslem.test(label, pred, g=bin)

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
#' riskCalibrationPlot.lrm(m)
riskCalibrationPlot.lrm <- function(rms.model, cox=FALSE){
  plot( rms::calibrate(rms.model,
                       method=c("boot"),
                       B=100,
                       smoother="lowess" ),
        xlab = "Predicted probability"
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
















