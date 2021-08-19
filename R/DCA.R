#' Decision curve analysis
#'
#' @param data.frame.list row is sample. List should inlcude names
#' @param label
#' @param rms Default FALSE
#' @param validation.df shold include all variable in data.frame.list
#' @param palette
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
#'
decisionCurveAnalysis <- function(data.frame.list=NULL, label = NULL, rms=FALSE, validation.df = NULL, palette="aaas"){
  # https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis
  # https://atm.amegroups.com/article/view/20389/pdf

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


#' Plotting the Net benefit function
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
#' data(LIRI)
#'
#' pred <- as.vector( unlist( LIRI[,c(3)] ) )
#' ntbft.boot(LIRI$status, pred)
#' plot.ntbft(ntbft.boot)
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


