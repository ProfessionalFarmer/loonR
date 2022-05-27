#' Survival analysis
#'
#' @param Event Required 0,1
#' @param Time Required
#' @param Group Required
#' @param group.prefix Default ""
#' @param ylab
#' @param title Default "Suvival analysis"
#' @param palette the color palette to be used.
#' @param conf.int logical value. If TRUE, plots confidence interval.
#' @param legend.position one of c("top", "bottom", "left", "right", "none"). Default is "none".
#' @param linetype line types. Allowed values includes i) "strata" for changing linetypes by strata (i.e. groups); ii) a numeric vector (e.g., c(1, 2)) or a character vector c("solid", "dashed").
#' @param calculate.pval If TRUE, just reture p value data.frame.
#' @param only.consider.group Groups to consider
#' @param not.consider.group Groups to exclude
#' @param risk.table Default TRUE. Show the strata Table
#' @param remove.na If remove NA samples
#' @param best.point If best cut point wanted, pls input values by Group variable
#' @param cut.quantile Love this paramter, hansome
#' @param cut.label
#' @param surv.median.line character vector for drawing a horizontal/vertical line at median survival. Allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal.
#' @param pval logical value, a numeric or a string or show "HR" or "PHR" or "HRCI" or "PHRCI"
#' @param returnSurv.fit If TRUE, reture the fit object. Easy to check median survival time
#' @param maximum.days The maximum days or months
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' loonR::survivaly_analysis(LIRI$status, LIRI$time, LIRI$ANLN > mean(LIRI$ANLN), legend.position="right", risk.table = F )
#' loonR::survivaly_analysis(LIRI$status, LIRI$time, LIRI$ANLN, best.point = T)
survivaly_analysis <- function(Event = NULL, Time = NULL, Group = NULL, group.prefix = NA, ylab = "Survival probability",
                               title = "", palette = "lancet", conf.int = FALSE, legend.position="none",
                               linetype = 1, calculate.pval = FALSE, remove.na = FALSE,
                               only.consider.group = NULL, not.consider.group = NULL, risk.table = TRUE,
                               best.point = F, cut.quantile = NULL, cut.label = NULL,
                               surv.median.line = "none", pval = TRUE, returnSurv.fit = FALSE,
                               maximum.days=NULL){

  if(!require(survminer) | !require(survival)){BiocManager::install(c("survminer","survival"))}
  library(magrittr)

  if(is.null(Event) | is.null(Time) | is.null(Group)){
    stop("Please set Event and Time")
  }

  if(!is.null(maximum.days)){
    if(!is.numeric(maximum.days)){stop("Pls input days or months for maximum.days")}
    # There not event in the maximum point.
    Event[Time >= maximum.days & Event] = 0
    Time[Time >= maximum.days] = maximum.days
  }

  surv.analysis.df <- data.frame(Event=Event, Time=as.numeric(Time), Group = Group, stringsAsFactors = F)


  if(sum(is.na(surv.analysis.df)) != 0){
    if(remove.na){
      surv.analysis.df = surv.analysis.df[rowSums(is.na(surv.analysis.df))==0,]
    }else{
      stop("Pls set remove NA: TRUE")
    }
  }


  # try to find best point
  if(best.point){

    best.cut = survminer::surv_cutpoint(
      surv.analysis.df,
      time = "Time",
      event = "Event",
      variables = c("Group"),
      minprop = 0.1,
      progressbar = TRUE
    )$cutpoint$cutpoint

    cat("Best cut point is ", best.cut, "\n")

    surv.analysis.df$Variable = surv.analysis.df$Group
    surv.analysis.df$Group = loonR::splitGroupByCutoff(
      values = surv.analysis.df$Variable,
      cut.point = best.cut, cut.label = c("Low","High")
    )$New.Label

  }

  # if user specify
  if(!is.null(cut.quantile) & !is.null(cut.label) ){
    if(best.point){
        stop("Pls unselect best point")
    }
    surv.analysis.df$Variable = surv.analysis.df$Group
    surv.analysis.df$Group = loonR::splitGroupByCutoff(
      values = surv.analysis.df$Variable,
      quantile.cutoff = cut.quantile, cut.label = cut.label
    )$New.Label
  }else if(is.null(cut.quantile) & is.null(cut.label) ){

  }else{
    stop("Error, pls input cut point and cut label together")
  }

  surv.analysis.df$Event <- as.numeric(surv.analysis.df$Event)
  surv.analysis.df$Time <- as.numeric(surv.analysis.df$Time)

  if(!is.null(only.consider.group)){
    surv.analysis.df %<>% filter(Group %in% c(only.consider.group))
    if(length(surv.analysis.df$Group) < 2){stop("Please make sure two or more groups will be analyzed")}
  }
  if(!is.null(not.consider.group)){
    surv.analysis.df %<>% filter(!Group %in% c(not.consider.group))
    if(length(surv.analysis.df$Group) < 2){stop("Please make sure two or more groups will be analyzed")}
  }



  if(!is.na(group.prefix)){
    surv.analysis.df$Group = paste(group.prefix, surv.analysis.df$Group, sep="")
  }

  surv.fit <- survminer::surv_fit(Surv(Time, Event) ~ Group, data = surv.analysis.df)

  names(surv.fit$strata) <- gsub("Group=", "", names(surv.fit$strata))

  if(calculate.pval){
    return(surv_pvalue(surv.fit))
  }

  if(TRUE){ # always TRUE, 20220309, calculate p and HR
    if(!require(gtsummary)){
      BiocManager::install("gtsummary")
      require(gtsummary)
    }
    coxph.fit = coxph(Surv(Time, Event) ~ Group, data = surv.analysis.df)

    gtsummary::tbl_regression(coxph.fit, exp = TRUE)
    hazard.ratio = round( exp(coef(coxph.fit)), 2)
    hazard.ratio.ci = round( confint(coxph.fit,level = 0.95), 2 )
    ptext = surv_pvalue(surv.fit)$pval.txt

    cat("\n", ptext)
  }
  if(pval=="HR"){
    ptext = paste("HR = ",hazard.ratio, sep ="")
  }else if(pval=="HRCI"){
    ptext = paste("HR = ",hazard.ratio, " (",hazard.ratio.ci[1],"-",hazard.ratio.ci[2],")", sep ="")
  }else if(pval=="PHRCI"){
    ptext = paste(ptext,"\nHR = ",hazard.ratio, " (",hazard.ratio.ci[1],"-",hazard.ratio.ci[2],")", sep ="")
  }else if(pval=="PHR"){
    ptext = paste(ptext,"\nHR = ",hazard.ratio, sep ="")
  }
  if(pval %in% c("HR","HRCI","PHR", "PHRCI")){
    cat("\n", ptext)
    pval = ptext
  }


  if(returnSurv.fit){
    return(surv.fit)
  }

  p = ggsurvplot(surv.fit,
             risk.table = risk.table, ylab = ylab,
             risk.table.y.text.col = TRUE,
             risk.table.height = 0.4, legend.title = "",
             pval = pval, conf.int = conf.int, risk.table.y.text = TRUE,
             tables.y.text = FALSE, legend = legend.position,
             palette = palette, title = title,
             surv.median.line = surv.median.line)
  rm(surv.analysis.df, surv.fit)

  p

}


#' Determine the Optimal Cutpoint for Continuous Variables
#'
#' @param values
#' @param event
#' @param time
#'
#' @return
#' @export
#'
#' @examples
#'
#' data(LIRI)
#' res = loonR::findSurvivalCutPoint(LIRI[,3],LIRI$status, LIRI$time)
#' res$pval
#' res$surv.plot
#' res$estimated.point
findSurvivalCutPoint <- function(values = NULL, event = NULL, time = NULL, plot.Surv = FALSE){
  if(is.null(values) | is.null(event) | is.null(time) ) {
    stop("Pls input value, event and time")
  }
  library(dplyr)
  if(sum( ! unique(event) %in% c(1,0,TRUE,FALSE)  ) != 0 ){
    stop("pls input 0, 1, TRUE, FALSE in event variable")
  }
  df = data.frame(Variable = values,
                  Event = event,
                  Time = as.numeric(time),
                  stringsAsFactors = F)

  df$Time = as.numeric(df$Time)
  df$Variable = as.numeric(df$Variable)

  if(is.na(df$Event)){
    warning("We will remove NA, pls take care")
  }


  res = list()
  res$survminer.res = survminer::surv_cutpoint(
    df,
    time = "Time",
    event = "Event",
    variables = c("Variable"),
    minprop = 0.1,
    progressbar = TRUE
  )
  res$rawData = df

  res$estimated.point = res$survminer.res$cutpoint$cutpoint

  if(plot.Surv){
    surv.plot = loonR::survivaly_analysis(
      Event = df$Event, Time = df$Time,
      Group = loonR::splitGroupByCutoff(
        values = df$Variable,
        cut.point = res$estimated.point, cut.label = c("L","H")
      )$New.Label, remove.na = T, best.point = FALSE, palette = c("#ED0000", "#00468B") # here make sure best.point is FALSE
    )
    res$surv.plot = surv.plot
  }

  res$rawData$group = loonR::splitGroupByCutoff(
    values = df$Variable,
    cut.point = res$estimated.point, cut.label = c("L","H")
  )$New.Label

  p.val = loonR::survivaly_analysis(
    Event = df$Event, Time = df$Time,
    Group = loonR::splitGroupByCutoff(
      values = df$Variable,
      cut.point = res$estimated.point, cut.label = c("L","H")
    )$New.Label, remove.na = T
  )$pval
  res$pval = p.val

  res

}

