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
#'
#' @return
#' @export
#'
#' @examples
#' data(LIRI)
#' loonR::survivaly_analysis(LIRI$status, LIRI$time, LIRI$ANLN > mean(LIRI$ANLN), legend.position="right", risk.table = F )
survivaly_analysis <- function(Event = NULL, Time = NULL, Group = NULL, group.prefix = NA, ylab = "Survival probability",
                               title = "", palette = "lancet", conf.int = FALSE, legend.position="none",
                               linetype = 1, calculate.pval = FALSE, remove.na = FALSE,
                               only.consider.group = NULL, not.consider.group = NULL, risk.table = TRUE){

  if(!require(survminer) | !require(survival)){BiocManager::install(c("survminer","survival"))}
  library(magrittr)

  if(is.null(Event) | is.null(Time) | is.null(Group)){
    stop("Please set Event and Time")
  }

  surv.analysis.df <- data.frame(Event=Event, Time=Time, Group = Group, stringsAsFactors = F)

  if(sum(is.na(surv.analysis.df)) != 0){
    if(remove.na){
      surv.analysis.df = surv.analysis.df[rowSums(is.na(surv.analysis.df))==0,]
    }else{
      stop("Pls set remove NA: TRUE")
    }
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

  p = ggsurvplot(surv.fit,
             risk.table = risk.table, ylab = ylab,
             risk.table.y.text.col = TRUE,
             risk.table.height = 0.4, legend.title = "",
             pval = TRUE, conf.int = conf.int, risk.table.y.text = TRUE,
             tables.y.text = FALSE, legend = legend.position,
             palette = palette, title = title)
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
findSurvivalCutPoint <- function(values = NULL, event = NULL, time = NULL){
  if(is.null(values) | is.null(event) | is.null(time) ) {
    stop("Pls input value, event and time")
  }
  library(dplyr)
  if(sum( ! unique(event) %in% c(1,0,TRUE,FALSE)  ) != 0 ){
    stop("pls input 0, 1, TRUE, FALSE in event variable")
  }
  df = data.frame(Variable = values,
                  Event = event,
                  Time = time, stringsAsFactors = F)

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

  res$estimated.point = res$survminer.res$Variable$estimate

  surv.plot = loonR::survivaly_analysis(
      Event = df$Event, Time = df$Time,
      Group = loonR::splitGroupByCutoff(
        values = df$Variable,
        cut.point = res$estimated.point, cut.label = c("L","H")
      )$New.Label, remove.na = T
  )

  res$surv.plot = surv.plot

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

