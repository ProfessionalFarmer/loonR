

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
#'
#' @return
#' @export
#'
#' @examples
survivaly_analysis <- function(Event = NULL, Time = NULL, Group = NULL, group.prefix = NA, ylab = "Survival probability",
                               title = "", palette = "lancet", conf.int = FALSE, legend.position="none",
                               linetype = 1, calculate.pval = FALSE, only.consider.group = NULL, not.consider.group = NULL, risk.table = TRUE){

  if(!require(survminer) | !require(survival)){BiocManager::install(c("survminer","survival"))}
  library(magrittr)

  if(is.null(Event) | is.null(Time) | is.null(Group)){
    stop("Please set Event and Time")
  }

  surv.analysis.df <- data.frame(Event=Event, Time=Time, Group = Group, stringsAsFactors = F)

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
             tables.y.text = FALSE, legend = "none",
             palette = palette, title = title)
  rm(surv.analysis.df, surv.fit)

  p

}




