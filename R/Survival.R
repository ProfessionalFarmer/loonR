

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
#'
#' @return
#' @export
#'
#' @examples
survivaly_analysis <- function(Event = NULL, Time = NULL, Group = NULL, group.prefix = NA, ylab = "Time",
                               title = "Suvival analysis", palette = "lancet", conf.int = FALSE, legend.position="none",
                               linetype = 1, calculate.pval = FALSE){

  if(!require(survminer) | !require(survival)){BiocManager::install(c("survminer","survival"))}

  if(is.null(Event) | is.null(Time) | is.null(Group)){
    stop("Please set Event and Time")
  }

  surv.analysis.df <- data.frame(Event=Event, Time=Time, stringsAsFactors = F)

  surv.analysis.df$Event <- as.numeric(surv.analysis.df$Event)
  surv.analysis.df$Group <- Group
  surv.analysis.df$Time <- as.numeric(surv.analysis.df$Time)

  if(!is.na(group.prefix)){
    surv.analysis.df$Group = paste(group.prefix, surv.analysis.df$Group, sep="")
  }

  surv.fit <- survminer::surv_fit(Surv(Time, Event) ~ Group, data = surv.analysis.df)

  names(surv.fit$strata) <- gsub("Group=", "", names(surv.fit$strata))

  if(calculate.pval){
    return(surv_pvalue(surv.fit))
  }

  p = ggsurvplot(surv.fit,
             risk.table = TRUE, ylab = ylab,
             risk.table.y.text.col = TRUE,
             risk.table.height = 0.4,
             pval = TRUE, conf.int = conf.int, risk.table.y.text = TRUE,
             tables.y.text = FALSE, legend = "none",
             palette = palette)
  rm(surv.analysis.df, surv.fit)

  p

}




