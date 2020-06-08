#' Count each event type and draw pie chart
#'
#' This function is evry similar to SuppaEventPie. Here the pie is draw based on ggplot2.
#'
#' @param data A data.frame or list object.
#' @param color ggsci color palette
#' @param colid If provide data.frame, column id need to be set
#' @param alpha Alpha value in plot
#' @param title Pie title
#' @param border Border color
#'
#' @return ggplot2 object
#' @export
#'
#' @examples DrawMutationType(res$Func.refGene, title = "# of mutation")
#' or DrawMutationType(res, col = 6, title = "# of mutation")
DrawMutationType <- function(data, color = "ucscgb", colid = 2, alpha = 0.7 , title = "", border="white"){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))

  library(ggsci)
  f <- parse(text=paste("pal_", color, sep = "")  )
  myPalette <- eval(f)(alpha = alpha)(n.color)

  Prop <- unclass(table(data))

  lbls <- names(  unclass(table(data))  )
  pct <- round(Prop/sum(Prop)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # add % to labels
  lbls <- paste(lbls,paste(" (",Prop,")",sep=""),sep="") # add value

  if (title == ""){
    title <- paste(" Total ", sum(Prop), sep = "" )
  }else{
    title <- paste(title, " (Total ", sum(Prop),")",sep = "" )
  }

  data <- data.frame(Count = Prop, Lable = lbls, stringsAsFactors = FALSE)

  # Basic piechart
  ggplot(data, aes(x="", y=Count, fill=Lable)) +
    geom_bar(stat="identity", width=1, color=border) +
    coord_polar("y", start=0) +
    #theme(legend.position="none") +
    #geom_text(aes(y = Count, label = Lable), color = "white", size=6) +
    scale_fill_manual(values = myPalette) + cowplot::theme_cowplot(font_family = "Arial") +
    theme_void()
}
