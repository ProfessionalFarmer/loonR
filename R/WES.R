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
DrawMutationTypePie <- function(data, color = "jco", colid = 2, alpha = 0.8 , title = "", border="white"){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))

  library(ggsci)
  library(ggplot2)
  f <- parse(text=paste("pal_", color, sep = "")  )

  if(n.color >= 10 ){
    myPalette <- pal_d3("category20",alpha = 0.8)(20)
    cat("Please note and check, too many colors (More than 9). Palette is forced to pal_d3 category20")
  }else{
    myPalette <- eval(f)(alpha = alpha)(n.color)
  }

  Prop <- unclass(table(data))

  lbls <- names(  unclass(table(data))  )
  pct <- round(Prop/sum(Prop)*100,4)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # add % to labels
  lbls <- paste(lbls,paste(" (",Prop,")",sep=""),sep="") # add value

  if (title == ""){
    title <- paste(" Total ", sum(Prop), sep = "" )
  }else{
    title <- paste(title, " (Total ", sum(Prop),")",sep = "" )
  }

  data <- data.frame(Count = Prop, Type = lbls, stringsAsFactors = FALSE)

  # Basic piechart
  ggplot(data, aes(x="", y=Count, fill=Type)) +
    geom_bar(stat="identity", width=1, color=border) +
    coord_polar("y", start=0) +
    labs(title=title,  x ="", y = "") +
    #theme(legend.position="none") +
    #geom_text(aes(y = Count, label = Type), color = "white", size=6) +
    scale_fill_manual(values = myPalette) + cowplot::theme_cowplot(font_family = "Arial") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}



#' Draw bar plot
#'
#' @param data A data.frame. Column names for ggplot2 :x=Sample, y=as.numeric(value), fill=variable
#' @param color ggsci color palette
#' @param alpha Color alpha
#' @param title Titile of plot
#' @param x X label
#' @param y Y label
#' @param position
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' data <- data.frame( variable = c("exonic","exonic","exonic","intronic","intronic","intronic"),
#'                     Sample = c("S1","S2","S3","S1","S2","S3"),
#'                     value = c(1,2,3,4,5,6)
#' )
#' DrawBarPlot(data)
#'
DrawMutationTypeBarPlot <- function(data, color = "jco", alpha = 0.8 ,
                        title = "", x = "Sample", y = "Number of mutations (SNV+INDEL)",
                        position = "stack"){

  library(ggplot2)
  library(ggsci)

  f <- parse(text=paste("pal_", color, sep = "")  )

  if(n.color >= 10 ){
    myPalette <- pal_d3("category20",alpha = 0.8)(20)
    cat("Please note and check, too many colors (More than 9). Palette is forced to pal_d3 category20")
  }else{
    myPalette <- eval(f)(alpha = alpha)(n.color)
  }

  ggplot(data, aes(x=Sample, y=as.numeric(value), fill=variable, label=value)) +
    geom_bar(position=position, stat="identity") + # position 参数为fill的时候，所有bar 等高
    labs(title=title,  x=x, y = y, fill="Class") +
    scale_fill_manual(values = myPalette) + cowplot::theme_cowplot(font_family = "Arial") +
    #geom_text(size = 3, position = position_stack(vjust = 0.5))+  # 控制显示数目
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

}





