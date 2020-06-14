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
  pct <- round(Prop/sum(Prop)*100,1)
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
                        position = "stack", flip=FALSE){

  library(ggplot2)
  library(ggsci)

  n.color <- length(unique(data$variable))

  f <- parse(text=paste("pal_", color, sep = "")  )

  if(n.color >= 10 ){
    myPalette <- pal_d3("category20",alpha = 0.8)(20)
    cat("Please note and check, too many colors (More than 9). Palette is forced to pal_d3 category20")
  }else{
    myPalette <- eval(f)(alpha = alpha)(n.color)
  }

  p <- ggplot(data, aes(x=Sample, y=as.numeric(value), fill=variable, label=value)) +
    geom_bar(position=position, stat="identity") + # position 参数为fill的时候，所有bar 等高
    labs(title=title,  x=x, y = y, fill="Class") +
    scale_fill_manual(values = myPalette) + cowplot::theme_cowplot(font_family = "Arial")
    #geom_text(size = 3, position = position_stack(vjust = 0.5))+  # 控制显示数目

  if (flip){
    p <- p + coord_flip()
  }else{
    p <- p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }

  p

}


#' Draw target coverage plot
#' Clone from https://gist.github.com/stephenturner/9396409   Thanks
#' input file from: bedtools coverage -hist -b  samp.01.bam -a target_regions.bed | grep ^all > samp.01.bam.hist.all.txt  or bedtools genomecov -ibam ../merge.bam -g genome.file.txt | grep ^genome > genome.cov
#' bedtools -g and -sorted should be noticed
#' Ref: http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
#' @param dir bedtools coverage output file should be in the same directory. File should suffix with "hist.all.txt$"
#' @return
#' @export
#'
#' @examples
#'
DrawTargetCoveragePlot <- function(dir){

  # Get a list of the bedtools output files you'd like to read in
  print(files <- list.files(path=dir,pattern="hist.all.txt$"))

  # Optional, create short sample names from the filenames.
  # For example, in this experiment, my sample filenames might look like this:
  # prefixToTrash-01.pe.on.pos.dedup.realigned.recalibrated.bam
  # prefixToTrash-02.pe.on.pos.dedup.realigned.recalibrated.bam
  # prefixToTrash-03.pe.on.pos.dedup.realigned.recalibrated.bam
  # This regular expression leaves me with "samp01", "samp02", and "samp03" in the legend.
  print(labs <- paste("", gsub(".hist.all.txt", "", files, perl=TRUE), sep=""))

  # Add path to file
  files <- paste(dir,files,sep="/")


  # Create lists to hold coverage and cumulative coverage for each alignment,
  # and read the data into these lists.
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
  }

  # Pick some colors
  # Ugly:
  # cols <- 1:length(cov)
  # Prettier:
  # ?colorRampPalette
  # display.brewer.all()
  library(RColorBrewer)
  library(ggsci)
  if(length(cov)<7){
    cols <- brewer.pal(length(cov), "Dark2")
  }else{
    cols <- pal_d3("category20",alpha = 0.8)(20)
  }
  #cols <- colorRampPalette(brewer.pal(8, "Accent"))(length(cov))

  # Save the graph to a file
  #png("exome-coverage-plots.png", h=1000, w=1000, pointsize=20)

  # Create plot area, but do not plot anything. Add gridlines and axis labels.
  plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
  abline(v = 20, col = "gray60")
  abline(v = 50, col = "gray60")
  abline(v = 80, col = "gray60")
  abline(v = 100, col = "gray60")
  abline(h = 0.50, col = "gray60")
  abline(h = 0.90, col = "gray60")
  axis(1, at=c(20,50,80), labels=c(20,50,80))
  axis(2, at=c(0.90), labels=c(0.90))
  axis(2, at=c(0.50), labels=c(0.50))

  # Actually plot the data for each of the alignments (stored in the lists).
  # 2018-02-26 Jason
  # for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])
  for (i in 1:length(cov)) points(c(0,cov[[i]][1:400, 2]), c(1,cov_cumul[[i]][1:400]), type='l', lwd=3, col=cols[i])

  # Add a legend using the nice sample labeles rather than the full filenames.
  legend("topright", legend=labs, col=cols, lty=1, lwd=4)


}






