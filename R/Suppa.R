

#' Count each event type and draw pie chart
#'
#' @param data A data.frame or list object.
#' @param color ggsci color palette
#' @param colid If provide data.frame, column id need to be set
#' @param alpha Alpha value in plot
#' @param title Pie title
#' @param border Border color
#'
#' @return
#' @export
#'
#' @examples
#' suppa.event.pie(ioe.events.df$Type, title = "# of events")
#' or suppa.event.pie(ioe.events.df, col = 2, title = "# of events")


suppa.event.pie <- function(data, color = "jco", colid = 2, alpha =1 , title = "", border="white"){

  if( inherits(data, "data.frame")  ){
    data <- unique(data)
    data <- as.vector(data[,colid]) # now data is a vector class
  }
  n.color <- length(unique(data))
  if(n.color >= 10 ){
    stop("Please check, too many colors (More than 9)")
  }


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

  # ggpubr::ggpie(data.frame(Value=Prop, Type =  names(  unclass(table(data))  ) ),"Value",label = lbls, fill = "Type", color = border, palette =myPalette, lab.pos = "out", lab.adjust = 100 )
  # draw
  # You can change the border of each area with the classical parameters:
  pie(Prop , labels =lbls , border=border, col=myPalette, main = title )

}






#' Combine dpis, psi and tpm file together
#' I modified the code from https://github.com/comprna/SUPPA/blob/master/scripts/Volcano_MA_plot.R
#'
#' @param sample.names A vector specify the column name of each psi in psivec file.
#' @param psi Psivec file. Generate by Suppa2 diffSplice
#' @param dpsi Dpsi file. Generate by Suppa2 diffSplice file
#' @param event.tpm Event TPM file. Should be specified in Suppa analysis by --save_tpm_events
#' @param pval.cutoff P-value cutoff
#' @param dpsi.cutoff ΔPSI cutoff
#' @param tpm.cutoff  Mean TPM cutoff
#'
#' @return A data.frame
#' @export
#'
#' @examples
#'
#' suppa.get.final.table(sample.names = smp.names,
#'          dpsi = "analysis/05suppa/diffSplice.events.dpsi.temp.0",
#'          psi  = "analysis/05suppa/diffSplice.events.psivec",
#'          event.tpm = "analysis/05suppa/diffSplice.events_avglogtpm.tab",
#'          dpsi.cutoff = dpsi.cutoff, pval.cutoff = pval.cutoff, tpm.cutoff = tpm.cutoff)
#'
suppa.get.final.table <- function(sample.names = "", psi = "", dpsi = "", event.tpm = "", pval.cutoff = 0.05, dpsi.cutoff = 0.3, tpm.cutoff = 0){

  if (sample.names=="" | psi=="" | dpsi=="" | event.tpm==""){
    stop("NA found in parameter")
  }

  #Load the dspi file (output from SUPPA diffSplice)
  dpsi <- read.table(file = dpsi, sep = "\t", header = TRUE, , row.names = 1, stringsAsFactors = FALSE)

  colnames(dpsi) <- c("dPSI","p_value")

  #Load the psi file (output from SUPPA diffSplice)
  psi_events <- read.table(file = psi, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  colnames(psi_events) <- sample.names

  #Load the tpm of the events (output from SUPPA diffSplice, activating flag --save_tpm_events) events_avglogtpm.tab
  event_TPM <- read.table(file=event.tpm,sep="\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(event_TPM) <- c("event","mean_TPM")


  #Merge dpsi and psi
  merge1 <- merge(dpsi,psi_events,by="row.names")
  merge2 <- merge(merge1,event_TPM,by.x="Row.names",by.y="event")

  rownames(merge2) <- merge2$Row.names

  final_table <- merge2
  rm(psi_events,event_TPM, dpsi, merge2, merge1)

  final_table <- final_table[,-1]
  #
  # final_table <- final_table[!is.nan(final_table$dPSI),]
  # correct p value, do not use
  final_table$cpval <- p.adjust(final_table$p_value, method = "BH")

  final_table$log10pval <- -log10(final_table$p_value)

  # inf -> 4
  if( sum(is.infinite(final_table$log10pval)) > 0 ){
    final_table[is.infinite(final_table$log10pval),]$p_value <- 0.0001
    final_table[is.infinite(final_table$log10pval),]$log10pval <- 4
  }


  # suppa performs log10 transformation. I recover the real average TPM and convert to log2 transformation
  # this is different from suppa tutorial
  final_table$mean_TPM <- 10^final_table$mean_TPM
  final_table$logRNAc <- log2(final_table$mean_TPM+1) # avoid minus expression

  final_table$sig <- "not sig"
  final_table[final_table$p_value < pval.cutoff &
                final_table$logRNAc > tpm.cutoff &
                ( final_table$dPSI  > dpsi.cutoff | final_table$dPSI < -dpsi.cutoff ) ,]$sig <- "sig"


  final_table <- cbind(final_table,rownames(final_table))
  colnames(final_table)[ncol(final_table)] <- "Name"

  final_table$Gene <- sapply( strsplit(as.character(final_table$Name), ";"), function(x) paste(x[1])  )
  # if isoform is event type in event analysis
  final_table$Isoform <- sapply( strsplit(as.character(final_table$Name), ";") , function(x) paste(x[2])  )

  return(final_table)

}

















