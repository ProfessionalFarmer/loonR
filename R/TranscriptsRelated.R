#' Remove isforms with extremly long interval (isoform.End - isoform.Start)
#'
#' @param gff gff file path
#' @param interval.length the maximum length between isoform start and isoform end
#'
#' @return
#' @export
#'
#' @examples remove.verylong.isoform.gff("/data/home2/Zhongxu/work/3.filter.gtf")
remove.verylong.isoform.gff <- function(gff, interval.length = 1500000){

  library(GenomicRanges)
  library(dplyr)
  library(plyranges)

  if(class(gff)=="character"){
    gff <- read_gff(gff)
  }else{
    gff <- gff
  }

  gff.tmp <- gff %>% select(transcript_id, gene_id, gene_name)
  gff.tmp <- as.data.frame(gff.tmp)

  pos.min <- aggregate(start ~ transcript_id, data = gff.tmp, min)
  pos.max <- aggregate(end   ~ transcript_id, data = gff.tmp, max)

  # 这个是在基因组上的长度，不是外显子的长度
  transcript.genomic.length <- merge(pos.min,pos.max,by="transcript_id")
  transcript.genomic.length$length <- transcript.genomic.length$end - transcript.genomic.length$start
  candidate.remove <- transcript.genomic.length[transcript.genomic.length$length > interval.length,]$transcript_id

  gff <- gff %>% filter(! transcript_id %in% candidate.remove)
  rm(gff.tmp, pos.min, pos.max)

  gff
}


#' Retrieve all transcripts related to target genes from gtf/gff file
#'
#' @param gff.Path Annotation file path or gff GRanges object
#' @param gene.vector Target genes
#' @param max.transcript.length Transcript start - end, default 1500000
#'
#' @return plyranges object
#' @export
#'
#' @examples
getGeneRelatedTranscripts.Gff <- function(gff.Path, gene.vector, max.transcript.length=1500000){


  library(dplyr)
  library(GenomicRanges)
  library(plyranges)


  if(class(gff.Path)=="character"){
    gr <- read_gff(gff.Path)
  }else{
    gr <- gff.Path
  }


  gr <- gr %>% filter(gene_id %in% gene.vector )

  # 查看是否有特别长的isoform
  gr <- loonR::remove.verylong.isoform.gff(gr, interval.length = max.transcript.length)
  gr

}



#' Retrieve protein sequence by biomaRt
#'
#' @param IDs
#' @param type Default: ensembl_transcript_id, Available: ensembl_transcript_id, ensembl_transcript_id_version, refseq_mrna, Possible filters are given by the listFilters function. biomaRt::listFilters(mart)
#'
#' @return
#' @export
#'
#' @examples
getTranscriptProteinSequence <- function(IDs,  type="ensembl_transcript_id"){

  library("biomaRt")
  # listMarts()
  mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "asia.ensembl.org")
  # listFilters(mart)
  getSequence(id = IDs,
              type = type, # refseq_mrna
              seqType = "peptide",
              mart = mart)


}



#' Convert gff file to data frame
#'
#' @param gff.Path
#'
#' @return
#' @export
#'
#' @examples
loadGff <- function(gff.path){
  library(dplyr)
  library(GenomicRanges)
  library(plyranges)

  gr <- read_gff(gff.path)
  gr

}

#' Write to gff file
#'
#' @param gff.path Output file path.
#' @param x Granges object.
#'
#' @return
#' @export
#'
#' @examples
writeGff <- function(x, gff.path){
  library(plyranges)
  write_gff(x, gff.path, index = FALSE)
}


#' Obtain transcript length from gff file or GRanges object
#'
#' @param gff Path or GRanges object
#' @param returnTranscriptLength Real transcripted length
#' @param returnEndStartLength Not the real transcript length. Transcript end - transcript start
#'
#' @return
#' @export
#'
#' @examples
getTranscriptLength <- function(gff, returnTranscriptLength=FALSE, returnEndStartLength=FALSE){

  if(returnTranscriptLength & returnEndStartLength){
    stop("Ccan't both")
  }
  if(!returnTranscriptLength & !returnEndStartLength){
    warning("Set default by column")
    returnTranscriptLength = TRUE
  }

  library(GenomicRanges)
  library(dplyr)
  library(plyranges)

  if(class(gff)=="character"){
    gff <- read_gff(gff)
  }else{
    gff <- gff
  }

  gff.tmp <- gff %>% select(transcript_id, gene_id, gene_name)
  gff.tmp <- as.data.frame(gff.tmp)

  if(returnEndStartLength){
    pos.min <- aggregate(start ~ transcript_id, data = gff.tmp, min)
    pos.max <- aggregate(end   ~ transcript_id, data = gff.tmp, max)

    # 这个是在基因组上的长度，不是外显子的长度
    transcript.genomic.length <- merge(pos.min,pos.max,by="transcript_id")
    transcript.genomic.length$Length <- transcript.genomic.length$end - transcript.genomic.length$start
    result = transcript.genomic.length
  }else if(returnTranscriptLength){
    gff.tmp$Length = gff.tmp$end - gff.tmp$start

    transcripted.length = gff.tmp %>% group_by(transcript_id) %>% summarise("exon_count"=n(), Length=sum(Length))
    result = transcripted.length

  }
  result

}

#' Obtain sequence length from fasta file
#'
#' @param fa.file
#' @param cores Default 20. Speed up
#'
#' @return A data.frame object
#' @export
#'
#' @examples
#' getFastaSummary(fa.path)
getFastaSummary <- function(fa.file, cores=20){
 if(!require(seqinr)){
   install.packages(seqinr)
 }
 fa <- seqinr::read.fasta(fa.file)

 n.seq <- length(fa)

 library(foreach)
 cl <- parallel::makeCluster(cores) #not to overload your computer
 doParallel::registerDoParallel(cl)


 res <- foreach(sequence=fa, .combine = rbind) %dopar% {
    name = seqinr::getName(sequence)
    description = seqinr::getAnnot(sequence)
    length = seqinr::getLength(sequence)
    gc = seqinr::GC(sequence)
    composition <- seqinr::count(sequence, wordsize = 1)
    seq = seqinr::c2s(sequence)

    c(name, description, length, gc, composition, seq)
 }


 #stop cluster
 parallel::stopCluster(cl)

 colnames(res) <- c("Name", "Annotation", "Length", "GC", "A", "C", "G", "T", "Sequence")
 res <- data.frame(res, stringsAsFactors = F)
 res$Length <- as.numeric(res$Length)

 res
}






