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
#' @param gff.Path Annotation file path
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

  gr <- read_gff(gff.Path)
  gr <- gr %>% filter(gene_id %in% gene.vector )

  # 查看是否有特别长的isoform
  gr <- loonR::remove.verylong.isoform.gff(gr, interval.length = max.transcript.length)
  gr

}



#' Retrieve protein sequence by biomaRt
#'
#' @param IDs
#' @param type Default: ens_hs_transcript. Available: Possible filters are given by the listFilters function. biomaRt::listFilters(ensembl)
#'
#' @return
#' @export
#'
#' @examples
getTranscriptProteinSequence <- function(IDs,  type="ens_hs_transcript"){

  library("biomaRt")
  # listMarts()
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  # listFilters(ensembl)
  getSequence(id = IDs,
              type = type, # refseq_mrna
              seqType = "peptide",
              mart = ensembl)


}






