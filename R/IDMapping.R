

#' Convert RefSeq ID to Ensembl ID by biomaRt
#'
#' RefSeq identifiers differ in format according to the type of record the identifiers are for as shown below:
#'
#' NG_XXXXX: RefSeq accessions for genomic region (nucleotide) records
#'
#' NM_XXXXX: RefSeq accessions for mRNA records
#'
#' NC_XXXXX: RefSeq accessions for chromosome records
#'
#' NP_XXXXX: RefSeq accessions for protein records
#'
#' XR_XXXXX: RefSeq accessions for model RNAs that are not associated with protein products
#'
#' XM_XXXXX: RefSeq accessions for model mRNA records
#'
#' XP_XXXXX: RefSeq accessions for model protein records
#'
#' Where XXXXX is a sequence of integers.
#'
#' NCBI https://www.ncbi.nlm.nih.gov/RefSeq/ allows users to query the RefSeq database using RefSeq identifiers.
#'
#'
#' @param nm.accessions
#'
#' @return
#' @export
#'
#' @examples
refSeq_Ensembl <- function(nm.accessions){
  library("biomaRt")
  #listMarts()
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  # listFilters(ensembl)
  # listAttributes(ensembl)
  getBM(  filters="refseq_mrna",
    attributes=c("refseq_mrna", "ensembl_transcript_id", "entrezgene_id", "hgnc_symbol"),
    values=nm.accessions, mart=ensembl      )
}


ensembl_EntrezID <- function(ensembl.ids){
  library("org.Hs.eg.db")
  library(dplyr)
  df <- as.data.frame(org.Hs.egENSEMBL2EG) %>%  filter(gene_id %in% candidate.tar$`Target Gene (Entrez ID)`)
  df
}



