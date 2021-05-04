

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


#' Title
#'
#' @param ensembl.ids
#'
#' @return
#' @export
#'
#' @examples
ensembl_EntrezID <- function(ensembl.ids){
  library("org.Hs.eg.db")
  library(dplyr)
  df <- as.data.frame(org.Hs.egENSEMBL2EG) %>%  filter(gene_id %in% candidate.tar$`Target Gene (Entrez ID)`)
  df
}


#' Title
#'
#' @param IDS
#' @param key
#' @param column
#'
#' @return
#' @export
#'
#' @examples
#'
#' keytypes(org.Hs.eg.db)
#'  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"
#'  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
#'  [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"
#'  [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
#'  [25] "UNIGENE"      "UNIPROT"
#'
id_mapping <- function(IDS, key = "ENSEMBL", column = c("SYMBOL") ){

  library("org.Hs.eg.db")

  res <- clusterProfiler::bitr(IDS, fromType = key, toType = column, OrgDb = "org.Hs.eg.db")
  res <- data.frame(Ref = IDS,
                    res[match(IDS, res[,1]),]
  )

  res

}

#' Retrieve full ID mapping table
#'
#' @return A full annotation table
#' @export
#'
#' @examples get_full_mapping_table
get_full_mapping_table <- function(){

  # https://www.biostars.org/p/384296/
  require('biomaRt')

  mart <- useMart('ENSEMBL_MART_ENSEMBL')
  mart <- useDataset('hsapiens_gene_ensembl', mart)

  # Check that it is indeed GRCh38:
  searchDatasets(mart = mart, pattern = 'hsapiens')

  # Now generate the table:
  transcriptLookup <- getBM(
      mart = mart,
      attributes = c(
        'ensembl_transcript_id',
        'ensembl_gene_id',
        'transcript_mane_select',
        'transcript_biotype'),
      uniqueRows = TRUE)

  geneLookup <- getBM(
    mart = mart,
    attributes = c(
      'ensembl_transcript_id',
      'ensembl_gene_id',
      'hgnc_symbol',
      'hgnc_id',
      'external_gene_name',
      "entrezgene_id",
      'gene_biotype'),
    uniqueRows = TRUE)

  annotLookup <- dplyr::full_join(transcriptLookup, geneLookup, by = "ensembl_transcript_id")

  annotLookup

}


