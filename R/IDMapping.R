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
#' @param clean If set TRUE, will only return mapped ID without duplicates
#'
#' @return
#' @export
#'
#' @examples
#' keytypes(org.Hs.eg.db)
#'  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"
#'  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
#'  [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"
#'  [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
#'  [25] "UNIGENE"      "UNIPROT"
id_mapping <- function(IDS, key = "ENSEMBL", column = c("SYMBOL"), clean = FALSE ){

  library("org.Hs.eg.db")

  res <- clusterProfiler::bitr(IDS, fromType = key, toType = column, OrgDb = "org.Hs.eg.db")
  res <- data.frame(Ref = IDS,
                    res[match(IDS, res[,1]),]
  )

  library(dplyr)
  if(clean){
    # remove NA
    res = res[!is.na(as.character(unlist(res[,3]))),]
    # remove duplicates
    res = res[!duplicated(unlist(res[,3])),]
  }
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

  # list attributes
  # biomaRt::listAttributes(mart)

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
      'gene_biotype',
      "mirbase_accession",
      "mirbase_trans_name"),
    uniqueRows = TRUE)

  annotLookup <- dplyr::full_join(transcriptLookup, geneLookup, by = "ensembl_transcript_id")

  annotLookup

}



#' Convert Ensembl data frame to gene symbol dataframe
#'
#' @param ensembl.df Row is ENSG gene, col is sample
#'
#' @return
#' @export
#'
#' @examples
convertEnsemblDF2SymbolDF <- function(ensembl.df){

  tmp.df = ensembl.df

  library(dplyr)

  idmapping.res = loonR::id_mapping(rownames(tmp.df), key = "ENSEMBL", column  = "SYMBOL")

  # filter NA
  idmapping.res = idmapping.res %>% filter(!is.na(ENSEMBL)) %>% unique()
  # filter duplicate
  idmapping.res = idmapping.res[!duplicated(idmapping.res$SYMBOL),]

  # select
  tmp.df = tmp.df[idmapping.res$ENSEMBL,]
  rownames(tmp.df) = idmapping.res$SYMBOL

  tmp.df

}


