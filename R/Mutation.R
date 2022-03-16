
#' Identify signature from sample data
#'
#' @param chr
#' @param Start_Position
#' @param End_Position
#' @param hg
#' @param pca
#' @param nmf
#' @param ref_base
#' @param alt_base
#' @param Sample_Barcode
#' @param project_name
#' @param n_sigs the number of signature will be identified
#'
#' @return
#' @export
#'
#' @examples
findSomaticSignature <- function(chr=NULL, Start_Position = NULL, End_Position = NULL, hg=38, pca = FALSE, nmf = FALSE,
                                 ref_base = NULL, alt_base = NULL, Sample_Barcode = NULL, project_name = "Project", n_sigs = 5){

  # https://bioconductor.org/packages/release/bioc/vignettes/SomaticSignatures/inst/doc/SomaticSignatures-vignette.html

  if( is.null(chr) | is.null(Start_Position) | is.null(End_Position) | is.null(ref_base) | is.null(alt_base) | is.null(Sample_Barcode) ){
    stop("Chr,  star, end, barcode, ref and alt should not be NULL")
  }

  if(!nmf & !pca){
    warning("You didn't specify the method to decompose, use NMF as default")
    nmf = TRUE
  }

  if(nmf==TRUE){
    m.decompose = nmfDecomposition
  }else if(pca==TRUE){
    m.decompose = pcaDecomposition
  }


  if(!require(SomaticSignatures)){
    BiocManager::install("SomaticSignatures")
  }
  if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)){
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
  }
  if(!require(BSgenome.Hsapiens.UCSC.hg38)){
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
  }
  if(!require(SomaticCancerAlterations)){
    BiocManager::install("SomaticCancerAlterations")
  }


  if(hg==38){
    g = BSgenome.Hsapiens.UCSC.hg38
  }else if(hg==37 | hg==19){
    g = BSgenome.Hsapiens.1000genomes.hs37d5
  }else{
    stop("Please set hg: 38, 19 or 37")
  }

  sca_vr = VRanges(
    seqnames =  stringr::str_remove_all(chr, "chr"),
    ranges = IRanges(Start_Position, End_Position),
    ref = ref_base,
    alt = alt_base,
    sampleNames = Sample_Barcode,
    #seqinfo = seqinfo(sca_data),
    study = project_name)

  # Motifs: Extracting the Sequence Context of Somatic Variants
  sca_motifs = mutationContext(sca_vr, g)

  # To continue with the estimation of the somatic signatures, the matrix \(M\) of the form {motifs × studies} is constructed. The normalize argument specifies that frequencies rather than the actual counts are returned.
  sca_mm = motifMatrix(sca_motifs, group = "study", normalize = TRUE)

  # The observed occurrence of the motifs, also termed somatic spectrum, can be visualized across studies, which gives a first impression of the data. The distribution of the motifs clearly varies between the studies.
  mutattionSpectrum <- plotMutationSpectrum(sca_motifs, "study")

  # Decomposition: Inferring Somatic Signatures
  sigs = identifySignatures(sca_mm, n_sigs, m.decompose)


  # Assessment: Number of Signatures
  n_sigs = 2:10
  gof = assessNumberSignatures(sca_mm, n_sigs)
  signature.number <- plotNumberSignatures(gof)

  # Visualization: Exploration of Signatures and Samples
  sample.signature.plot <- plotSignatureMap(sigs) + ggtitle("Somatic Signatures - Heatmap")
  signature.plot <- plotSignatures(sigs) + ggtitle("Somatic Signatures - Barchart")


  res = list(raw.motifs = sca_motifs, normalzied.motifs = sca_mm,
             mutation.spectrum.plot = mutattionSpectrum, signatures = sigs,
             assess.sig.number.plot = signature.number, sample.signature.plot = sample.signature.plot,
             signature.plot = signature.plot)



}







#' Estimate signature percents
#'
#' @param chr
#' @param Start_Position
#' @param End_Position
#' @param ref_base
#' @param alt_base
#' @param Sample_Barcode
#' @param signatures.ref 注意列名Default is deconstructSigs::signatures.cosmic  Signature * 96 context, colnames should be the same as signatures.nature2013 dataset or signatures.cosmic
#' @param ref.genome Default is hg38
#' @param tri.counts.method  'default', 'genome' or 'exome' or 'exome2genome'. Default is does not result in further normalization. genome is : the input data frame is normalized by number of times each trinucleotide context is observed in the genome
#'
#' @return
#' @export
#'
#' @examples
extractSignaturePercent <- function(chr=NULL, Start_Position = NULL, End_Position = NULL,  ref_base = NULL, alt_base = NULL, Sample_Barcode = NULL, signatures.ref = NULL, ref.genome = 'hg38', tri.counts.method = 'genome'){
  # https://github.com/raerose01/deconstructSigs
  if(!require(deconstructSigs)){
    BiocManager::install("deconstructSigs")
  }
  library(deconstructSigs)
  if(is.null(signatures.ref)){
    warning("Use signatures.cosmic as signatures reference")
    signatures.ref = deconstructSigs::signatures.cosmic
  }
  if( is.null(chr) | is.null(Start_Position) | is.null(End_Position) | is.null(ref_base) |
      is.null(alt_base) | is.null(Sample_Barcode) ){
    stop("Chr,  star, end, barcode, ref and alt, and signatures.ref should not be NULL")
  }
  if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)){
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
  }
  if(!require(BSgenome.Hsapiens.UCSC.hg38)){
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
  }


  warning(ref.genome, " used as reference genome")

  if(ref.genome==38){
    ref.genome = BSgenome.Hsapiens.UCSC.hg38
  }else if(ref.genome==37 | ref.genome==19){
    ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5
  }else{
    stop("Please set ref.genome: 38, 19 or 37")
  }

  mutation.df <- data.frame(Chromosome = chr,
                            Start_Position = as.numeric(Start_Position),
                            End_Position = as.numeric(End_Position),
                            Ref = ref_base,
                            Alt = alt_base,
                            Tumor_Sample_Barcode = Sample_Barcode,
                            check.names = FALSE)


  # Convert to deconstructSigs input
  sigs.input <- mut.to.sigs.input(mut.ref = mutation.df,
                                  sample.id = "Tumor_Sample_Barcode",
                                  chr = "Chromosome",
                                  pos = "Start_Position",
                                  ref = "Ref",
                                  alt = "Alt",
                                  bsg = ref.genome
                                  )


  # colnames(sig.ref) = colnames(signatures.nature2013)

 # todo calculate all the samples' signature
 warning("Pls note: the input data frame is normalized by number of times each trinucleotide context is observed in the ", tri.counts.method,"\n")
 warning("Default will not perform normalization while genome, exome, exome2genome will do\n")

  samples.signature.res = lapply(unique(mutation.df$Tumor_Sample_Barcode), function(sample){

    sample_res = whichSignatures(tumor.ref = sigs.input,
                               signatures.ref = signatures.ref,
                               #signatures.ref = signatures.nature2013,
                               sample.id = sample,
                               contexts.needed = TRUE,
                               tri.counts.method = tri.counts.method)

    sample_res
  } )
  names(samples.signature.res) = unique(mutation.df$Tumor_Sample_Barcode)

  ###### weights
  library(foreach)
  weights = foreach(sample=unique(mutation.df$Tumor_Sample_Barcode), .combine = rbind) %do% {
    samples.signature.res[[sample]]$weights
  }
  ##### tumor
  tumor = foreach(sample=unique(mutation.df$Tumor_Sample_Barcode), .combine = rbind) %do% {
    samples.signature.res[[sample]]$tumor
  }

  ##### product
  product = foreach(sample=unique(mutation.df$Tumor_Sample_Barcode), .combine = rbind) %do% {
    samples.signature.res[[sample]]$product
  }

  ##### diff
  diff = foreach(sample=unique(mutation.df$Tumor_Sample_Barcode), .combine = rbind) %do% {
    samples.signature.res[[sample]]$diff
  }

  ##### unknown
  unknown = foreach(sample=unique(mutation.df$Tumor_Sample_Barcode), .combine = rbind) %do% {
    samples.signature.res[[sample]]$unknown
  }
  rownames(unknown) = unique(mutation.df$Tumor_Sample_Barcode)

  res = list(
    sigs.input = sigs.input,
    weights = weights,
    tumor = tumor,
    product = product,
    diff = diff,
    unknown = unknown
  )

  res

}


#' Convert MAF format to a binary matrix
#'
#' @param maf File path or data.frame
#'
#' @return
#' @export
#'
#' @examples
maf2binaryTable <- function(maf){

  uuid = uuid::UUIDgenerate(use.time = T, n = 1L)

    # if(!require("SMDIC")){
    #   install.packages("SMDIC")
    #   library("SMDIC")
    # }
    # if(is.vector(maf)){
    #    res = SMDIC::maf2matrix(maf, percent = 0.01, nonsynonymous = F) # 提取sysnonymous和nonsynonymous的突变
    # }else if(is.matrix(maf) | is.data.frame(maf)){
    #   write.table(maf, paste0(uuid,".tsv"), row.names = F, sep="\t", quote = F)
    #   res = SMDIC::maf2matrix(paste0(uuid,".tsv"), percent = 0.01, nonsynonymous = F)
    #   file.remove(paste0(uuid,".tsv"))
    # }

    # The following method is more flexible
    if(!require("gnomeR")){
      devtools::install_github("AxelitoMartin/gnomeR")
      library("gnomeR")
    }
    if(is.vector(maf)){
      maf = read.table(maf, header = T, row.names = F, sep = "\t", quote = "")
    }
    res <- gnomeR::binmat(
      maf = maf, mut.type = "SOMATIC", SNP.only = FALSE,
      include.silent = FALSE, specify.plat = FALSE, recode.aliases = FALSE
    )

  res

}



