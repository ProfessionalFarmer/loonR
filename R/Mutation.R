

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
#' @param signatures.ref Signature * 96 context, colnames should be the same as signatures.nature2013 dataset
#'
#' @return
#' @export
#'
#' @examples
extractSignaturePercent <- function(chr=NULL, Start_Position = NULL, End_Position = NULL,  ref_base = NULL, alt_base = NULL, Sample_Barcode = NULL, signatures.ref = NULL){

  if(!require(deconstructSigs)){
    BiocManager::install("deconstructSigs")
  }
  if( is.null(chr) | is.null(Start_Position) | is.null(End_Position) | is.null(ref_base) |
      is.null(alt_base) | is.null(Sample_Barcode) | is.null(signatures.ref) ){
    stop("Chr,  star, end, barcode, ref and alt, and signatures.ref should not be NULL")
  }

  mutation.df <- data.frame(chr = chr, Start_Position = Start_Position,
                            End_Position = End_Position,
                            Ref = ref, Alt = alt,
                            Tumor_Sample_Barcode = Sample_Barcode)


  # Convert to deconstructSigs input
  sigs.input <- mut.to.sigs.input(mut.ref = mutation.df,
                                  sample.id = "Tumor_Sample_Barcode",
                                  chr = "Chromosome",
                                  pos = "Start_Position",
                                  ref = "Ref",
                                  alt = "Alt")


  colnames(sig.ref) = colnames(signatures.nature2013)


#  # todo calculate all the samples' signature
#
#   sample_1 = whichSignatures(tumor.ref = sigs.input,
#                              signatures.ref = sig.ref,
#                              #signatures.ref = signatures.nature2013,
#                              sample.id = "FP1705100059DN01",
#                              contexts.needed = TRUE,
#                              tri.counts.method = 'genome')
#


}


