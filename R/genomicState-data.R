#' Genomic State for Hsapiens.UCSC.hg19.knownGene
#'
#' Pre-computed genomic state for Hsapiens UCSC hg19 knownGene annotation built 
#' using \link[derfinder]{makeGenomicState} for 
#' TxDb.Hsapiens.UCSC.hg19.knownGene version 2.14.0. The object has been subset 
#' for chr21 only.
#'
#'
#' @name genomicState
#' @docType data
#' @format  A GRangesList with two components.
#' \describe{
#' \item{fullGenome }{ classifies each region as either being exon, intron or 
#' intergenic.}
#' \item{codingGenome }{ classfies the regions as being promoter, exon, intro, 
#' 5UTR, 3UTR or intergenic.}
#' }
#' @keywords datasets
#' @seealso \link[derfinder]{makeGenomicState}
NULL 
