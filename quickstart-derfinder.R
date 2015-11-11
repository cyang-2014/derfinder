## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-----------
## Track time spent on making the vignette
startTime <- Sys.time()

## Bib setup
library('knitcitations')

## Load knitcitations with a clean bibliography
cleanbib()
cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
# Note links won't show for now due to the following issue
# https://github.com/cboettig/knitcitations/issues/63

## Write bibliography information
bibs <- c(knitcitations = citation('knitcitations'),
    derfinder = citation('derfinder')[1], 
    BiocStyle = citation('BiocStyle'),
    knitrBootstrap = citation('knitrBootstrap'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown'),
    brainspan = RefManageR::BibEntry(bibtype = 'Unpublished', key = 'brainspan', title = 'Atlas of the Developing Human Brain [Internet]. Funded by ARRA Awards 1RC2MH089921-01, 1RC2MH090047-01, and 1RC2MH089929-01.', author = 'BrainSpan', year = 2011, url = 'http://developinghumanbrain.org'),
    originalder = citation('derfinder')[2],
    R = citation(),
    IRanges = citation('IRanges'),
    devtools = citation('devtools'),
    testthat = citation('testthat'),
    GenomeInfoDb = citation('GenomeInfoDb'),
    GenomicRanges = citation('GenomicRanges'),
    ggplot2 = citation('ggplot2'),
    biovizBase = citation('biovizBase'),
    bumphunter = citation('bumphunter')[1],
    TxDb.Hsapiens.UCSC.hg19.knownGene = citation('TxDb.Hsapiens.UCSC.hg19.knownGene'),
    AnnotationDbi = citation('AnnotationDbi'),
    BiocParallel = citation('BiocParallel'),
    derfinderHelper = citation('derfinderHelper'),
    GenomicAlignments = citation('GenomicAlignments'),
    GenomicFeatures = citation('GenomicFeatures'),
    GenomicFiles = citation('GenomicFiles'),
    Hmisc = citation('Hmisc'),
    qvalue = citation('qvalue'),
    Rsamtools = citation('Rsamtools'),
    rtracklayer = citation('rtracklayer'),
    S4Vectors = citation('S4Vectors'),
    bumphunterPaper = RefManageR::BibEntry(bibtype = 'article', key = 'bumphunterPaper', title = 'Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies', author = 'Jaffe, Andrew E and Murakami, Peter and Lee, Hwajin and Leek, Jeffrey T and Fallin, M Daniele and Feinberg, Andrew P and Irizarry, Rafael A', year = 2012, journal = 'International Journal of Epidemiology'),
    derfinderData = citation('derfinderData')
)

write.bibtex(bibs,
    file = 'quickstartRef.bib')
bib <- read.bibtex('quickstartRef.bib')

## Assign short names
names(bib) <- names(bibs)

## Working on Windows?
windowsFlag <- .Platform$OS.type == 'windows'

## ----'installDer', eval = FALSE------------------------------------------
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("derfinder")

## ----'citation'----------------------------------------------------------
## Citation info
citation('derfinder')

## ----'start', message=FALSE----------------------------------------------
## Load libraries
library('derfinder')
library('derfinderData')
library('GenomicRanges')

## ----'locateAMYfiles'----------------------------------------------------
## Determine the files to use and fix the names
files <- rawFiles(system.file('extdata', 'AMY', package = 'derfinderData'), samplepatt = 'bw', fileterm = NULL)
names(files) <- gsub('.bw', '', names(files))

## ----'getData', eval = !windowsFlag--------------------------------------
## Load the data from disk
fullCov <- fullCoverage(files = files, chrs = 'chr21', verbose = FALSE)

## ----'getDataWindows', eval = windowsFlag, echo = FALSE------------------
## Load data in Windows case
foo <- function() { 
    load(system.file('extdata', 'fullCov', 'fullCovAMY.RData', package = 'derfinderData'))
    return(fullCovAMY) 
}
fullCov <- foo()

## ----'regionMatrix'------------------------------------------------------
## Use regionMatrix()
regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76, verbose = FALSE)

## ----'exploreRegMatRegs'-------------------------------------------------
## regions output
regionMat$chr21$regions

## ----'exploreRegMatBP'---------------------------------------------------
## Base-level coverage matrices for each of the regions
## Useful for plotting
lapply(regionMat$chr21$bpCoverage[1:2], head, n = 2)

## Check dimensions. First region is 565 long, second one is 138 bp long.
## The columns match the number of samples (12 in this case).
lapply(regionMat$chr21$bpCoverage[1:2], dim)

## ----'exploreRegMatrix'--------------------------------------------------
## Dimensions of the coverage matrix
dim(regionMat$chr21$coverageMatrix)

## Coverage for each region. This matrix can then be used with limma or other pkgs
head(regionMat$chr21$coverageMatrix)

## ----'phenoData'---------------------------------------------------------
## Get pheno table
pheno <- subset(brainspanPheno, structure_acronym == 'AMY')

## ----'identifyDERsDESeq2'------------------------------------------------
## Required
library('DESeq2')

## Round matrix
counts <- round(regionMat$chr21$coverageMatrix)

## Round matrix and specify design
dse <- DESeqDataSetFromMatrix(counts, pheno, ~ group + gender)

## Perform DE analysis
dse <- DESeq(dse, test = 'LRT', reduced = ~ gender, fitType = 'local')

## Extract results
mcols(regionMat$chr21$regions) <- c(mcols(regionMat$chr21$regions), results(dse))

## Save info in an object with a shorter name
ers <- regionMat$chr21$regions
ers

## ----'vennRegions'-------------------------------------------------------
## Find overlaps between regions and summarized genomic annotation
annoRegs <- annotateRegions(ers, genomicState$fullGenome, verbose = FALSE)

library('derfinderPlot')
venn <- vennRegions(annoRegs, counts.col = 'blue')

## ----'nearestGene'-------------------------------------------------------
## Load database of interest
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
txdb <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene, 'chr21')

## Find nearest feature
library('bumphunter')
genes <- annotateTranscripts(txdb)
annotation <- matchGenes(ers, subject = genes)

## Extract the region coverage
regionCov <- regionMat$chr21$bpCoverage

## ----'overviewNearest'---------------------------------------------------
## Overview of the candidate DERs in the genome
plotOverview(regions = ers, annotation = annotation, type = 'annotation')

## ----'first5'------------------------------------------------------------
plotRegionCoverage(regions = ers, regionCoverage = regionCov, 
    groupInfo = pheno$group, nearestAnnotation = annotation, 
    annotatedRegions = annoRegs, whichRegions = seq_len(5), txdb = txdb, scalefac = 1, 
    ask = FALSE)

## ----createVignette, eval=FALSE------------------------------------------
## ## Create the vignette
## library('rmarkdown')
## system.time(render('quickstart-derfinder.Rmd', 'BiocStyle::html_document'))
## 
## ## Extract the R code
## library('knitr')
## knit('quickstart-derfinder.Rmd', tangle = TRUE)

## ----createVignette2-----------------------------------------------------
## Clean up
file.remove('quickstartRef.bib')

## ----reproducibility1, echo=FALSE----------------------------------------
## Date the vignette was generated
Sys.time()

## ----reproducibility2, echo=FALSE----------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)

## ----reproducibility3, echo=FALSE----------------------------------------
## Session info
library('devtools')
options(width = 120)
session_info()

## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE-----
## Print bibliography
bibliography()

