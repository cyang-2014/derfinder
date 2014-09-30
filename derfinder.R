
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
    derfinder = citation('derfinder'), 
    knitrBootstrap = citation('knitrBootstrap'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown'),
    brainspan = RefManageR::BibEntry(bibtype = 'Unpublished', key = 'brainspan', title = 'Atlas of the Developing Human Brain [Internet]. Funded by ARRA Awards 1RC2MH089921-01, 1RC2MH090047-01, and 1RC2MH089929-01.', author = 'BrainSpan', year = 2011, url = 'http://developinghumanbrain.org'),
    originalder = RefManageR::BibEntry(bibtype = 'article', key = 'originalder', title = 'Differential expression analysis of RNA-seq data at single-base resolution', author = 'Alyssa C. Frazee and Sarven Sabunciyan and Kasper D. Hansen and Rafael A. Irizzary and Jeffrey T. Leek', year = 2013, journal = 'Biostatistics'),
    R = citation(),
    IRanges = citation('IRanges'),
    devtools = citation('devtools'),
    testthat = citation('testthat'),
    GenomeInfoDb = citation('GenomeInfoDb'),
    GenomicRanges = citation('GenomicRanges'),
    ggplot2 = citation('ggplot2'),
    biovizBase = citation('biovizBase'),
    bumphunter = citation('bumphunter'),
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
    bumphunterPaper = RefManageR::BibEntry(bibtype = 'article', key = 'bumphunterPaper', title = 'Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies', author = 'Jaffe, Andrew E and Murakami, Peter and Lee, Hwajin and Leek, Jeffrey T and Fallin, M Daniele and Feinberg, Andrew P and Irizarry, Rafael A', year = 2012, journal = 'International Journal of Epidemiology')
)

write.bibtex(bibs,
    file = 'derfinderRef.bib')
bib <- read.bibtex('derfinderRef.bib')

## Assign short names
names(bib) <- names(bibs)


## ----'start', message=FALSE----------------------------------------------
## Load libraries
library('derfinder')


## ----'phenoData', bootstrap.show.code=FALSE, results = 'asis'------------
library('knitr')
## Construct pheno table
pheno <- data.frame(
    gender = c('F', 'M', 'M', 'M', 'F', 'F', 'F', 'M', 'F', 'M', 'M', 'F'),
    lab = c('HSB97.AMY', 'HSB92.AMY', 'HSB178.AMY', 'HSB159.AMY', 'HSB153.AMY', 'HSB113.AMY', 'HSB130.AMY', 'HSB136.AMY', 'HSB126.AMY', 'HSB145.AMY', 'HSB123.AMY', 'HSB135.AMY'),
    Age = c(-0.547619047619048, -0.452380952380952, -0.571428571428571, -0.380952380952381, -0.666666666666667, -0.666666666666667, 21, 23, 30, 36, 37, 40),
    RIN = c(9.1, 9.2, 9.8, 9.9, 9.3, 9.4, 8.6, 8.1, 8.4, 7.4, 7.5, 8.5)
)
pheno$structure_acronym <- 'AMY'
pheno$structure_name <- 'amygdaloid complex'
pheno$file <- paste0('http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/bigwig/', pheno$lab, '.bw')
pheno$group <- factor(ifelse(pheno$Age < 0, 'fetal', 'adult'), levels = c('fetal', 'adult'))

## Display the main information
p <- pheno[, -which(colnames(pheno) %in% c('structure_acronym', 'structure_name', 'file'))]
kable(p, format = 'html', row.names = TRUE)


## ----'getData'-----------------------------------------------------------
## Determine the files to use and fix the names
files <- pheno$file
names(files) <- gsub('.AMY', '', pheno$lab)

## Load the data
system.time(fullCov <- fullCoverage(files = files, chrs = 'chr21'))


## ----'exploreFullCov'----------------------------------------------------
## Lets explore it
fullCov


## ----'filterCov'---------------------------------------------------------
## Filter coverage
filteredCov <- lapply(fullCov, filterData, cutoff = 2)


## ----'exploreFilteredCov'------------------------------------------------
## Similar to fullCov but with $position
filteredCov


## ----'compareCov'--------------------------------------------------------
## Compare the size in Mb
round(c(fullCov = object.size(fullCov), filteredCov = object.size(filteredCov)) / 1024^2, 1)


## ----'libSize'-----------------------------------------------------------
## Get some idea of the library sizes
sampleDepths <- sampleDepth(collapseFullCoverage(fullCov), 1)
sampleDepths


## ----'makeModels'--------------------------------------------------------
## Define models
models <- makeModels(sampleDepths, testvars = pheno$group, adjustvars = pheno[, c('gender', 'RIN')]) 

## Explore the models
lapply(models, head)


## ----'analyze'-----------------------------------------------------------
## Create a analysis directory
dir.create('analysisResults')
originalWd <- getwd()
setwd(file.path(originalWd, 'analysisResults'))

## Perform differential expression analysis
system.time(results <- analyzeChr(chr = 'chr21', filteredCov$chr21, models, groupInfo = pheno$group, writeOutput = TRUE, cutoffFstat = 1e-02, nPermute = 20, seeds = 20140923 + seq_len(20), returnOutput = TRUE))


## ----'exploreResults'----------------------------------------------------
## Explore
names(results)


## ----'exploreOptionsStats'-----------------------------------------------
## Explore optionsStats
names(results$optionsStats)

## Call used
results$optionsStats$analyzeCall


## ----'exploreCovPrep'----------------------------------------------------
## Explore coveragePrep
names(results$coveragePrep)

## Group means
results$coveragePrep$groupMeans


## ----'exploreFstats'-----------------------------------------------------
## Explore optionsStats
results$fstats

## Note that the length matches the number of bases used
identical(length(results$fstats), sum(results$coveragePrep$position))


## ----'exploreRegs'-------------------------------------------------------
## Explore regions
names(results$regions)


## ----'exploreRegs2'------------------------------------------------------
## Permutation summary information
results$regions[2:4]


## ----'exploreRegs3'------------------------------------------------------
## Candidate DERs
results$regions$regions


## ----'sensitivity'-------------------------------------------------------
## Width of potential DERs
summary(width(results$regions$regions))
sum(width(results$regions$regions) > 50)

## Width of candidate DERs
sig <- as.logical(results$regions$regions$significant)
summary(width(results$regions$regions[ sig ]))
sum(width(results$regions$regions[ sig ]) > 50)


## ----'exploreAnnotation'-------------------------------------------------
## Nearest annotation
head(results$annotation)


## ----'exploreTime'-------------------------------------------------------
## Time spent
results$timeinfo

## Use this information to make a plot
timed <- diff(results$timeinfo)
timed.df <- data.frame(Seconds = as.numeric(timed), Step = factor(names(timed),
    levels = rev(names(timed))))
library('ggplot2')
ggplot(timed.df, aes(y = Step, x = Seconds)) + geom_point()


## ----'mergeResults'------------------------------------------------------
## Go back to the original directory
setwd(originalWd)

## Merge results from several chromosomes. In this case we only have one.
mergeResults(chrs='chr21', prefix="analysisResults",
    genomicState = genomicState$fullGenome, 
    optionsStats = results$optionsStats)

## Files created by mergeResults()
dir('analysisResults', pattern = '.Rdata')


## ----'optionsMerge'------------------------------------------------------
## Options used to merge
load(file.path('analysisResults', 'optionsMerge.Rdata'))

## Contents
names(optionsMerge)

## Merge call
optionsMerge$mergeCall


## ----'exploreFullRegs'---------------------------------------------------
## Load all the regions
load(file.path('analysisResults', 'fullRegions.Rdata'))

## Metadata columns
names(mcols(fullRegions))


## ----'exploreFullAnnoRegs'-----------------------------------------------
## Load annotateRegions() output
load(file.path('analysisResults', 'fullAnnotatedRegions.Rdata'))

## Information stored
names(fullAnnotatedRegions)

## Take a peak
lapply(fullAnnotatedRegions, head)


## ----'extra'-------------------------------------------------------------
## Find overlaps between regions and summarized genomic annotation
annoRegs <- annotateRegions(fullRegions, genomicState$fullGenome)

## Indeed, the result is the same because we only used chr21
identical(annoRegs, fullAnnotatedRegions)

## Get the region coverage
regionCov <- getRegionCoverage(fullCov, fullRegions)

## Explore the result
head(regionCov[[1]])


## ----'derfinderPlot', eval = FALSE---------------------------------------
## library('derfinderPlot')
## 
## ## Overview of the candidate DERs in the genome
## plotOverview(regions = fullRegions, annotation = results$annotation, type = 'fwer')
## 
## suppressPackageStartupMessages(library('TxDb.Hsapiens.UCSC.hg19.knownGene'))
## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## 
## ## Base-levle coverage plots for the first 10 regions
## plotRegionCoverage(regions = fullRegions, regionCoverage = regionCov,
##     groupInfo = pheno$group, nearestAnnotation = results$annotation,
##     annotatedRegions = annoRegs, whichRegions=1:10, txdb = txdb, scalefac = 1,
##     ask = FALSE)
## 
## ## Cluster plot for the first region
## plotCluster(idx = 1, regions = fullRegions, annotation = results$annotation, coverageInfo = fullCov$chr21, txdb = txdb, groupInfo = pheno$group, titleUse = 'fwer')


## ----'regionMatrix'------------------------------------------------------
## Use regionMatrix()
system.time(regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76, totalMapped = 2^(sampleDepths), targetSize = 40e6))

## Explore results
names(regionMat$chr21)


## ----'exploreRegMatRegs'-------------------------------------------------
## regions output
regionMat$chr21$regions

## Number of regions
length(regionMat$chr21$regions)


## ----'exploreRegMatBP'---------------------------------------------------
## Base-level coverage matrices for each of the regions
## Useful for plotting
lapply(regionMat$chr21$bpCoverage[1:2], head, n = 2)

## Check dimensions. First region is 123 long, second one is 2 bp long.
## The columns match the number of samples (12 in this case).
lapply(regionMat$chr21$bpCoverage[1:2], dim)


## ----'exploreRegMatrix'--------------------------------------------------
## Dimensions of the coverage matrix
dim(regionMat$chr21$coverageMatrix)

## Coverage for each region. This matrix can then be used with limma or other pkgs
head(regionMat$chr21$coverageMatrix)


## ----'featureLevel'------------------------------------------------------
## Get the exon level matrix
system.time(exonCov <- coverageToExon(fullCov, genomicState$fullGenome, L = 76))

## Dimensions of the matrix
dim(exonCov)

## Explore a little bit
tail(exonCov)


## ----'regionMatAnnotate'-------------------------------------------------
## Annotate regions as exonic, intronic or intragenic
system.time(annoGenome <- annotateRegions(regionMat$chr21$regions, genomicState$fullGenome))
## Note that the genomicState object included in derfinder only has information
## for chr21 (hg19).

## Identify closest genes to regions
suppressPackageStartupMessages(library('bumphunter'))
system.time(annoNear <- annotateNearest(regionMat$chr21$regions, subject = 'hg19'))


## ----'static-vis', eval = FALSE------------------------------------------
## ## Identify the top regions by highest total coverage
## top <- order(regionMat$chr21$regions$area, decreasing = TRUE)[1:100]
## 
## ## Base-level plots for the top 100 regions with transcript information
## suppressPackageStartupMessages(library('TxDb.Hsapiens.UCSC.hg19.knownGene'))
## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## 
## library('derfinderPlot')
## plotRegionCoverage(regionMat$chr21$regions, regionCoverage = regionMat$chr21$bpCoverage, groupInfo = pheno$group, nearestAnnotation = annoNear, annotatedRegions = annoGenome, whichRegions = top, scalefac = 1, txdb = txdb, ask = FALSE)


## ----'epivizr', eval = FALSE---------------------------------------------
## ## Load epivizr, it's available from Bioconductor
## library('epivizr')
## 
## ## Load data to your browser
## mgr <- startEpiviz()
## ders_dev <- mgr$addDevice(
##     fullRegions[as.logical(fullRegions$significantFWER) ], "Candidate DERs")
## ders_potential_dev <- mgr$addDevice(
##     fullRegions[!as.logical(fullRegions$significantFWER) ], "Potential DERs")
## regs_dev <- mgr$addDevice(regionMat$chr21$regions, "Region Matrix")
## 
## ## Go to a place you like in the genome
## mgr$navigate("chr21", start(regionMat$chr21$regions[top[1]]) - 100, end(regionMat$chr21$regions[top[1]]) + 100)
## 
## ## Stop the navigation
## mgr$stopServer()


## ----'exportBigWig'------------------------------------------------------
## Subset only the first sample
fullCovSmall <- lapply(fullCov, '[', 1)

## Export to BigWig
bw <- createBw(fullCovSmall)

## See the file. Note that the sample name is used to name the file.
dir(pattern = '.bw')

## Internally createBw() coerces each sample to a GRanges object before 
## exporting to a BigWig file. If more than one sample was exported, the
## GRangesList would have more elements.
bw


## ----'advancedArg'-------------------------------------------------------
## URLs to advanced arguemtns
sapply(c('analyzeChr', 'loadCoverage'), advancedArg, browse = FALSE)
## Set browse = TRUE if you want to open them in your browser


## ----'citation'----------------------------------------------------------
## Citation info
citation('derfinder')


## ----createVignette, eval=FALSE, bootstrap.show.code=FALSE---------------
## ## Create the vignette
## library('knitrBootstrap')
## 
## knitrBootstrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'
## if(knitrBootstrapFlag) {
##     ## CRAN version
##     library('knitrBootstrap')
##     system.time(knit_bootstrap('derfinder.Rmd', chooser=c('boot', 'code'), show_code = TRUE))
##     unlink('derfinder.md')
## } else {
##     ## GitHub version
##     library('rmarkdown')
##     system.time(render('derfinder.Rmd', 'knitrBootstrap::bootstrap_document'))
## }
## ## Note: if you prefer the knitr version use:
## # library('rmarkdown')
## # system.time(render('derfinder.Rmd', 'html_document'))
## ## Clean up
## unlink('analysisResults', recursive = TRUE)
## file.remove('HSB97.bw')
## file.remove('derfinderRef.bib')
## 
## ## Extract the R code
## library('knitr')
## knit('derfinder.Rmd', tangle = TRUE)


## ----reproducibility1, echo=FALSE, bootstrap.show.code=FALSE-------------
## Date the vignette was generated
Sys.time()


## ----reproducibility2, echo=FALSE, bootstrap.show.code=FALSE-------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)


## ----reproducibility3, echo=FALSE, bootstrap.show.code=FALSE, bootstrap.show.message=FALSE----
## Session info
library('devtools')
session_info()


## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE-----
## Print bibliography
bibliography()

