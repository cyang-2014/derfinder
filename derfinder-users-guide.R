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
    file = 'derfinderUsersGuideRed.bib')
bib <- read.bibtex('derfinderUsersGuideRed.bib')

## Assign short names
names(bib) <- names(bibs)

## Working on Windows?
windowsFlag <- .Platform$OS.type == 'windows'

## ----'start', message=FALSE----------------------------------------------
## Load libraries
library('derfinder')
library('derfinderData')
library('GenomicRanges')

## ----'phenoData', results = 'asis'---------------------------------------
library('knitr')
## Get pheno table
pheno <- subset(brainspanPheno, structure_acronym == 'AMY')

## Display the main information
p <- pheno[, -which(colnames(pheno) %in% c('structure_acronym', 'structure_name', 'file'))]
rownames(p) <- NULL
kable(p, format = 'html', row.names = TRUE)

## ----'getData', eval = !windowsFlag--------------------------------------
## Determine the files to use and fix the names
files <- rawFiles(system.file('extdata', 'AMY', package = 'derfinderData'), samplepatt = 'bw', fileterm = NULL)
names(files) <- gsub('.bw', '', names(files))

## Load the data from disk
system.time(fullCov <- fullCoverage(files = files, chrs = 'chr21'))

## ----'getDataWindows', eval = windowsFlag, echo = FALSE------------------
## Load data in Windows case
foo <- function() { 
    load(system.file('extdata', 'fullCov', 'fullCovAMY.RData', package = 'derfinderData'))
    return(fullCovAMY) 
}
fullCov <- foo()

## ----'webData', eval = FALSE---------------------------------------------
## ## Determine the files to use and fix the names
## files <- pheno$file
## names(files) <- gsub('.AMY', '', pheno$lab)
## 
## ## Load the data from the web
## system.time(fullCov <- fullCoverage(files = files, chrs = 'chr21'))

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

## ----'regionMatrix'------------------------------------------------------
## Use regionMatrix()
system.time(regionMat <- regionMatrix(fullCov, cutoff = 30, L = 76))

## Explore results
class(regionMat)
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

## Check dimensions. First region is 565 long, second one is 138 bp long.
## The columns match the number of samples (12 in this case).
lapply(regionMat$chr21$bpCoverage[1:2], dim)

## ----'exploreRegMatrix'--------------------------------------------------
## Dimensions of the coverage matrix
dim(regionMat$chr21$coverageMatrix)

## Coverage for each region. This matrix can then be used with limma or other pkgs
head(regionMat$chr21$coverageMatrix)

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
deseq <- regionMat$chr21$regions
mcols(deseq) <- c(mcols(deseq), results(dse))

## Explore the results
deseq

## ----'buildLimmaModels'--------------------------------------------------
## Build models
mod <- model.matrix(~ pheno$group + pheno$gender)
mod0 <- model.matrix(~ pheno$gender)

## ----'transformForLimma'-------------------------------------------------
## Transform coverage
transformedCov <- log2(regionMat$chr21$coverageMatrix + 32)

## ----'identifyDERsLimma'-------------------------------------------------
## Example using limma
library('limma')

## Run limma
fit <- lmFit(transformedCov, mod)
fit0 <- lmFit(transformedCov, mod0)

## Determine DE status for the regions
getF  <- function(fit, fit0, theData) {
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)
	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1","df0")
	fout = data.frame(fout)
	return(fout)
}

ff <- getF(fit, fit0, transformedCov)

## Get the p-value and assign it to the regions
limma <- regionMat$chr21$regions
limma$fstat <- ff$fstat
limma$pvalue <- ff$f_pval
limma$padj <- p.adjust(ff$f_pval, 'BH')

## Explore the results
limma

## ----limmaVSdeseq2-------------------------------------------------------
table(limma$padj < 0.05, deseq$padj < 0.05)

## ----'createMeanBW', eval = !windowsFlag---------------------------------
## Calculate the mean: this step takes a long time with many samples
meanCov <- Reduce('+', fullCov$chr21) / ncol(fullCov$chr21)

## Save it on a bigwig file called meanChr21.bw
createBw(list('chr21' = DataFrame('meanChr21' = meanCov)), keepGR = 
    FALSE)

## ----'railMatrix', eval = !windowsFlag-----------------------------------
## Identify files to use
summaryFile <- 'meanChr21.bw'
## We had already found the sample BigWig files and saved it in the object 'files'
## Lets just rename it to sampleFiles for clarity.
sampleFiles <- files

## Get the regions
system.time( 
    regionMat.rail <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile, 
        sampleFiles = sampleFiles, L = 76, cutoff = 30, maxClusterGap = 3000L)
)

## ----'checkDifferences', eval = !windowsFlag-----------------------------
## Overall not identical due to small rounding errors
identical(regionMat, regionMat.rail)

## Actual regions are the same
identical(ranges(regionMat$chr21$regions), ranges(regionMat.rail$chr21$regions))

## When you round, the small differences go away
identical(round(regionMat$chr21$regions$value, 4), round(regionMat.rail$chr21$regions$value, 4))

identical(round(regionMat$chr21$regions$area, 4), round(regionMat.rail$chr21$regions$area, 4))

## ----'libSize'-----------------------------------------------------------
## Get some idea of the library sizes
sampleDepths <- sampleDepth(collapseFullCoverage(fullCov), 1)
sampleDepths

## ----'makeModels'--------------------------------------------------------
## Define models
models <- makeModels(sampleDepths, testvars = pheno$group, adjustvars = pheno[, c('gender')]) 

## Explore the models
lapply(models, head)

## ----'analyze'-----------------------------------------------------------
## Create a analysis directory
dir.create('analysisResults')
originalWd <- getwd()
setwd(file.path(originalWd, 'analysisResults'))

## Perform differential expression analysis
system.time(results <- analyzeChr(chr = 'chr21', filteredCov$chr21, models, groupInfo = pheno$group, writeOutput = TRUE, cutoffFstat = 5e-02, nPermute = 20, seeds = 20140923 + seq_len(20), returnOutput = TRUE))

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

## ----'featureLevel'------------------------------------------------------
## Get the exon level matrix
system.time(exonCov <- coverageToExon(fullCov, genomicState$fullGenome, L = 76))

## Dimensions of the matrix
dim(exonCov)

## Explore a little bit
tail(exonCov)

## ----'regionMatAnnotate'-------------------------------------------------
## Annotate regions as exonic, intronic or intergenic
system.time(annoGenome <- annotateRegions(regionMat$chr21$regions, genomicState$fullGenome))
## Note that the genomicState object included in derfinder only has information
## for chr21 (hg19).

## Identify closest genes to regions
suppressPackageStartupMessages(library('bumphunter'))
suppressPackageStartupMessages(library('TxDb.Hsapiens.UCSC.hg19.knownGene'))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- annotateTranscripts(txdb)
system.time(annoNear <- matchGenes(regionMat$chr21$regions, genes))

## ----'static-vis', eval = FALSE------------------------------------------
## ## Identify the top regions by highest total coverage
## top <- order(regionMat$chr21$regions$area, decreasing = TRUE)[1:100]
## 
## ## Base-level plots for the top 100 regions with transcript information
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

## ----'exportBigWig', eval = !windowsFlag---------------------------------
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

## ----'exampleNameStyle', eval = FALSE------------------------------------
## ## Set global species and chrsStyle options
## options(species = 'arabidopsis_thaliana')
## options(chrsStyle = 'NCBI')
## 
## ## Then proceed to load and analyze the data

## ----'analyzeNonHuman', eval = FALSE-------------------------------------
## ## Load transcript database information
## library('TxDb.Athaliana.BioMart.plantsmart28')
## 
## ## Set organism options
## options(species = 'arabidopsis_thaliana')
## options(chrsStyle = 'NCBI')
## 
## ## Run command
## analyzeChr(txdb = TxDb.Athaliana.BioMart.plantsmart28, --Other arguments--)

## ----'runExample'--------------------------------------------------------
## Find some regions to work with
example('loadCoverage', 'derfinder')
example('getRegionCoverage', 'derfinder')

## ----'loadWhich'---------------------------------------------------------
## Illustrate reading data from a set of regions
test <- loadCoverage(files = files, chr = '21', cutoff = NULL, which = regions, protectWhich = 0, fileStyle = 'NCBI')

## Some reads were ignored and thus the coverage is lower as can be seen below:
sapply(test$coverage, max) - sapply(genomeDataRaw$coverage, max)

## ----'loadWhich2'--------------------------------------------------------
## Illustrate reading data from a set of regions

test2 <- loadCoverage(files = files, chr = '21', cutoff = NULL, which = 
regions, protectWhich = 3e4, fileStyle = 'NCBI')

## Adding some padding to the regions helps get the same coverage
identical(sapply(test2$coverage, max), sapply(genomeDataRaw$coverage, max))

## A more detailed test reveals that the coverage matches at every base
all(mapply(function(x, y) { identical(x, y) }, test2$coverage, genomeDataRaw$coverage))

## ----'derfinder-analysis', eval = FALSE----------------------------------
## ## Run derfinder's analysis steps with timing info
## 
## ## Load libraries
## library("getopt")
## 
## ## Available at http:/bioconductor.org/packages/derfinder
## library("derfinder")
## 
## ## Specify parameters
## spec <- matrix(c(
## 	'DFfile', 'd', 1, "character", "path to the .Rdata file with the results from loadCoverage()",
## 	'chr', 'c', 1, "character", "Chromosome under analysis. Use X instead of chrX.",
## 	'mcores', 'm', 1, "integer", "Number of cores",
## 	'verbose' , 'v', 2, "logical", "Print status updates",
## 	'help' , 'h', 0, "logical", "Display help"
## ), byrow=TRUE, ncol=5)
## opt <- getopt(spec)
## 
## ## Testing the script
## test <- FALSE
## if(test) {
## 	## Speficy it using an interactive R session and testing
## 	test <- TRUE
## }
## 
## ## Test values
## if(test){
## 	opt <- NULL
## 	opt$DFfile <- "/ProjectDir/derCoverageInfo/chr21Cov.Rdata"
## 	opt$chr <- "21"
## 	opt$mcores <- 1
## 	opt$verbose <- NULL
## }
## 
## ## if help was asked for print a friendly message
## ## and exit with a non-zero error code
## if (!is.null(opt$help)) {
## 	cat(getopt(spec, usage=TRUE))
## 	q(status=1)
## }
## 
## ## Default value for verbose = TRUE
## if (is.null(opt$verbose)) opt$verbose <- TRUE
## 
## if(opt$verbose) message("Loading Rdata file with the output from loadCoverage()")
## load(opt$DFfile)
## 
## ## Make it easy to use the name later. Here I'm assuming the names were generated using output='auto' in loadCoverage()
## eval(parse(text=paste0("data <- ", "chr", opt$chr, "CovInfo")))
## eval(parse(text=paste0("rm(chr", opt$chr, "CovInfo)")))
## 
## ## Just for testing purposes
## if(test) {
## 	tmp <- data
## 	tmp$coverage <- tmp$coverage[1:1e6, ]
## 	library("IRanges")
## 	tmp$position[which(tmp$pos)[1e6 + 1]:length(tmp$pos)] <- FALSE
## 	data <- tmp
## }
## 
## ## Load the models
## load("models.Rdata")
## 
## ## Load group information
## load("groupInfo.Rdata")
## 
## 
## ## Run the analysis with lowMemDir
## analyzeChr(chr=opt$chr, coverageInfo=data, models=models, cutoffFstat=1e-06, cutoffType="theoretical", nPermute=1000, seeds=seq_len(1000), maxClusterGap=3000, groupInfo=groupInfo, subject="hg19", mc.cores=opt$mcores, lowMemDir=file.path(tempdir(), paste0("chr", opt$chr) , "chunksDir")), verbose=opt$verbose, chunksize=1e5)
## 
## ## Done
## if(opt$verbose) {
## 	print(proc.time())
## 	print(sessionInfo(), locale=FALSE)
## }
## 

## ----createVignette, eval=FALSE------------------------------------------
## ## Create the vignette
## library('rmarkdown')
## system.time(render('derfinder-users-guide.Rmd', 'BiocStyle::html_document'))
## 
## ## Extract the R code
## library('knitr')
## knit('derfinder-users-guide.Rmd', tangle = TRUE)

## ----createVignette2-----------------------------------------------------
## Clean up
file.remove('derfinderUsersGuideRed.bib')
unlink('analysisResults', recursive = TRUE)
file.remove('HSB113.bw')
file.remove('meanChr21.bw')

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

