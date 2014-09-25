
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
write.bibtex(c(knitcitations = citation('knitcitations'),
    derfinder = citation('derfinder'), 
    knitrBootstrap = citation('knitrBootstrap'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown')),
    file = 'derfinderAdvRef.bib')
bib <- read.bibtex('derfinderAdvRef.bib')

## Fix some names to work with CRAN and GitHub versions
names(bib)[names(bib) == 'hester2013knitrbootstrap'] <- 'hester2014knitrbootstrap'


## ----'start', message=FALSE----------------------------------------------
## Load libraries
library('derfinder')


## ----'advancedArg'-------------------------------------------------------
## URLs to advanced arguemtns
sapply(c('analyzeChr', 'loadCoverage'), advancedArg, browse = FALSE)
## Set browse = TRUE if you want to open them in your browser


## ----'runExample', bootstrap.show.output=FALSE, bootstrap.show.message=FALSE----
## Find some regions to work with
example('loadCoverage', 'derfinder')
example('getRegionCoverage', 'derfinder')


## ----'loadWhich', bootstrap.show.message=FALSE---------------------------
## Illustrate reading data from a set of regions
test <- loadCoverage(files = files, chr = '21', cutoff = NULL, which = regions, protectWhich = 0)

## Some reads were ignored and thus the coverage is lower as can be seen below:
sapply(test$coverage, max) - sapply(genomeDataRaw$coverage, max)


## ----'loadWhich2', bootstrap.show.message=FALSE--------------------------
## Illustrate reading data from a set of regions
test2 <- loadCoverage(files = files, chr = '21', cutoff = NULL, which = regions, protectWhich = 3e4)

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
## ## Available at https://github.com/lcolladotor/derfinder
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
## analyzeChr(chr=opt$chr, coverageInfo=data, models=models, cutoffFstat=1e-06, cutoffType="theoretical", nPermute=1000, seeds=seq_len(1000), maxClusterGap=3000, groupInfo=groupInfo, subject="hg19", mc.cores=opt$mcores, lowMemDir=paste0("chr", opt$chr, "/chunksDir"), verbose=opt$verbose, chunksize=1e5)
## 
## ## Done
## if(opt$verbose) {
## 	print(proc.time())
## 	print(sessionInfo(), locale=FALSE)
## }
## 


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
##     system.time(knit_bootstrap('derfinderAdvanced.Rmd', chooser=c('boot', 'code'), show_code = TRUE))
##     unlink('derfinder.md')
## } else {
##     ## GitHub version
##     library('rmarkdown')
##     system.time(render('derfinderAdvanced.Rmd', 'knitrBootstrap::bootstrap_document'))
## }
## ## Note: if you prefer the knitr version use:
## # library('rmarkdown')
## # system.time(render('derfinderAdvanced.Rmd', 'html_document'))
## ## Clean up
## file.remove('derfinderAdvRef.bib')
## 
## ## Extract the R code
## library('knitr')
## knit('derfinderAdvanced.Rmd', tangle = TRUE)


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


## ----vignetteBiblio, results='asis', echo=FALSE--------------------------
## Print bibliography
bibliography()

