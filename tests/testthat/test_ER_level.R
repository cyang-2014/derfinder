context('ER-level approach')

## Create some toy data
library('IRanges')
set.seed(20160606)
x <- Rle(round(runif(1e4, max=10)))
y <- Rle(round(runif(1e4, max=10)))
z <- Rle(round(runif(1e4, max=10)))
fullCov <- list('chr21' = DataFrame(x, y, z))

## Calculate a proxy of library size
libSize <- sapply(fullCov$chr21, sum)

## Run region matrix normalizing the coverage
regionMat <- regionMatrix(fullCov = fullCov, maxRegionGap = 10L,
    maxClusterGap = 300L, L = 36, totalMapped = libSize, targetSize = 4e4)

## You can alternatively use filterData() on fullCov to reduce the required
## memory before using regionMatrix(). This can be useful when mc.cores > 1
filteredCov <- lapply(fullCov, filterData, returnMean=TRUE, filter='mean', 
    cutoff=5, totalMapped = libSize, targetSize = 4e4)
regionMat2 <- regionMatrix(filteredCov, maxRegionGap = 10L,
    maxClusterGap = 300L, L = 36, runFilter=FALSE)
    
regionMat3 <- regionMatrix(fullCov = fullCov, maxRegionGap = 10L, 
    maxClusterGap = 30000L, L = 36, totalMapped = libSize, targetSize = 4e4)
    
## fullCov with position info
fullCov2 <- list('chr21' = 
    list('coverage' = DataFrame(x, y, z), 'position' = NULL))
regionMat4 <- regionMatrix(fullCov = fullCov2, maxRegionGap = 10L,
    maxClusterGap = 300L, L = 36, totalMapped = libSize, targetSize = 4e4)

## fullCov with non-NULL position info
fullCov3 <- list('chr21' = 
    list('coverage' = DataFrame(x, y, z), 'position' = Rle(rep(TRUE, 1e4))))
regionMat5 <- regionMatrix(fullCov = fullCov3, maxRegionGap = 10L,
    maxClusterGap = 300L, L = 36, totalMapped = libSize, targetSize = 4e4)    
    
    
## Filter
filtered <- filterData(fullCov3[[1]]$coverage, index = fullCov3[[1]]$position, returnMean = TRUE, cutoff = 5)

test_that('regionMatrix', {
    expect_equal(regionMat, regionMat2)
    expect_equal(sum(regionMat3$chr21$regions$cluster == 1),
        length(regionMat3$chr21$regions))
    expect_equal(regionMat, regionMat4)
    expect_equal(regionMat, regionMat5)
    expect_lt(length(filtered$meanCoverage), 1e4)
    expect_equal(length(filtered$meanCoverage), sum(filtered$position))
})

if(.Platform$OS.type != 'windows') {
    ## Get data
    library('derfinderData')
    
    ## Identify sample files
    sampleFiles <- rawFiles(system.file('extdata', 'AMY', package = 
        'derfinderData'), samplepatt = 'bw', fileterm = NULL)
    names(sampleFiles) <- gsub('.bw', '', names(sampleFiles))
    
    ## Create the mean bigwig file. This file is normally created by Rail
    ## but in this example we'll create it manually.
    library('GenomicRanges')
    railCov <- fullCoverage(files = sampleFiles, chrs = 'chr21')
    meanCov <- Reduce('+', railCov$chr21) / ncol(railCov$chr21)
    createBw(list('chr21' = DataFrame('meanChr21' = meanCov)), keepGR = 
        FALSE)
    
    summaryFile <- 'meanChr21.bw'
    
    ## Get the regions
    railMat <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile, 
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1, maxClusterGap = 3000L)
        
    ## Changing chunksize
    railMat2 <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile, 
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1, maxClusterGap = 3000L,
        chunksize = 100)
    
    ## Changing maxClusterGap
    railMat3 <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile, 
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1, maxClusterGap = 1e6)
        
    ## Reproducing results with regionMatrix
    railMat4 <- regionMatrix(fullCov = railCov, L = 76, cutoff = 5.1,
        maxClusterGap = 3000L)
    
    ## Smoothing
    meanCov2 <- meanCov[9820000:9830000, ]
    createBw(list('chr21' = DataFrame('mean2Chr21' = meanCov2)), keepGR = 
        FALSE)
    summaryFile2 <- 'mean2Chr21.bw'
    
    ## First with the less memory intensive of the two methods
    railMat5 <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile2, 
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1, maxClusterGap = 3000L,
        smoothMean = TRUE, smoothFunction = bumphunter::runmedByCluster,
        k = 299)
    
    ## Next with the more memory intensive one
    railMat6 <- railMatrix(chrs = 'chr21', summaryFiles = summaryFile2, 
        sampleFiles = sampleFiles, L = 76, cutoff = 5.1, maxClusterGap = 3000L,
        smoothMean = TRUE, minInSpan = 76, minNum = 76, bpSpan = 300)
    
    ## Reproducing results with regionMatrix
    railCov2 <- railCov
    railCov2$chr21 <- railCov2$chr21[9820000:9830000, ]
    
    ## Note that railMat 5 and 6 don't have the correct coverageMatrix since
    ## the sampleFiles were not changed, only the mean was changed.
    railMat7 <- regionMatrix(fullCov = railCov2, L = 76, cutoff = 5,
        maxClusterGap = 3000L, smoothMean = TRUE,
        smoothFunction = bumphunter::runmedByCluster, k = 299, returnBP = FALSE)
        
        
    test_that('railMatrix', {
        expect_equal(railMat, railMat2)
        expect_lt(max(railMat3$chr21$regions$cluster), max(railMat$chr21$regions$cluster))
        expect_equivalent(railMat$chr21$regions, railMat4$chr21$regions)
        expect_equal(railMat$chr21$coverageMatrix, railMat4$chr21$coverageMatrix, tolerance = 0.005)
        expect_lt(railMat5$chr21$regions$value, railMat$chr21$regions[1]$value)
        expect_lt(railMat6$chr21$regions$value, railMat5$chr21$regions$value)
        expect_equal(abs(railMat6$chr21$regions$value * width(railMat6$chr21$regions)), railMat6$chr21$regions$area)
        expect_equal(railMat5$chr21$regions, railMat7$chr21$regions, tolerance = 0.0000001)
    })
}
