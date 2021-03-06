<a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/derfinder.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="https://bioconductor.org/packages/stats/bioc/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/downloads/derfinder.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/posts/derfinder.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/derfinder.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

Status: Travis CI [![Build Status](https://travis-ci.org/lcolladotor/derfinder.svg?branch=master)](https://travis-ci.org/lcolladotor/derfinder),
Bioc-release <a href="http://www.bioconductor.org/packages/release/bioc/html/derfinder.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/release/derfinder.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/build/release/bioc/derfinder.svg" title="build results; click for full report"></a>,
Bioc-devel <a href="http://www.bioconductor.org/packages/devel/bioc/html/derfinder.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/devel/derfinder.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/derfinder/"><img border="0" src="http://www.bioconductor.org/shields/build/devel/bioc/derfinder.svg" title="build results; click for full report"></a>.

Bioc-release <a href="https://bioconductor.org/developers/how-to/unitTesting-guidelines/#coverage"><img border="0" src="http://www.bioconductor.org/shields/coverage/release/derfinder.svg" title="Test coverage percentage, or 'unknown'"></a>, Bioc-devel <a href="https://codecov.io/github/Bioconductor-mirror/derfinder?branch=master"><img border="0" src="http://www.bioconductor.org/shields/coverage/devel/derfinder.svg" title="Test coverage percentage, or 'unknown'"></a>, Codecov [![codecov.io](https://codecov.io/github/lcolladotor/derfinder/coverage.svg?branch=master)](https://codecov.io/github/lcolladotor/derfinder?branch=master)

derfinder
=========

Annotation-agnostic differential expression analysis of RNA-seq data at base-pair resolution via the DER Finder approach. This package contains two different implementations of this approach. The first one is the single base-level F-statistics implementation and the second one is via identifying expressed regions. For more information about `derfinder` check the vignettes [here](http://www.bioconductor.org/packages/derfinder).


# Further documentation

You can generate HTML reports from the results using __regionReport__ 
available [here](https://github.com/lcolladotor/regionReport).

# Installation instructions

Get R 3.2.0 from [CRAN](http://cran.r-project.org/).

```R
## From Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite('derfinder')

## Suggested:
biocLite(c('derfinderPlot', 'regionReport'))
```

# Vignettes

The vignettes for this package can be viewed [here](http://lcolladotor.github.io/derfinder/) or via [Bioconductor's website](http://www.bioconductor.org/packages/derfinder).

# 'Watch' for updates

This software is in development, so we highly recommend 'watching' the 
repository: Click on the top right under `Watch`. You will then receive 
notifications for issues, comments, and pull requests as described 
[here](https://help.github.com/articles/notifications).

You will need a GitHub account to be able to `Watch` the repository.

# Citation

Below is the citation output from using `citation('derfinder')` in R. Please 
run this yourself to check for any updates on how to cite __derfinder__.

To cite package __derfinder__ in publications use:

Collado-Torres L, Nellore A, Frazee AC, Wilks C, Love MI, Irizarry RA, Leek JT and Jaffe AE (2016). "Flexible expressed region analysis for RNA-seq with derfinder". _bioRxiv_. <URL: http://dx.doi.org/10.1101/015370>, <URL:
http://biorxiv.org/content/early/2016/05/07/015370>.

Frazee AC, Sabunciyan S, Hansen KD, Irizarry RA and Leek JT (2014). “Differential expression analysis of RNA-seq data at
single-base resolution.” _Biostatistics_, *15 (3)*, pp. 413-426. <URL: http://dx.doi.org/10.1093/biostatistics/kxt053>, <URL:
http://biostatistics.oxfordjournals.org/content/15/3/413.long>.

A BibTeX entry for LaTeX users is

@Manual{,
    title = {Flexible expressed region analysis for RNA-seq with derfinder},
    author = {Leonardo Collado-Torres and Abhinav Nellore and Alyssa C. Frazee and Christopher Wilks and Michael I. Love and Rafael A. Irizarry and Jeffrey T. Leek and Andrew E. Jaffe},
    year = {2016},
    journal = {bioRxiv},
    doi = {10.1101/015370},
    url = {http://biorxiv.org/content/early/2016/05/07/015370}
}

@Article{,
    title = {Differential expression analysis of RNA-seq data at single-base resolution},
    author = {Alyssa C. Frazee and Sarven Sabunciyan and Kasper D. Hansen and Rafael A. Irizarry and Jeffrey T. Leek},
    year = {2014},
    journal = {Biostatistics},
    volume = {15 (3)},
    pages = {413-426},
    doi = {10.1093/biostatistics/kxt053},
    url = {http://biostatistics.oxfordjournals.org/content/15/3/413.long},
}

# DER Finder versions

* The original implementation of the DER Finder approach as published in Frazee et al, Biostatistics 2014 is available via GitHub at [derfinder](https://github.com/leekgroup/derfinder).
* The version implementing the single base-level approach via calculating F-stastics as described in the pre-print Collado-Torres et al, bioRxiv 2015 is available via Bioconductor at [derfinder](http://bioconductor.org/packages/derfinder). The same package has the functions required for the expressed regions-level approach.

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/) as well as Bioconductor's nightly build.
