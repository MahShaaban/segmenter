![Basic Checks](https://github.com/MahShaaban/segmenter/workflows/Basic%20Checks/badge.svg)
![Build and Push Docker Image](https://github.com/MahShaaban/segmenter/workflows/Build%20and%20Push%20Docker%20Image/badge.svg)

# segmenter

Perform Chromatin Segmentation Analysis in R by Calling ChromHMM

## Description

Call *chromHMM* from within *R*, capture the output files in an `S4`  object and 
interface to other relevant Bioconductor analysis tools. In addition `segmenter` 
provides functions to test, select and visualize the output of the segmentation.

## Installation

You can install the released version of `segmenter` from [Bioconductor](https://bioconductor.org/) with:

``` r
BiocManager::install("segmenter")
```

or the development version from GitHub with:

```r
remotes::install_github('MahShaaban/segmenter@devel')
```

## Getting started

```r
# load required libraries
library(segmenter)
```

```r
# locate input and annotation files
inputdir <- system.file('extdata/SAMPLEDATA_HG18',
                        package = 'segmenter')
                        
coordsdir <- system.file('extdata/COORDS',
                         package = 'chromhmmData')
anchorsdir <- system.file('extdata/ANCHORFILES',
                          package = 'chromhmmData')
chromsizefile <- system.file('extdata/CHROMSIZES',
                             'hg18.txt',
                             package = 'chromhmmData')
```

```r
# run command
obj <- learn_model(inputdir = inputdir,
                   coordsdir = coordsdir,
                   anchorsdir = anchorsdir,
                   chromsizefile = chromsizefile,
                   numstates = 3,
                   assembly = 'hg18',
                   cells = c('K562', 'GM12878'),
                   annotation = 'RefSeq',
                   binsize = 200)
```

```r
# show the object
show(obj)
```

```r
# access object slots
emission(obj)
```

```r
# operate on the object
get_frequency(segments = segment(obj))
```
