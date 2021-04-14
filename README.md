![Basic Checks](https://github.com/MahShaaban/segmentr/workflows/Basic%20Checks/badge.svg)
![Build and Push Docker Image](https://github.com/MahShaaban/segmentr/workflows/Build%20and%20Push%20Docker%20Image/badge.svg)

# segmentr

Call chromHMM from within R, capture the output files in an R object and
interface to other relevant analysis tools in R.

## Installation

You can install the released version of segmentr from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("segmentr")
```

## Example

```r
# load required libraries
library(segmentr)
```

```r
# locate input and output files
inputdir <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18',
                        package = 'segmentr')
outputdir <- tempdir()
coordsdir <- system.file('extdata/ChromHMM/COORDS',
                         package = 'segmentr')
anchorsdir <- system.file('extdata/ChromHMM/ANCHORFILES',
                          package = 'segmentr')
chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
                             'hg18.txt',
                             package = 'segmentr')

# run command
obj <- learn_model(inputdir = inputdir,
                   outputdir = outputdir,
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


