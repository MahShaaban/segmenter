## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
#  # install from bioconductor
#  BiocManager::install('segmenter')
#  
#  # install from github
#  remotes::install_github('MahShaaban/segmenter@devel')

## ----load_library, eval=FALSE-------------------------------------------------
#  # load library
#  library(segmenter)

## ----prepare_directories, message = FALSE, eval=FALSE-------------------------
#  # locate input and annotation files
#  inputdir <- system.file('extdata/SAMPLEDATA_HG18',
#                          package = 'segmenter')
#  coordsdir <- system.file('extdata/COORDS',
#                           package = 'chromhmmData')
#  anchorsdir <- system.file('extdata/ANCHORFILES',
#                            package = 'chromhmmData')
#  chromsizefile <- system.file('extdata/CHROMSIZES',
#                               'hg18.txt',
#                               package = 'chromhmmData')

## ----getting_stated, eval=FALSE-----------------------------------------------
#  # run command
#  obj <- learn_model(inputdir = inputdir,
#                     coordsdir = coordsdir,
#                     anchorsdir = anchorsdir,
#                     chromsizefile = chromsizefile,
#                     numstates = 3,
#                     assembly = 'hg18',
#                     cells = c('K562', 'GM12878'),
#                     annotation = 'RefSeq',
#                     binsize = 200)

## ----show, eval=FALSE---------------------------------------------------------
#  # show the object
#  show(obj)

## ----setup--------------------------------------------------------------------
# load required libraries
library(segmenter)
library(Gviz)
library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)

## ----genomic_annotations------------------------------------------------------
# coordinates
coordsdir <- system.file('extdata/COORDS',
                         package = 'chromhmmData')

list.files(file.path(coordsdir, 'hg18'))

# anchors
anchorsdir <- system.file('extdata/ANCHORFILES',
                          package = 'chromhmmData')

list.files(file.path(anchorsdir, 'hg18'))

# chromosomes' sizes
chromsizefile <- system.file('extdata/CHROMSIZES',
                             'hg18.txt',
                              package = 'chromhmmData')

readLines(chromsizefile, n = 3)

## ----cellmarkfiletable--------------------------------------------------------
# a table to assign marker and cell names to the bam files
cellmarkfiletable <- system.file('extdata',
                                 'cell_mark_table.tsv',
                                 package = 'segmenter')

readLines(cellmarkfiletable, n = 3)

## ----binary_inputs------------------------------------------------------------
# locate input and output
inputdir <- system.file("extdata", package = "bamsignals")
outputdir <- tempdir()

