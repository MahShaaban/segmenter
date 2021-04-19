context("Test System Calls")

test_that("LearnModel works", {
    inputdir <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18',
                            package = 'segmenter')
    outputdir <- tempdir()
    coordsdir <- system.file('extdata/ChromHMM/COORDS',
                             package = 'segmenter')
    anchorsdir <- system.file('extdata/ChromHMM/ANCHORFILES',
                              package = 'segmenter')
    chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
                                 'hg18.txt',
                                  package = 'segmenter')

    expect_null(
        .LearnModel(inputdir = inputdir,
                    outputdir = outputdir,
                    coordsdir = coordsdir,
                    anchorsdir = anchorsdir,
                    chromsizefile = chromsizefile,
                    numstates = 3,
                    assembly = 'hg18')
    )
})
