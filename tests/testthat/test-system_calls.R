context("Test System Calls")

test_that("LearnModel works", {
    inputdir <- system.file('extdata/SAMPLEDATA_HG18',
                            package = 'segmenter')
    outputdir <- tempdir()
    coordsdir <- system.file('extdata/COORDS',
                             package = 'chromhmmData')
    anchorsdir <- system.file('extdata/ANCHORFILES',
                              package = 'chromhmmData')
    chromsizefile <- system.file('extdata/CHROMSIZES',
                                 'hg18.txt',
                                  package = 'chromhmmData')

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
