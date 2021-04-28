context("Test Wrappers")

test_that("learn_model works", {
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

    expect_s4_class(obj, 'segmentation')
    
    expect_true(file.exists(file.path(outputdir, model_file(3))))
    expect_true(file.exists(file.path(outputdir, emissions_file(3))))
    expect_true(file.exists(file.path(outputdir, transitions_file(3))))
    expect_true(file.exists(file.path(outputdir, overlap_files(3, 'K562'))))
    expect_true(file.exists(file.path(outputdir, segments_files(3, 'K562'))))
    expect_true(file.exists(file.path(outputdir, enrichment_files(3, 'K562'))))
})
