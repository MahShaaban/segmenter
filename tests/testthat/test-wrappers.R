context("Test Wrappers")

test_that("learn_model works", {
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
})

test_that("can make file names.", {
    outputdir <- system.file('extdata/output', package = 'segmentr')

    expect_true(file.exists(file.path(outputdir, model_file(3))))
    expect_true(file.exists(file.path(outputdir, emissions_file(3))))
    expect_true(file.exists(file.path(outputdir, transitions_file(3))))
    expect_true(file.exists(file.path(outputdir, overlap_files(3, 'K562'))))
    expect_true(file.exists(file.path(outputdir, segments_files(3, 'K562'))))
    expect_true(file.exists(file.path(outputdir, enrichment_files(3, 'K562'))))
})
