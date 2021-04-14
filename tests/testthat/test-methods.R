context("Test S4 Methods")

test_that("accessors works", {
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
    numstates <- 3
    cells <- c('K562', 'GM12878')

    obj <- learn_model(inputdir = inputdir,
                       outputdir = outputdir,
                       coordsdir = coordsdir,
                       anchorsdir = anchorsdir,
                       chromsizefile = chromsizefile,
                       numstates = 3,
                       assembly = 'hg18',
                       cells = c('K562', 'GM12878'),
                       annotation = 'RefSeq',
                       binsize = 200,
                       read_bins = TRUE)

    expect_s4_class(obj, 'segmentation')

    expect_true(is.list(model(obj)))

    expect_true(is.matrix(emission(obj)))
    expect_equal(nrow(emission(obj)), numstates)

    expect_true(is.matrix(transition(obj)))
    expect_equal(nrow(transition(obj)), numstates)
    expect_equal(ncol(transition(obj)), numstates)

    expect_true(is.list(overlap(obj)))
    expect_equal(length(overlap(obj)), length(cells))
    expect_equal(names(overlap(obj)), cells)

    expect_true(is.list(TSS(obj)))
    expect_equal(length(TSS(obj)), length(cells))
    expect_equal(names(TSS(obj)), cells)

    expect_true(is.list(TES(obj)))
    expect_equal(length(TES(obj)), length(cells))
    expect_equal(names(TES(obj)), cells)

    expect_true(is.list(segment(obj)))
    expect_equal(length(segment(obj)), length(cells))
    expect_equal(names(segment(obj)), cells)

    expect_true(is.list(bins(obj)))
    expect_equal(length(bins(obj)), length(cells))
    expect_equal(names(bins(obj)), cells)

    expect_true(is.character(cells(obj)))
    expect_true(is.character(markers(obj)))
    expect_true(is.character(states(obj)))

    expect_true(is.numeric(likelihood(obj)))
})
