context("Test Visualization")

test_that("annotate_segments works.", {
    numstates <- 3
    cells <- c('K562', 'GM12878')

    library(TxDb.Hsapiens.UCSC.hg18.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene

    segs <- segment(test_obj)
    segs_annotated <- annotate_segments(segs, TxDb = txdb, verbose = FALSE)

    expect_s4_class(segs_annotated[[1]], 'GRanges')
    expect_true(ncol(GenomicRanges::mcols(segs_annotated[[1]])) > 1)
    expect_equal(length(segs_annotated[[1]]), length(segs[[1]]))
})

test_that("get_frequency works.", {
    numstates <- 3
    cells <- c('K562', 'GM12878')

    segs <- segment(test_obj)
    freq <- get_frequency(segs, normalize = TRUE, tidy = TRUE)

    expect_true(is.data.frame(freq))
    expect_equal(ncol(freq), 3)
    expect_equal(nrow(freq), numstates * length(cells))
    expect_true(all(cells %in% freq$cell))
    expect_true(all(segs[[1]]$state %in% freq$state))
    expect_equal(sum(freq$frequency), length(cells))
})

test_that("get_width works.", {
    numstates <- 3
    cells <- c('K562', 'GM12878')

    segs <- segment(test_obj)
    wdth <- get_width(segs)

    expect_true(is.data.frame(wdth))
    expect_equal(ncol(wdth), 3)
    expect_equal(nrow(wdth), sum(lengths(segs)))
    expect_true(all(cells %in% wdth$cell))
    expect_true(all(segs[[1]]$state %in% wdth$state))

    wdth2 <- get_width(segs, average = TRUE)

    expect_true(is.data.frame(wdth2))
    expect_equal(ncol(wdth2), 3)
    expect_equal(nrow(wdth2), numstates * length(cells))
    expect_true(all(cells %in% wdth2$cell))
    expect_true(all(segs[[1]]$state %in% wdth2$state))
})
