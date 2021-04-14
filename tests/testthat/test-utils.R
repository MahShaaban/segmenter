context("Test Utils")

test_that("tidy_ranges works.", {
    binfile <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18/',
                           'K562_chr11_binary.txt.gz',
                           package = 'segmentr')
    chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
                                 'hg18.txt',
                                 package = 'segmentr')

    bins <- read_bins_file(binfile)
    gr <- range_bins(list(bins), chromsizefile, 200, tidy = FALSE)
    tidy_gr <- tidy_ranges(gr)

    expect_s4_class(tidy_gr, 'GRanges')
    expect_equal(bins$seqname,
                 unique(as.character(GenomicRanges::seqnames(tidy_gr))))

    nms <- colnames(S4Vectors::values(gr))[1:3]
    tidy_gr2 <- tidy_ranges(gr, nms)

    expect_s4_class(tidy_gr2, 'GRanges')
    expect_equal(length(unique(tidy_gr2$group)), length(nms))
})

test_that("range_bins works.", {
    binfile <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18/',
                           'K562_chr11_binary.txt.gz',
                           package = 'segmentr')
    chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
                                 'hg18.txt',
                                 package = 'segmentr')

    bins <- read_bins_file(binfile)
    gr <- range_bins(list(bins), chromsizefile, 200, tidy = FALSE)

    expect_s4_class(gr, 'GRanges')
    expect_equal(bins$seqname,
                 unique(as.character(GenomicRanges::seqnames(gr))))
    expect_equal(nrow(bins$binaries), length(gr))

    # tidy works
    gr_tidy <- range_bins(list(bins), chromsizefile, 200, tidy = TRUE)

    expect_s4_class(gr_tidy, 'GRanges')
    expect_equal(names(GenomicRanges::mcols(gr_tidy)),
                 c('id', 'group', 'value'))

    # returns summarized experiment object
    se <- range_bins(list(bins), chromsizefile, 200,
                     return = 'SummarizedExperiment')

    expect_s4_class(se, 'SummarizedExperiment')
    expect_equal(nrow(se), length(gr))
})

test_that("range_counts works.", {
    bam_file <- system.file("extdata", "randomBam.bam", package = "bamsignals")
    rand_anno <- system.file("extdata", "randomAnnot.Rdata", package = "bamsignals")
    features <- GenomicRanges::promoters(get(load(rand_anno)))

    reads <- cbind(read_bam_file(bam_file, features),
                   read_bam_file(bam_file, features),
                   read_bam_file(bam_file, features))

    gr <- range_counts(reads, features)

    expect_s4_class(gr, 'GRanges')
    expect_equal(length(gr), length(features))
    expect_equal(ncol(S4Vectors::values(gr)), ncol(reads))

    expect_error(range_counts(reads, features, average = TRUE))

    marks <- paste0('mark', c(1, 1, 2))
    gr2 <- range_counts(reads, features, average = TRUE, marks = marks)

    expect_s4_class(gr2, 'GRanges')
    expect_equal(length(gr2), length(features))
    expect_equal(ncol(S4Vectors::values(gr2)), length(unique(marks)))
    expect_true(all(names(S4Vectors::values(gr2)) == unique(marks)))

    gr3 <- range_counts(reads, features, tidy = TRUE, average = TRUE,
                        marks = marks)

    expect_s4_class(gr3, 'GRanges')
    expect_equal(length(gr3), length(features) * length(unique(marks)))
    expect_true(all(gr3$group %in% marks))

    se <- range_counts(reads, features, return = 'SummarizedExperiment',
                       tidy = TRUE, average = TRUE, marks = marks)

    expect_s4_class(se, 'SummarizedExperiment')
})

test_that("count_reads_ranges works.", {
    inputbamdir <- system.file("extdata", package="bamsignals")

    cellmarkfiletable <- system.file('extdata/input',
                                     'cell_mark_table.tsv',
                                     package = 'segmentr')

    rand_anno <- system.file("extdata",
                             "randomAnnot.Rdata",
                             package = "bamsignals")
    segments <- GenomicRanges::promoters(get(load(rand_anno)))
    segments <- split(segments, rep(c('cell1', 'cell2')), 10)

    gr <- count_reads_ranges(segments, cellmarkfiletable, inputbamdir)

    expect_true(is.list(gr))
    expect_equal(length(gr), length(segments))
    expect_true(all(sapply(gr, nrow) == lengths(segments)))
})

test_that("merge_segments_bins works.", {
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
                       binsize = 200,
                       read_bins = TRUE)

    se <- merge_segments_bins(segment(obj)[[1]],
                              bins(obj)[[1]])

    expect_s4_class(se, 'SummarizedExperiment')
})
