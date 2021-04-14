#' Tidy the metadata of a \code{GRanges} object
#'
#' @param gr A \code{GRanges} object
#' @param columns A \code{character} vectors. The names of columns to be tidied.
#' @param low An \code{integer}. All values <= this \code{integer} will be
#' removed.
#'
#' @return A \code{GRanges} object
#'
#' @importFrom S4Vectors values 'values<-'
#' @importFrom utils stack
#'
#' @export
tidy_ranges <- function(gr, columns, low = 0) {
    # get metadata into a data.frame
    df <- as.data.frame(values(gr))

    # index columns, otherwise, use all
    if (missing(columns)) ind <- names(df)
    else ind <- intersect(names(df), columns)

    # stack columns
    df2 <- stack(subset(df, select = ind))
    df2$id <- rep(nrow(df), length(ind))
    df2 <- subset(df2, select = c('id', 'ind', 'values'))
    names(df2) <- c('id', 'group', 'value')

    # filter out zeros
    df2 <- df2[df2$value > low, ]

    # make a new gr and add the gathered column column
    res <- gr[df2$id]

    ind2 <- setdiff(names(values(res)), ind)
    values(res) <- cbind(values(res)[, ind2], df2)

    # return
    return(res)
}

#' Format the loaded binarized data
#'
#' The function takes the \code{data.frame}s of the loaded binarized data files
#' and format them into \code{GRanges} or \code{SummarizedExperiment} objects.
#'
#' @param bins A \code{list} of the \code{read_bins_file} output.
#' @param chromsizefile A string. The path to the chromosomes sizes file.
#' @param binsize An integer. The number in bp used to generate binarized files.
#' @param return A string. Possible values are \code{GRanges} (default) or
#' \code{SummarizedExperiment}.
#' @param tidy A \code{logical}. Default is \code{TRUE}. Whether to tidy the
#' metadata columns of the \code{GRanges} object.
#'
#' @return \code{GRanges} (default) or \code{SummarizedExperiment}.
#'
#' @importFrom GenomicRanges tileGenome mcols bindROWS GRanges granges
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors values 'values<-'
#'
#' @export
range_bins <- function(bins, chromsizefile, binsize, return = 'GRanges',
                       tidy = TRUE) {
    # make a named vector of chrom lengths
    chromsizes <- read_chromsize_file(chromsizefile)
    seqname_vec <- chromsizes$width
    names(seqname_vec) <- chromsizes$seqnam

    # loop over bins and add as metadata to a gr of the same chromosome
    grs <- lapply(bins, function(x) {
        # get the chrom size and name
        lngth <- seqname_vec[[x$seqname]]
        names(lngth) <- x$seqname

        # tile genome into bins of binsize
        tiles <- tileGenome(
            seqlengths = lngth,
            tilewidth = binsize,
            cut.last.tile.in.chrom = TRUE
        )

        # subset tiles to the number of bins
        tiles <- tiles[seq_len(nrow(x$binaries))]

        # add bins to tiles
        values(tiles) <- x$binaries

        tiles
    })

    # bind granges together in one big range
    res <- bindROWS(GRanges(), objects = grs)

    # make a return object
    if (return == 'SummarizedExperiment') {
        res <- SummarizedExperiment(
            assays = list(binaries = values(res)),
            rowRanges = granges(res)
        )
        return(res)
    } else if (return == 'GRanges') {
        # return tiles
        if (tidy) res <- tidy_ranges(res)

        return(res)
    } else {
        stop()
    }
}

#' Format the loaded counts data
#'
#' The function takes the \code{data.frame}s of the loaded counts data and
#' format them into \code{GRanges} or \code{SummarizedExperiment} objects.
#'
#' @param counts A \code{matrix} of the \code{read_bam_file} output.
#' @param features A \code{GRanges}. That was used to count the bam files.
#' @param return A string. Possible values are \code{GRanges} (default) or
#' \code{SummarizedExperiment}.
#' @param tidy A \code{logical}. Default is \code{TRUE}. Whether to tidy the
#' metadata columns of the \code{GRanges} object.
#' @param average A \code{logical}. Default is \code{FALSE}. Whether to average
#' the counts by \code{marks} before building the object.
#' @param marks A \code{character} vector. The length shoud equal the numbe of
#' columns in \code{counts} and is used for averaging and renaming the matrix
#' columns.
#'
#' @return \code{GRanges} (default) or \code{SummarizedExperiment}.
#'
#' @importFrom stats aggregate
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors values 'values<-'
#' @importFrom GenomicRanges mcols
#'
#' @export
range_counts <- function(counts, features, return = 'GRanges', tidy = FALSE,
                         average = FALSE, marks) {
    # make a matrix of counts
    stopifnot(is.matrix(counts))

    # average and change names if true
    if (average) {
        # stop if marks are not provided or not equal the nrow m
        stopifnot(!missing(marks))
        stopifnot(length(marks) == ncol(counts))

        # average by mark
        df <- aggregate(t(counts),
                        by = list(mark = marks),
                        FUN = function(x) round(mean(x)))

        # make a matrix
        counts <- as.matrix(t(df[, -1]))
        colnames(counts) <- unique(marks)
        rownames(counts) <- NULL
    }

    if (return == 'SummarizedExperiment') {
        res <- SummarizedExperiment(
            assays = list(counts = counts),
            rowRanges = features
        )

        return(res)

    } else if (return == 'GRanges') {
        # add to grange
        gr <- features
        values(gr) <- cbind(mcols(gr), as.data.frame(counts))

        if (tidy) {
            stopifnot(!missing(marks))
            stopifnot(average)

            gr <- tidy_ranges(gr)
        }

        # return
        return(gr)
    }
}

#' Count reads in \code{GRanges} objects from bam files
#'
#' @param ranges A \code{GRanges} to count in.
#' @param cellmarkfiletable A string. The path to the input files table.
#' @param inputbamdir A \code{string}. The path to the input bam files
#' directory.
#'
#' @return A \code{SummarizedExperiment} object with \code{ranges} as its
#' \code{rowRanges} and the counts as the \code{assay}.
#'
#' @export
count_reads_ranges <- function(ranges, cellmarkfiletable, inputbamdir) {
    # read input_file
    input <- read_cellmark_file(cellmarkfiletable)

    # make a list of inputs
    input_list <- split(input, input$cell)
    stopifnot(all(names(ranges) %in% names(input_list)))

    # reorder input_list using ranges names
    input_list <- input_list[names(ranges)]

    # loop over lists
    reads <- mapply(function(x, y) {
        # make files names from table and inputbamdir
        fls <- file.path(inputbamdir, x$file)

        # read each file using corresponding ranges
        m <- lapply(fls, function(f) {
            read_bam_file(f, y)
        })

        # bind as matrix
        counts <- do.call(cbind, m)

        # make granges objects
        range_counts(counts, y, return = 'SummarizedExperiment',
                     average = TRUE, marks = x$mark)
    }, input_list, ranges, SIMPLIFY = FALSE)

    # return
    return(reads)
}

#' Merge segments and bins objects
#'
#' @param segments A \code{GRanges} object. Usually the output of calling
#' \code{segment} on the the output object of \code{lean_model}.
#' @param bins A \code{SummarizedExperiment} object. Usually the output of
#'  calling \code{bins} on the the output object of \code{lean_model}.
#'
#' @return A \code{SummarizedExperiment} object with the segment assignment
#' added to the metadata of the \code{rowRanges}.
#'
#' @importFrom SummarizedExperiment rowRanges SummarizedExperiment assay
#' @importFrom IRanges mergeByOverlaps
#' @importFrom S4Vectors values 'values<-'
#'
#' @export
merge_segments_bins <- function(segments, bins) {
    # extract row ranges of bins
    bins_ranges <- rowRanges(bins)
    bins_ranges$id <- seq_len(length(bins_ranges))

    # merge segments and bins row ranges of bins
    df <- mergeByOverlaps(
        bins_ranges,
        segments
    )

    # make a new granges with the overlaps
    gr <- df$bins_ranges
    values(gr) <- cbind(
        values(df$bins_ranges),
        values(df$segments)
    )

    # make a new summarized experiment object
    se <- SummarizedExperiment(
        assays = list(binaries = assay(bins)[df$bins_ranges$id,]),
        rowRanges = gr
    )

    # return
    return(se)
}
