#' Annotate segments
#'
#' Annotate the \code{GRanges} objects of the segments using
#' \code{\link[ChIPseeker]{annotatePeak}} (see for details)
#' @param segments A \code{GRanges} object. Usually the output of calling
#' \code{segment} on the the output object of \code{lean_model}.
#' @param ... Other arguments passed to \code{\link[ChIPseeker]{annotatePeak}}
#'
#' @return A \code{GRanges} object which is identical to the input in addition
#' to the annotations as metadata columns.
#'
#' @importFrom ChIPseeker annotatePeak as.GRanges
#'
#' @export
annotate_segments <- function(segments, ...) {
    # apply annotate peaks to items of segments list
    grs <- lapply(segments,
                  function(x) {
                      gr <- annotatePeak(x, ...)
                      as.GRanges(gr)
                  })

    # return
    return(grs)
}

#' Get the frequency of the segments in each cell type
#'
#' @inheritParams annotate_segments
#' @param normalize A logical. Whether the frequency should be normalized by
#' the total number of segments
#' @param tidy A logical.
#' @param plot A logical.
#' @param ... Other arguments passed to barplot
#'
#' @return A \code{data.frame} when tidy is TRUE otherwise a matrix or a plot
#'
#' @importFrom graphics barplot
#'
#' @export
get_frequency <- function(segments, normalize = FALSE, tidy = FALSE,
                          plot = FALSE, ...) {
    # calculate the frequency of each state in each cell
    freqs <- lapply(segments, function(x) {
        # tabularize states
        tab <- table(x$state)

        # normalize if true
        if (normalize) {
            tab <- tab / sum(tab)
        }

        # make a data.frame
        data.frame(
            state = names(tab),
            frequency = as.numeric(tab)
        )
    })

    # bind rows
    res <- do.call(rbind, freqs)

    # add cells columns and remove row names
    res$cell <- rep(names(segments), sapply(freqs, nrow))
    rownames(res) <- NULL

    # return
    if (tidy) {
        return(res)
    } else {
        mat <- matrix(res$frequency,
                      nrow = length(unique(res$state)),
                      dimnames = list(unique(res$state),
                                      unique(res$cell)))
        if (plot) return(barplot(mat, ...))
        else return(mat)
    }
}

#' Get the width of the segments in each cell type
#'
#' @inheritParams annotate_segments
#' @param average A logical. Whether the width should be averaged across cells.
#'
#' @return A \code{data.frame}
#'
#' @importFrom GenomicRanges width
#'
#' @export
get_width <- function(segments, average = FALSE) {
    # get width of each state in each cell
    wdth <- lapply(segments, function(x) {
        df <- data.frame(
            state = x$state,
            width = GenomicRanges::width(x)
        )

        # average if true
        if (average) {
            # calculate the mean and round up to nearest integer
            df <- aggregate(df$width,
                            by = list(state = df$state),
                            FUN = function(x) round(mean(x)))

            # rename columns
            names(df) <- c('state', 'width')
        }

        # return
        df
    })

    # bind rows
    res <- do.call(rbind, wdth)

    # add cells columns and remove row names
    res$cell <- rep(names(segments), sapply(wdth, nrow))
    rownames(res) <- NULL

    # return
    return(res)
}

#' Get the distances of the segments to the TSS in each cell type
#'
#' @inheritParams annotate_segments
#' @param average A logical. Whether the width should be averaged across cells.
#'
#' @return A \code{data.frame}
#'
#' @importFrom GenomicRanges mcols
#'
#' @export
get_distance <- function(segments, average = FALSE) {
    # get distance of each state in each cell
    dist <- lapply(segments, function(x) {
        # stop if not annotated properly
        stopifnot('distanceToTSS' %in% names(mcols(x)))

        # continue if present
        df <- data.frame(
            state = x$state,
            distance = x$distanceToTSS
        )

        # average if true
        if (average) {
            # calculate the mean and round up to nearst integer
            df <- aggregate(df$distanceToTSS,
                            by = list(state = df$state),
                            FUN = function(x) round(mean(x)))
        }

        # return
        df
    })

    # bind rows
    res <- do.call(rbind, dist)

    # add cells columns and remove row names
    res$cell <- rep(names(segments), sapply(dist, nrow))
    rownames(res) <- NULL

    # return
    return(res)
}

#' Extract binariezed data for a GRange
#'
#' @param bins A SummarizedExperiment object
#' @param which A GRanges object
#' @param tidy A logical
#'
#' @return A GRanges object
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors values subjectHits
#' @importFrom SummarizedExperiment assay rowRanges
#'
#' @export
extract_bins <- function(bins, which, tidy = FALSE) {
    # find overlaps between the selected range and row ranges of se
    ol <- findOverlaps(which, rowRanges(bins))

    # subset the se object
    hits <- bins[subjectHits(ol)]

    # extract ranges
    gr <- rowRanges(hits)

    # add assay data as metadata to ranges
    values(gr) <- cbind(values(gr), as.data.frame(assay(hits)))

    # tidy if true
    if (tidy) tidy_ranges(gr, columns = colnames(hits))
    else gr
}

#' Compare two or more models
#'
#' @param objs A list of segmentation items
#' @param type A string. What to compare. Default to 'emission'
#' @param plot A logical.
#' @param ... Other arguments passed to plot
#'
#' @return A numeric vector or a plot with the same values.
#'
#' @importFrom stats cor
#'
#' @export
compare_models <- function(objs, type = 'emission', plot = FALSE, ...) {
    if (type == 'emission') {
        # compare based on the emission parameters
        e <- lapply(objs, function(x) {
            emission(x)
        })

        c <- lapply(e, function(x) {
            m <- cor(t(e[[length(e)]]), t(x))
            apply(m, 1, max)
        })

        l <- apply(do.call(cbind, c), 2, max)
    } else if (type == 'likelihood') {
        # compare based on the likelihood
        l <- lapply(objs, function(x) {
            d <- model(x)
            d$likelihood
        })
        l <- unlist(l)
    }

    # plot or return vector
    if (plot) plot(l, ...)
    else return(l)
}

#' Visualize the model output
#'
#' @param obj A segmentation object
#' @param type A string. Which kind of parameter to print. Default is 'emission'
#' and possible values are 'emission', 'transition', 'overlap', 'TSS' or 'TES'
#' @param ... Other arguments to path to Heatmap
#'
#' @return A heatmap
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapList add_heatmap
#'
#' @export
plot_heatmap <- function(obj, type = 'emission', ...) {
    mat <- slot(obj, type)
    if (is.matrix(mat)) {
        Heatmap(mat, ...)
    } else if (is.list(mat)) {
        hh <- lapply(mat, function(x) Heatmap(x, ...))
        hml <- HeatmapList()
        for(i in hh) hml <- add_heatmap(hml, i)
        hml
    }
}
