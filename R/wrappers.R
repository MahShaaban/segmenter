#' Lean a multi-state model from chromatin data
#'
#' Integrate multiple ChIP-seq chromatin datasets of histone modifications,
#' transcription factors or other DNA binding proteins to build a multi-state
#' model of the combinatorial and spatial frequently occuring patterns.
#' The function uses as an input binarized ChIP-seq data and the genome
#' annotations on which the states will be discovered.
#'
#' @param inputdir A string. The path to binarized files.
#' @param outputdir A string. The path to a directory where output will be
#' written.
#' @param numstates An integer. The number of desired states in the model.
#' @param coordsdir A string. The path to genomic coordiantes files.
#' @param anchorsdir A string. The path to the genomic anchors files.
#' @param chromsizefile A string. The path to the chromosomes sizes file.
#' @param cells A \code{character} vector. The names of the cells as they occur
#' in the binarized files (first line).
#' @param annotation A string. The name of the type of annotation as it occurs
#' in the genomic annotation files.
#' @param binsize An integer. The number in bp used to generate binarized files.
#' @param inputbamdir A string. The path to the input bam files. Only used when
#' \code{count = TRUE}.
#' @param cellmarkfiletable A string. The path to the input files table. Only
#' used when \code{bins = TRUE}.
#' @param read_only A logical. Default is \code{FALSE}. Whether to look for and
#' load output files or generate the model from scratch.
#' @param read_bins A logical. Default is \code{FALSE}. Whether to load the
#' binarized data into the output object.
#' @param counts A logical. Default is \code{FALSE}. Whether to load the
#' reads counts in bins data into the output object.
#' @param assembly A string. The name of the genomic assembely.
#'
#' @return An object of class \code{\link{segmentation}} (see for details)
#' and the files written to the output directory.
#'
#' @details By default, this functions runs the analysis commands, writes the
#' output to files and loads it into an object of class
#' \code{\link{segmentation}}. In addition, the binarized data and the reads
#' counts in the bins can be loaded. When \code{read_only} is \code{TRUE}.
#' The functions looks for previously generated files in the \code{output}
#' directory and load them without rerunning the commands.
#'
#' @examples
#' # locate input and output files
#' inputdir <- system.file('extdata/ChromHMM/SAMPLEDATA_HG18',
#'                         package = 'segmenter')
#' outputdir <- tempdir()
#' coordsdir <- system.file('extdata/ChromHMM/COORDS',
#'                          package = 'segmenter')
#' anchorsdir <- system.file('extdata/ChromHMM/ANCHORFILES',
#'                           package = 'segmenter')
#' chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
#'                              'hg18.txt',
#'                              package = 'segmenter')
#'
#' # run command
#' obj <- learn_model(inputdir = inputdir,
#'                    outputdir = outputdir,
#'                    coordsdir = coordsdir,
#'                    anchorsdir = anchorsdir,
#'                    chromsizefile = chromsizefile,
#'                    numstates = 3,
#'                    assembly = 'hg18',
#'                    cells = c('K562', 'GM12878'),
#'                    annotation = 'RefSeq',
#'                    binsize = 200)
#'
#' # show the output
#' obj
#'
#' @seealso LearnModel
#'
#' @importFrom GenomicRanges granges makeGRangesFromDataFrame
#'
#' @export
learn_model <- function(inputdir, outputdir, numstates, coordsdir, anchorsdir,
                        chromsizefile, assembly, cells, annotation, binsize,
                        inputbamdir, cellmarkfiletable,
                        read_only = FALSE, read_bins = FALSE, counts = FALSE) {

    # make a temporary directory when missing
    if (missing(outputdir)) {
        if (read_only) stop("When reading only outputdir should be provided.")
        outputdir <- tempdir()
    }

    # call LearnModel from R
    if (!read_only) {
        LearnModel(inputdir,
                   outputdir,
                   numstates,
                   coordsdir,
                   anchorsdir,
                   chromsizefile,
                   assembly)
    }

    # capture output
    res <- list()

    ## model file
    res$model <- tryCatch(
        read_model_file(
            file.path(outputdir,
                      model_file(numstates))
        )
    )

    # emissions file
    res$emission <- tryCatch(
        read_emissions_file(
            file.path(outputdir,
                      emissions_file(numstates))
        )
    )

    # transitions file
    res$transition <- tryCatch(
        read_transitions_file(
            file.path(outputdir,
                      transitions_file(numstates))
        )
    )

    # overlap files
    res$overlap <- tryCatch({
        fls <- file.path(outputdir, overlap_files(numstates, cells))
        ol <- lapply(fls, read_overlap_file)
        names(ol) <- cells
        ol
    })

    # segment files
    res$segment <- tryCatch({
        fls <- file.path(outputdir, segments_files(numstates, cells))
        segs <- lapply(fls, function(x) {
            df <- read_segements_file(x)
            df <- df[df$start < df$end, ]
            makeGRangesFromDataFrame(df,
                                     keep.extra.columns = TRUE)
        })
        names(segs) <- cells
        segs
    })

    # TSS files
    res$TSS <- tryCatch({
        fls <- file.path(outputdir,
                         enrichment_files(numstates, cells, annotation, 'TSS'))
        tss <- lapply(fls, read_enrichment_file)
        names(tss) <- cells
        tss
    })

    # TES files
    res$TES <- tryCatch({
        fls <- file.path(outputdir,
                         enrichment_files(numstates, cells, annotation, 'TES'))
        tes <- lapply(fls, read_enrichment_file)
        names(tes) <- cells
        tes
    })

    # read bins
    if (read_bins) {
        # read bins from input director
        fls <- list.files(inputdir, pattern = 'binary.txt', full.names = TRUE)
        bins <- lapply(fls, read_bins_file)

        # index bins by cell
        ind <- unlist(lapply(bins, function(x) x$cell))

        # stop if names of cells don't match
        stopifnot(all(cells %in% unique(ind)))

        # make an empty list
        res$bins <- list()

        # loop over cells to read bins and put them in SummarizedExperiment
        for (i in cells) {
            stopifnot(i %in% names(res$segment))

            # make bins object
            se <- range_bins(bins[ind == i],
                             chromsizefile,
                             binsize = binsize,
                             return = 'SummarizedExperiment')

            # merge bins and segments
            res$bins[[i]] <- merge_segments_bins(res$segment[[i]], se)
        }
    } else {
        # otherwise, return an empty list
        res$bins <- list()
    }

    # read counts when true
    if (counts) {
        stopifnot(read_only)
        stopifnot(!missing(cellmarkfiletable))
        stopifnot(!missing(inputbamdir))

        # read counts in bins
        res$counts <- count_reads_ranges(
            lapply(res$bins, granges),
            cellmarkfiletable = cellmarkfiletable,
            inputbamdir = inputbamdir
        )
    } else {
        # otherwise, return an empty list
        res$counts <- list()
    }

    # make class
    obj <- new(
        'segmentation',
        model = res$model,
        emission = res$emission,
        transition = res$transition,
        overlap = res$overlap,
        segment = res$segment,
        TSS = res$TSS,
        TES = res$TES,
        bins = res$bins,
        counts = res$counts
    )

    # return
    return(obj)
}

#' Binarize the bam files
#'
#' Transform the aligned reads into a binary format.
#'
#' @param inputdir A string. The dirctory of the bam files.
#' @inheritParams learn_model
#' @param cellmarkfiletable A string. The path to the input files table. Only
#'
#' @return NULL. Write files to the outputdir
#'
#' @examples
#' # locate input and output files
#' inputdir <- system.file("extdata", package = "bamsignals")
#' cellmarkfiletable <- system.file('extdata/input','cell_mark_table.tsv', package = 'segmenter')
#' chromsizefile <- system.file('extdata/ChromHMM/CHROMSIZES',
#'                              'hg18.txt',
#'                               package = 'segmenter')
#' outputdir <- tempdir()
#'
#' # run command
#' binarize_bam(inputdir,
#'              chromsizefile = chromsizefile,
#'              cellmarkfiletable = cellmarkfiletable,
#'              outputdir = outputdir)
#'
#' # show output files
#' list.files(outputdir, pattern = '*_binary.txt')
#'
#' @seealso Binarize binarize_bed
#'
#' @export
binarize_bam <- function(inputdir, cellmarkfiletable, chromsizefile, binsize = 200,
                         outputdir) {

    # make a temporary directory when missing
    if (missing(outputdir)) {
        stop("outputdir should be specified.")
    }

    # call Binarize from R
    Binarize(inputdir,
             cellmarkfiletable,
             chromsizefile,
             binsize,
             outputdir,
             type = 'bam')
}

#' Binarize the bed files
#'
#' Transform the aligned reads into a binary format.
#'
#' @inheritParams binarize_bam
#'
#' @return NULL. Write files to the outputdir
#'
#' @seealso Binarize binarize_bam
#'
#' @export
binarize_bed <- function(inputdir, cellmarkfiletable, chromsizefile, binsize = 200,
                         outputdir) {

    # make a temporary directory when missing
    if (missing(outputdir)) {
        stop("outputdir should be specified.")
    }

    # call Binarize from R
    Binarize(inputdir,
             cellmarkfiletable,
             chromsizefile,
             binsize,
             outputdir,
             type = 'bed')
}

#' Make model file name
#'
#' @param numstates An integer
#'
#' @return A string
#'
#' @export
model_file <- function(numstates) {
    paste0('model_', numstates, '.txt')
}

#' Make emissions file name
#'
#' @param numstates An integer
#'
#' @return A string
#'
#' @export
emissions_file <- function(numstates) {
    paste0('emissions_', numstates, '.txt')
}

#' Make transitions file name
#'
#' @param numstates An integer
#'
#' @return A string
#'
#' @export
transitions_file <- function(numstates) {
    paste0('transitions_', numstates, '.txt')
}

#' Make overlap file names
#'
#' @param numstates An integer
#' @param cells A character vector
#'
#' @return A character vector
#'
#' @export
overlap_files <- function(numstates, cells) {
    paste(cells, numstates, 'overlap.txt', sep = '_')
}

#' Make segments file names
#'
#' @param numstates An integer
#' @param cells A character vector
#'
#' @return A character vector
#'
#' @export
segments_files <- function(numstates, cells) {
    paste(cells, numstates, 'segments.bed', sep = '_')
}

#' Make enrichment file names
#'
#' @param numstates An integer
#' @param cells A character vector
#' @param table A string
#' @param annotation A string
#'
#' @return A character vector
#'
#' @export
enrichment_files <- function(numstates, cells, table = 'RefSeq',
                             annotation = 'TSS') {
    paste(cells, numstates,
          paste0(table, annotation),
          'neighborhood.txt',
          sep = '_')
}
