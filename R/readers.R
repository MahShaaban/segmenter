#' Read \code{chromsizefile}
#'
#' The file should contain exactly two columns. One for the name of the
#' chromosome and the other for its length.
#'
#' @param file A string. The path to the file.
#'
#' @return A \code{data.frame}
#'
#' @examples
#' # locate the file
#' chromsizefile <- system.file('extdata/CHROMSIZES',
#'                              'hg18.txt',
#'                              package = 'chromhmmData')
#'
#' # read the file
#' read_chromsize_file(chromsizefile)
#'
#' @importFrom utils count.fields read.delim
#'
#' @export
read_chromsize_file <- function(file) {
    # stop if the number of columns is not two
    stopifnot(max(count.fields(file, sep = '\t')) == 2)

    # otherwise, read the file into a data.frame
    df <- read.delim(file,
                     col.names = c('seqname', 'width'),
                     header = FALSE)

    # return
    return(df)
}

#' Read \code{segments} files
#'
#' The segments files are the output of running \code{learn_model} and named
#' \code{<cell>_3_segment.bed}
#'
#' @param file A string. The path to the file.
#' @param states A \code{character} vector. The names of the states.
#'
#' @return A \code{data.frame}
#'
#' @examples
#' # locate the file
#' segmentfile <- file.path(tempdir(), 'GM12878_3_segments.bed')
#'
#' # read the file
#' segs <- read_segements_file(segmentfile)
#' head(segs)
#'
#' @importFrom utils count.fields read.delim
#'
#' @export
read_segements_file <- function(file, states) {
    # stop if the number of columns is not two
    stopifnot(max(count.fields(file)) == 4)

    # otherwise, read the file into a data.frame
    df <- read.delim(file,
                     col.names = c('seqname', 'start', 'end', 'state'),
                     header = FALSE)

    # replace state names when provided
    if (!missing(states)) {
        # stop if the number of unique states doesn't match
        stopifnot(length(states) == length(unique(df$state)))

        # otherwise, match ordered values of the two vectors
        ind <- as.numeric(gsub("[^\\d]+", "", df$state, perl = TRUE))
        df$state <- states[ind]
    }

    # return
    return(df)
}

#' Read \code{modelfile}
#'
#' The model file is the output of running \code{learn_model} and named
#' \code{model_#.txt}
#'
#' @param file A string. The path to the file.
#'
#' @return A \code{data.frame}
#'
#' @examples
#' # locate the file
#' modelfile <- file.path(tempdir(), 'model_3.txt')
#'
#' # read the file
#' read_model_file(modelfile)
#'
#' @export
read_model_file <- function(file) {
    # read lines from file
    lns <- readLines(file)

    # extract info in the first line
    first_line <- unlist(strsplit(lns[1], '\t'))

    res <- list(
        number_states = as.integer(first_line[1]),
        number_marks = as.integer(first_line[2]),
        likelihood = as.numeric(first_line[4])
    )

    # index the first word in the rest of the lines
    ind <- sub('\t.*', '', lns[-1])

    # loop over different words and read
    for (i in unique(ind)) {
        # split the lines by \t and simplify into a matrix
        mat <- strsplit(lns[-1][ind == i], '\t')
        mat <- do.call(rbind, mat)

        # subset and transform into a data.frame
        df <- as.data.frame(mat[, -1])

        # add the data.frame to the return object
        res[[i]] <- df
    }

    # TODO: add column names
    # TODO: ensure correct types of columns

    # return
    return(res)
}

#' Read \code{emissions} file
#'
#' The segments files are the output of running \code{learn_model} and named
#' \code{emissions_3_segment.bed}
#'
#' @param file A string. The path to the file.
#' @param states A \code{character} vector. The names of the states.
#' @param marks A \code{character} vector. The names of the marks
#'
#' @return A \code{matrix}
#'
#' @examples
#' # locate the file
#' fl <- file.path(tempdir(), 'emissions_3.txt')
#'
#' # read the file
#' read_emissions_file(fl)
#'
#' @importFrom utils read.delim
#'
#' @export
read_emissions_file <- function(file, states, marks) {
    # read the file into a data.frame
    df <- read.delim(file)

    # format into a matrix
    mat <- as.matrix(df[, -1])

    if (!missing(marks))    colnames(mat) <- marks
    if (!missing(states))   rownames(mat) <- states

    # return
    return(mat)
}

#' Read \code{transitions} file
#'
#' The segments files are the output of running \code{learn_model} and named
#' \code{transitions_3_segment.bed}
#'
#' @param file A string. The path to the file.
#' @param states A \code{character} vector. The names of the states.
#'
#' @return A \code{matrix}
#'
#' @examples
#' # locate the file
#' fl <- file.path(tempdir(), 'transitions_3.txt')
#' 
#' # read the file
#' read_transitions_file(fl)
#'
#' @importFrom utils read.delim
#'
#' @export
read_transitions_file <- function(file, states) {
    # read the file into a data.frame
    df <- read.delim(file)

    # format into a matrix
    mat <- as.matrix(df[, -1])

    if (!missing(states))   rownames(mat) <- states
    if (!missing(states))   colnames(mat) <- states

    # return
    return(mat)
}

#' Read \code{segments} files
#'
#' The segments files are the output of running \code{learn_model} and named
#' \code{<cell>_3_overlap.txt}
#'
#' @param file A string. The path to the file.
#' @param states A \code{character} vector. The names of the states.
#' @param regions A \code{character} vector. The names of the regions.
#'
#' @return A \code{matrix}
#'
#' @examples
#' # locate the file
#' fl <- file.path(tempdir(), 'GM12878_3_overlap.txt')
#' 
#' # read the file
#' read_overlap_file(fl)
#'
#' @importFrom utils read.delim
#'
#' @export
read_overlap_file <- function(file, states, regions) {
    # read the file into a data.frame
    df <- read.delim(file)

    # format into a matrix
    mat <- as.matrix(df[-nrow(df), -1])

    if (!missing(regions))  colnames(mat) <- regions
    if (!missing(states))   rownames(mat) <- states

    # return
    return(mat)
}

#' Read \code{enrichment} files
#'
#' The segments files are the output of running \code{learn_model} and named
#' \code{<cell>_3_TSS.txt} or \code{<cell>_3_TES.txt}.
#'
#' @param file A string. The path to the file.
#' @param states A \code{character} vector. The names of the states.
#' @param regions A \code{character} vector. The names of the regions.
#'
#' @return A \code{matrix}
#'
#' @examples
#' # locate the file
#' fl <- file.path(tempdir(), 'GM12878_3_RefSeqTSS_neighborhood.txt')
#'
#' # read the file
#' read_enrichment_file(fl)
#'
#' @importFrom utils read.delim
#'
#' @export
read_enrichment_file <- function(file, states, regions) {
    # read the file into a data.frame
    df <- read.delim(file)

    # format into a matrix
    mat <- as.matrix(df[, -1])

    if (!missing(regions))  colnames(mat) <- regions
    if (!missing(states))   rownames(mat) <- states

    # return
    return(mat)
}

#' Read \code{bins} files
#'
#' The files contain the cell and the chromosome info in the first line and
#' the binarized data from all marks in the rest.
#'
#' @param file A string. The path to the file.
#'
#' @return A \code{list} of 3 items: cell, seqname and binaries.
#'
#' @examples
#' # locate the file
#' fl <- system.file('extdata/SAMPLEDATA_HG18/',
#'                   'GM12878_chr11_binary.txt.gz',
#'                   package = 'segmenter')
#'
#' # read the file
#' read_bins_file(fl)
#'
#' @importFrom utils read.delim
#'
#' @export
read_bins_file <- function(file) {
    # read the first line
    first_line <- readLines(file, n = 1)

    # split line by \t and unlist
    first_line <- unlist(strsplit(first_line, '\t'))

    # read the rest of the file into a data.frame
    binaries <- read.delim(file, skip = 1)

    # make a return object as lines
    res <- list(cell = first_line[1],
                seqname = first_line[2],
                binaries = binaries)

    # return
    return(res)
}

#' Read \code{cellmarktable} file
#'
#' The file should contain at least three columns: cell, mark and file for the
#' names of the cells/conditions, the available marks and binarized data files.
#'
#' @param file A string. The path to the file.
#'
#' @return A \code{data.frame}
#'
#' @examples
#' # locate the file
#' fl <- system.file('extdata',
#'                   'cell_mark_table.tsv',
#'                   package = 'segmenter')
#'
#' # read the file
#' read_cellmark_file(fl)
#'
#' @importFrom utils count.fields read.delim
#'
#' @export
read_cellmark_file <- function(file) {
    # stop if the number of columns is less than 3
    stopifnot(max(count.fields(file, sep = '\t')) >= 3)

    # read file into a data.frame
    df <- read.delim(file, header = FALSE)

    # make a column names vector
    col.names <- c('cell', 'mark', 'file', 'control_file')

    # rename columns
    names(df) <- col.names[seq_len(ncol(df))]

    # return
    return(df)
}

#' Read \code{bam} files
#'
#' Count the reads in each range of the \code{GRanges} object
#'
#' @param file A string. The path to the file.
#' @param features A \code{GRanges} object.
#' @param ... Other arguments passed to \code{\link[bamsignals]{bamCount}}.
#'
#' @return A \code{matrix}
#'
#' @examples
#' # locate the bam file
#' bam_file <- system.file("extdata", "randomBam.bam", package = "bamsignals")
#'
#' # load a granges object
#' rand_anno <- system.file("extdata",
#'                          "randomAnnot.Rdata",
#'                          package = "bamsignals")
#' features <- GenomicRanges::promoters(get(load(rand_anno)))
#'
#' # count reads in ranges
#' read_bam_file(bam_file, features)
#'
#' @importFrom bamsignals bamCount
#'
#' @export
read_bam_file <- function(file, features, ...) {
    # stop if no index file is in the same path as file
    stopifnot(file.exists(gsub(".bam$", ".bam.bai", file)))

    # count reads in features in file
    reads <- bamCount(file, features, ...)

    # return
    return(reads)
}
