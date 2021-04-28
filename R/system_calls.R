#' Call Java \code{LearnModel}
#'
#' Call the Java module \code{LearnModel} which learns a multi-state model
#' from ChIP-seq data.
#'
#' @param inputdir A string. The path to binarized files.
#' @param outputdir A string. The path to a directory where output will be
#' written.
#' @param numstates An integer. The number of desired states in the model.
#' @param coordsdir A string. The path to genomic coordiantes files.
#' @param anchorsdir A string. The path to the genomic anchors files.
#' @param chromsizefile A string. The path to the chromosomes sizes file.
#' @param assembly A string. The name of the genomic assembely.
#' @param optional A string. Other optional arguments passed to the Java
#' command.
#'
#' @return \code{NULL}. Output files are written to the output directory.
#'
#' @seealso learn_model
.LearnModel <- function(inputdir, outputdir, numstates, coordsdir, anchorsdir,
                       chromsizefile, assembly, optional) {
    # locate LearnModel module
    LearnModel <- paste(
        system.file("java",
                    "ChromHMM.jar",
                    package = "segmenter"),
        "LearnModel"
    )

    # optional passed as arguments
    genome_annotation <- paste(
        '-u', coordsdir,
        '-v', anchorsdir,
        '-l', chromsizefile
    )

    # make options string
    if (missing(optional)) {
        optional <- paste("-noautoopen", # change default
                          "-nobrowser",  # change default
                          "-noimage",    # change default
                          collapse = '')
    }

    # run cmd
    system2(
        command = "java",
        args = paste(
            '-jar',
            LearnModel,
            genome_annotation,
            optional,
            inputdir, outputdir, numstates, assembly),
        stdout = paste(outputdir, 'log.txt', sep = '/'),
        stderr = paste(outputdir, 'log_error.txt', sep = '/')
    )

    # return
    return(invisible(NULL))
}

#' Call Java \code{BinarizeBed}
#'
#' Call the Java module \code{BinarizeBed} which binarize a bed file of the
#' aligned reads.
#'
#' @param inputdir A string. The path to bed files.
#' @param outputdir A string. The path to a directory where output will be
#' written.
#' @param chromsizefile A string. The path to the chromosomes sizes file.
#' @param cellmarkfiletable A tab delimited files of three columns.
#' The columns contains the cell, mark and the name or the bed file.
#' @param binsize An integer. The bin size to use. Default is 200.
#' @param type A string. The file type 'bam' or 'bed'.
#'
#' @return \code{NULL}. Output files are written to the output directory.
#'
#' @seealso binarize_bed
.Binarize <- function(inputdir, cellmarkfiletable, chromsizefile, binsize,
                     outputdir, type) {
    # locate module
    type <- switch (type,
                    'bam' = 'BinarizeBam',
                    'bed' = 'BinarizeBed')
    module <- paste(
        system.file("java",
                    "ChromHMM.jar",
                    package = "segmenter"),
        type
    )

    # make options string
    optional <- paste("-b", binsize)

    # run cmd
    system2(
        command = "java",
        args = paste(
            '-jar',
            module,
            optional,
            chromsizefile, inputdir, cellmarkfiletable, outputdir),
        stdout = paste(outputdir, 'log.txt', sep = '/'),
        stderr = paste(outputdir, 'log_error.txt', sep = '/')
    )

    # return
    return(invisible(NULL))
}
