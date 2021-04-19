#' segmentation objects
#'
#' The \code{segmentation} class consists of matrices and lists. The components
#' contain the output of the chromatin segmentation analysis. Loading the input
#' data is optional. The object is returned as a result of calling
#' \code{\link{learn_model}} or reading its already existing output.
#'
#' @slot model list. The \code{list} consists of 6 items corresponding
#' to the contents of the \code{model_#.txt} file. These are
#' \code{number_states} and \code{number_marks} for the numbers of states
#' and marks in the model; \code{likelihood} and \code{probinit} for the
#' likelihood and the initial probabilities of the multi-state model;
#' \code{transitionprobs} and \code{emissionprobs} for the probabilities
#' of the transitions and emissions parameters of the model. Can be
#' accessed using \code{\link{model}}.
#' @slot emission matrix. The \code{matrix} contains the emission
#' parameters of n states (rows) for n marks (columns) corresponding to
#' the contents of the \code{emission_#.txt} file. Can be accessed using
#' \code{\link{emission}}.
#' @slot transition matrix. The \code{matrix} contains the transition
#' parameters of n by n states corresponding to the contents of the
#' \code{transition_#.txt} file. Can be accessed using
#' \code{\link{transition}}.
#' @slot overlap list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{matrix} of the overlap enrichment of n states
#' (rows) at n genomic annotations (columns) corresponding to the contents
#' of the \code{<cell>_#_overlap.txt} files. Can be accessed using
#' \code{\link{overlap}}.
#' @slot TSS list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{matrix} of the overlap enrichment of n states
#' (rows) at n locations around the transcription start site (TSS)
#' (columns) corresponding to the contents of the
#' \code{<cell>_#_TSS_neighborhood.txt} files. Can be accessed using
#' \code{\link{TSS}}.
#' @slot TES list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{matrix} of the overlap enrichment of n states
#' (rows) at n locations around the transcription end site (TES)
#' (columns) corresponding to the contents of the
#' \code{<cell>_#_TES_neighborhood.txt} files. Can be accessed using
#' \code{\link{TES}}.
#' @slot segment list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{\link[GenomicRanges]{GRanges}} object containing the
#' segmentation and assigned states as a metadata column 'state'. These
#' contents correspond to the \code{<cell>_#_segment.bed} files. Annotations
#' of the ranges are optional. Can be accessed using \code{\link{segment}}.
#' @slot bins list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object containing the binarized input data. The coordinates of the bins
#' are saved as the \code{\link[SummarizedExperiment]{rowRanges}} each
#' assigned to a state and the binary data itself is saved as
#' \code{\link[SummarizedExperiment]{assay}}. Can be accessed using
#' \code{\link{bins}}.
#' @slot counts list. A \code{list} of n number of cells/conditions items.
#' Each item is a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object containing the read counts in bins. The coordinates of the bins
#' are saved as the \code{\link[SummarizedExperiment]{rowRanges}} each
#' assigned to a state and the counts data itself is saved as
#' \code{\link[SummarizedExperiment]{assay}}. Can be accessed using
#' \code{\link{counts}}.
#'
#' @name segmentation
#'
#' @aliases class:segmentation segmentation-class
#'
#' @seealso \code{\link{learn_model}}
#'
#' @export
setClass('segmentation',
         representation(model = 'list',
                        emission = 'matrix',
                        transition = 'matrix',
                        overlap = 'list',
                        TSS = 'list',
                        TES = 'list',
                        segment = 'list',
                        bins = 'list',
                        counts = 'list'))

#' Accessors for the \code{segmentation} objects
#'
#' These functions can be used to access the contents of \code{segmentation}
#' objects as well as modifying them.
#'
#' @param object An object of class \code{segmentation}
#' @param cell A string
#' @param ... Other argument passed to the accessors
#'
#' @return The data in the corresponding slot or a subset of it.
#'
#' @seealso segmentation
#'
#' @name accessors
NULL

#' @rdname accessors
#'
#' @examples
#' model(test_obj)
#'
#' @export
setGeneric('model', function(object) standardGeneric('model'))

#' @rdname accessors
setMethod('model', 'segmentation', function(object) object@model)

#' @rdname accessors
#'
#' @examples
#' emission(test_obj)
#'
#' @export
setGeneric('emission', function(object) standardGeneric('emission'))

#' @rdname accessors
setMethod('emission', 'segmentation', function(object) object@emission)

#' @rdname accessors
#'
#' @examples
#' transition(test_obj)
#'
#' @export
setGeneric('transition', function(object) standardGeneric('transition'))

#' @rdname accessors
setMethod('transition', 'segmentation', function(object) object@transition)

#' @rdname accessors
#'
#' @examples
#' overlap(test_obj)
#' overlap(test_obj, cell = 'K562')
#'
#' @export
setGeneric('overlap', function(object, ...) standardGeneric('overlap'))

#' @rdname accessors
setMethod('overlap',
          'segmentation',
          function(object, cell) {
              ol <- object@overlap
              if (!missing(cell)) {
                  stopifnot(cell %in% names(ol))
                  ol[cell]
              } else {
                  ol
              }
          })

#' @rdname accessors
#'
#' @examples
#' TSS(test_obj)
#' TSS(test_obj, cell = 'K562')
#'
#' @export
setGeneric('TSS', function(object, ...) standardGeneric('TSS'))

#' @rdname accessors
setMethod('TSS',
          'segmentation',
          function(object, cell) {
              tss <- object@TSS
              if (!missing(cell)) {
                  stopifnot(cell %in% names(tss))
                  tss[cell]
              } else {
                  tss
              }
          })

#' @rdname accessors
#'
#' @examples
#' TES(test_obj)
#' TES(test_obj, cell = 'K562')
#'
#' @export
setGeneric('TES', function(object, ...) standardGeneric('TES'))

#' @rdname accessors
setMethod('TES',
          'segmentation',
          function(object, cell) {
              tes <- object@TES
              if (!missing(cell)) {
                  stopifnot(cell %in% names(tes))
                  tes[cell]
              } else {
                  tes
              }
          })

#' @rdname accessors
#'
#' @examples
#' segment(test_obj)
#' segment(test_obj, cell = 'K562')
#'
#' @export
setGeneric('segment', function(object, ...) standardGeneric('segment'))

#' @rdname accessors
setMethod('segment',
          'segmentation',
          function(object, cell) {
              segs <- object@segment
              if (!missing(cell)) {
                  stopifnot(cell %in% names(segs))
                  segs[cell]
              } else {
                  segs
              }
          })

#' @rdname accessors
#'
#' @examples
#' bins(test_obj)
#'
#' @export
setGeneric('bins', function(object, ...) standardGeneric('bins'))

#' @rdname accessors
setMethod('bins',
          'segmentation',
          function(object, cell) {
              bins <- object@bins
              if (!missing(cell)) {
                  stopifnot(cell %in% names(bins))
                  bins[cell]
              } else {
                  bins
              }
          })

#' @rdname accessors
#'
#' @examples
#' counts(test_obj)
#'
#' @export
setGeneric('counts', function(object, ...) standardGeneric('counts'))

#' @rdname accessors
setMethod('counts',
          'segmentation',
          function(object, cell) {
              counts <- object@counts
              if (!missing(cell)) {
                  stopifnot(cell %in% names(counts))
                  counts[cell]
              } else {
                  counts
              }
          })

#' @rdname accessors
#'
#' @examples
#' likelihood(test_obj)
#'
#' @export
setGeneric('likelihood', function(object) standardGeneric('likelihood'))

#' @rdname accessors
setMethod('likelihood', 'segmentation', function(object) {
    object@model$likelihood
})

#' @rdname accessors
#'
#' @examples
#' cells(test_obj)
#'
#' @export
setGeneric('cells', function(object) standardGeneric('cells'))

#' @rdname accessors
setMethod('cells',
          'segmentation',
          function(object) {
              names(object@overlap)
          })

#' @rdname accessors
#'
#' @examples
#' states(test_obj)
#'
#' @export
setGeneric('states', function(object) standardGeneric('states'))

#' @rdname accessors
setMethod('states',
          'segmentation',
          function(object) {
              ss <- rownames(object@emission)
              if (is.null(ss)) {
                  as.character(
                      seq_len(nrow(object@emission))
                  )
              } else {
                  ss
              }
          })

#' @rdname accessors
#'
#' @examples
#' markers(test_obj)
#'
#' @export
setGeneric('markers', function(object) standardGeneric('markers'))

#' @rdname accessors
setMethod('markers',
          'segmentation',
          function(object) {
              colnames(object@emission)
          })

#' Methods to interact with \code{segmentation} objects
#'
#' These functions can be used to interact with \code{segmentation} objects for
#' purposes other than accessing or modifying their contents.
#'
#' @param object An object of class \code{segmentation}
#'
#' @return Prints a summary of the \code{segmentation} object contents.
#'
#' @seealso segmentation
#' @seealso accessors
#'
#' @name methods
NULL

#' @rdname methods
#'
#' @examples
#' show(test_obj)
#'
#' @export
setMethod(
    'show',
    signature = 'segmentation',
    definition = function(object) {
        # describe object class
        cat("# An object of class '", class(object), "' \n", sep = "")

        # describe information contents
        cat("# Contains a chromatin segmentation model:", "\n", sep = "")
        cat("## States: ", paste(states(object), collapse = ' '),"\n", sep = "")
        cat("## Marks: ", paste(markers(object), collapse = ' '),"\n", sep = "")
        cat("## Cells: ", paste(cells(object), collapse = ' '), "\n", sep = "")

        # describe slots and accessors
        cat("# Contains nine slots: ", "\n", sep = "")
        cat("## model: use 'model(object)' to access", "\n", sep = "")
        cat("## emission: use 'emission(object)' to access", "\n", sep = "")
        cat("## transition: use 'transition(object)' to access", "\n", sep = "")
        cat("## overlap: use 'overlap(object)' to access", "\n", sep = "")
        cat("## TSS: use 'TSS(object)' to access", "\n", sep = "")
        cat("## TES: use 'TES(object)' to access", "\n", sep = "")
        cat("## segment: use 'segment(object)' to access", "\n", sep = "")
        cat("## bins: use 'bins(object)' to access", "\n", sep = "")
        cat("## counts: use 'counts(object)' to access", "\n", sep = "")

        # more info
        cat("# For more info about how to use the object, use ?accessors")

        # return NULL
        invisible(NULL)
    }
)

# TODO: implement update method
# TODO: implement validate method
