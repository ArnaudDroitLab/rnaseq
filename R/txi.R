#' Filter out txi object
#'
#' This function allow to filter the txi object by taking a list of sample to
#' keep. All the expected matrices will be filtered (counts, abundance, fpkm,
#' ruvg_counts and combat_counts.
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param samples The list of samples to keep.
#'
#' @return The filter txi object.
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#'
#' @examples
#' txi <- get_demo_txi() # sample names are "a", "b", "c" and "d"
#' txi_filter <- filter_txi(txi, c("a", "c")) # Only keep "a" and "c"
#'
#' @export
filter_txi <- function(txi, samples) {
    validate_txi(txi)
    filter_matrices <- function(txi, name) {
        stopifnot(all(samples %in% colnames(txi[[name]])))
        txi[[name]] <- txi[[name]][,colnames(txi[[name]]) %in% samples]
        txi
    }
    txi <- filter_matrices(txi, "counts")
    txi <- filter_matrices(txi, "abundance")
    txi <- filter_matrices(txi, "length")
    if (!is.null(txi$fpkm)) {
        txi <- filter_matrices(txi, "fpkm")
    }
    if (!is.null(txi$ruvg_counts)) {
        txi <- filter_matrices(txi, "ruvg_counts")
    }
    if (!is.null(txi$combat_counts)) {
        txi <- filter_matrices(txi, "combat_counts")
    }
    if (!is.null(txi$extra_count_matrix)) {
        txi <- filter_matrices(txi, "extra_count_matrix")
    }
    validate_txi(txi)
    txi
}

#' Create a dummy txi objects
#'
#' This function will create a dummy txi object that will have placeholders for
#' some of the tables.
#'
#' This function is useful if you only have access to raw counts or abundance
#' tables. Otherwise, it is recommended to use the \code{import_kallisto}
#' function as it will create a complete and valid \code{txi} object.
#'
#' The \code{txi} will contains extra elements to specify how it was created
#' and will provide extra checks to avoid using it incorrectly with the other
#' functions. For instance, to use the \code{produce_pca_df} function, the
#' \code{txi} must have a valid \code{abundance} table. If the \code{dummy_txi}
#' was created without such a table, it will return an error message to avoid
#' producing an invalid result.
#'
#' You have to provide at least the raw counts or the abundance table **and**
#' the corresponding annotation in the same format as the one produced by the
#' \code{anno:prepare_anno} function.
#'
#' @param tables A \code{list} of the tables to add to the \code{dummy_txi}
#' object. Must contains at \code{counts} (raw counts) or \code{abundance}
#' (usually TPM). The \code{anno} is also mandatory.
#' @param txOut Are the counts corresponding to transcripts abundance? Default:
#' \code{FALSE}.
#'
#' @return The \code{dummy_txi} object.
#'
#' @examples
#' txi <- get_demo_txi()
#'
#' lst <- list()
#' lst$counts <- txi$counts
#' lst$abundance <- txi$abundance
#' lst$anno <- txi$anno
#'
#' # We then create a dummy txi object
#' dummy_txi <- create_dummy_txi(txi)
#' # This dummy_txi contains a placeholder for the length matrix
#'
#' @export
create_dummy_txi <- function(tables, txOut = FALSE) {
    stopifnot(is(tables, "list"))
    stopifnot(any(c("counts", "abundance") %in% names(tables)))
    stopifnot("anno" %in% names(tables))
    stopifnot(is(tables$anno, "data.frame"))
    if (!is.null(tables$abundance)) {
        stopifnot(is(tables$abundance, "matrix"))
        m <- tables$abundance
    }
    if (!is.null(tables$counts)) {
        stopifnot(is(tables$counts, "matrix"))
        m <- tables$counts
    }
    stopifnot(is(txOut, "logical"))

    # Provided tables shoud not be identical
    if (!is.null(tables$counts) & !is.null(tables$abundance)) {
        stopifnot(!identical(tables$counts, tables$abundance))
    }
    if (!is.null(tables$counts) & !is.null(tables$length)) {
        stopifnot(!identical(tables$counts, tables$length))
    }
    if (!is.null(tables$abundance) & !is.null(tables$length)) {
        stopifnot(!identical(tables$abundance, tables$length))
    }

    txi <- tables
    add_matrix <- function(txi, m, n) {
        if (!n %in% names(txi)) {
            txi[[n]] <- m
        }
        txi
    }
    txi <- add_matrix(txi, m, "counts")
    txi <- add_matrix(txi, m, "abundance")
    txi <- add_matrix(txi, m, "length")
    txi$txOut <- txOut
    txi$dummy <- names(tables)
    validate_txi(txi)
    txi
}

#' Validate the txi object content
#'
#' @param txi a \code{txi} object produced by \code{import_kallisto}.
#'
#' @return Invisibly returns \code{TRUE} if all checks pass.
#'
#' @examples
#' txi <- get_demo_txi()
#' validate_txi(txi) 
#'
#' @export
validate_txi <- function(txi) {
    # Global
    stopifnot(is(txi, "list"))
    stopifnot(all(c("counts", "abundance", "length", "anno") %in% names(txi)))

    # Matrices
    stopifnot(is(txi$counts, "matrix"))
    stopifnot(is(txi$abundance, "matrix"))
    stopifnot(is(txi$length, "matrix"))
    stopifnot(is.numeric(txi$counts))
    stopifnot(is.numeric(txi$abundance))
    stopifnot(is.numeric(txi$length))
    stopifnot(identical(colnames(txi$counts), colnames(txi$abundance)))
    stopifnot(identical(colnames(txi$counts), colnames(txi$length)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$abundance)))
    stopifnot(identical(rownames(txi$counts), rownames(txi$length)))

    if (!is.null(txi$ruvg_counts)) {
        stopifnot(is(txi$ruvg_counts, "matrix"))
        stopifnot(is.numeric(txi$ruvg_counts))
        stopifnot(identical(colnames(txi$counts), colnames(txi$ruvg_counts)))
        stopifnot(identical(rownames(txi$counts), rownames(txi$ruvg_counts)))
    }

    if (!is.null(txi$combat_counts)) {
        stopifnot(is(txi$combat_counts, "matrix"))
        stopifnot(is.numeric(txi$combat_counts))
        stopifnot(identical(colnames(txi$counts), colnames(txi$combat_counts)))
        stopifnot(identical(rownames(txi$counts), rownames(txi$combat_counts)))
    }

    # Data.frame
    stopifnot(is(txi$anno, "data.frame"))
    expected_col <- c("id", "ensembl_gene", "symbol", "entrez_id", "transcript_type")
    stopifnot(expected_col %in% colnames(txi$anno))
    stopifnot(identical(txi$anno$id, rownames(txi$counts)))

    # txOut
    stopifnot(is(txi$txOut, "logical"))
    anno_not_ercc <- dplyr::filter(txi$anno, !stringr::str_detect(id, "ERCC"))
    if (txi$txOut) {
        stopifnot(all(anno_not_ercc$id != anno_not_ercc$ensembl_gene))
    } else {
        stopifnot(all(anno_not_ercc$id == anno_not_ercc$ensembl_gene))
    }

    invisible(TRUE)
}

produce_txi <- function(files, anno, ignoreTxVersion = TRUE, use_ruv = FALSE) {
    # TODO: replace use_ruv by normalize %in% c("ruvg", "combat", "both")
    txi <- list()
    txi$tx <- import_kallisto(files, anno = anno, txOut = TRUE, ignoreTxVersion = ignoreTxVersion)
    txi$gene <- summarize_to_gene(txi$tx, anno = anno, ignoreTxVersion = ignoreTxVersion)

    # RUV
    if (use_ruv) {
        txi$tx <- ruvg_normalization(txi$tx,
                                     housekeeping_genes = housekeeping_genes,
                                     ignoreTxVersion = ignoreTxVersion)
        txi$gene <- ruvg_normalization(txi$gene,
                                       housekeeping_genes = housekeeping_genes,
                                       ignoreTxVersion = ignoreTxVersion)
    }
    txi
}

