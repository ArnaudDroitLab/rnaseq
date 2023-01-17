#' Save the analysis results in odysse format
#'
#' The goal of this function is to produce de deliverables in the Odysse
#' format. The characteristics of this format are:
#'     
#' * counts: A single count file for raw counts, TPM and FPKM. The anno is
#' included in the file and each count cell contains the 3 types of counts in
#' the <raw_counts>/<tpm>/<fpkm> format.
#' * de: There is a file produced for each comparison. The results contains the
#' annotations and the extra analysis (mean TPM counts and splicing analysis).
#' There is a distinction between the fold change and the ratio where the ratio
#' is 2^log2FoldChange and the fold change is symetrical (i.e.: ratios smaller
#' than 1 are converted to their negative value).
#'
#' In both cases, the results at the transcripts levels are found in the first
#' part of the file and the results at the gene levels are found in the second
#' part of the file.
#'
#' @param res The results obtained from the \code{produce_deliverables}
#' function or in a similar format.
#' @param outdir The directory where the tables will be saved
#' @param prefix The value to add at the beginning of each lines
#' @param use_ruv Use RUVg normalization? Needs to be pre-computed using the
#'                \code{ruvg_normalization} function. Default: \code{FALSE}.
#' @param digits Integer indicating the number of decimal places
#' @param force Should the files be re-created if they already exists? Default:
#' \code{FALSE}.
#' @param ... Extra params for the format_de_results function.
#' @param ncores Number of cores to use for de analysis. Default \code{1}.
#'
#' @return A vector of kallisto abundance filenames
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#'
produce_odysse_format <- function(res, outdir, prefix, use_ruv = FALSE,
                                  digits = 4, force = FALSE, ..., ncores = 1) {
    # TODO: add de_infos to validate all the DE are present? YES!!1
    # TODO: verbose
    stopifnot(is(res, "list"))
    stopifnot(c("txi", "counts", "de") %in% names(res))
    stopifnot(c("gene", "tx") %in% names(res$txi))
    stopifnot(c("gene", "tx") %in% names(res$counts))
    stopifnot(c("gene", "tx") %in% names(res$de))

    stopifnot(identical(res$infos_tables$de_infos$id_de, names(res$de$gene)))
    stopifnot(identical(res$infos_tables$de_infos$id_de, names(res$de$gene)))

    stopifnot(is.character(outdir))
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
    stopifnot(is.numeric(ncores))
    stopifnot(identical(round(ncores), ncores))
    stopifnot(ncores > 0)

    # Counts
    output_counts <- paste0(outdir, "/", prefix, "_counts.csv")
    if (!file.exists(output_counts) | force) {
        counts <- list()
        if (ncores == 1) {
            counts$tx <- format_counts(res$txi$tx, digits = digits,
                                       use_ruv = use_ruv)
            counts$gene <- format_counts(res$txi$gene, digits = digits,
                                         use_ruv = use_ruv)
        } else {
            counts <- parallel::mclapply(res$txi, format_counts, use_ruv,
                                         mc.cores = ncores)
        }
        counts <- rbind(counts$tx, counts$gene)

        write_csv(counts, output_counts)
    }

    # DE
    format_de_res <- function(current_de, current_id_de, lvl) {
        current_info_de <- dplyr::filter(res$infos_tables$de_infos,
                                  id_de == current_id_de)
        current_contrast_1 <- current_info_de$contrast_1
        current_contrast_2 <- current_info_de$contrast_2
        current_group <- current_info_de$group
        design <- res$infos_tables$design_infos
        i <- design[[current_group]] == current_contrast_1
        j <- design[[current_group]] == current_contrast_2
        current_samples_1 <- design$sample[i]
        current_samples_2 <- design$sample[j]

        format_de_odysse(current_de, res$txi[[lvl]], current_samples_1,
                         current_samples_2, ...)

    }
    de <- list()
    if (ncores == 1) {
        de$gene <- purrr::imap(res$de$gene, format_de_res, "gene")
        de$tx <- purrr::imap(res$de$tx, format_de_res, "tx")
    } else {
        mc_fun <- function(n, lvl) {
            format_de_res(res$de[[lvl]][[n]], n)
        }
        de$gene <- parallel::mclapply(names(res$de$gene), mc_fun,
                                      lvl = "gene", mc.cores = ncores)
        names(de$gene) <- names(res$de$gene)
        de$tx <- parallel::mclapply(names(res$de$tx), mc_fun,
                                      lvl = "tx", mc.cores = ncores)
        names(de$tx) <- names(res$de$tx)
    }

    stopifnot(all(names(de$gene) %in% names(de$tx)))
    stopifnot(all(names(de$tx) %in% names(de$gene)))
    save_de <- function(x) {
        output_de <- paste0(outdir, "/", prefix, "_", x, ".csv")
        current_de <- rbind(de$tx[[x]], de$gene[[x]])
        write_csv(current_de, output_de)
        current_de
    } 
    res_de <- map(names(de$gene), save_de)

    invisible(list(counts = counts, de = res_de))
}

#' Prepare formated table counts in odysse format.
#'
#' The table contains annotation and the counts (raw_counts/tpm/fpkm).
#'
#' @param txi The txi object returned by the import_kallisto function.
#' @param digits Integer indicating the number of decimal places
#' @param use_ruv Use RUVg normalization? Needs to be pre-computed using the
#'                \code{ruvg_normalization} function. Default: \code{FALSE}.
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' counts <- format_counts(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr everything
#'
#' @export
format_counts <- function(txi, digits = 4, use_ruv = FALSE) {
    validate_txi(txi)
    if (use_ruv) {
        stopifnot(!is.null(txi$ruvg_counts))
        stopifnot(identical(rownames(txi$counts), rownames(txi$ruvg_counts)))
        stopifnot(identical(colnames(txi$counts), colnames(txi$ruvg_counts)))
    }

    # Extract values
    if (!use_ruv) {
        raw_counts <- round(txi$counts, digits)
    } else {
        raw_counts <- round(txi$ruvg_counts, digits)
    }
    tpm <- round(txi$abundance, digits)
    fpkm <- round(txi$fpkm, digits)

    # Merge values
    res <- paste(paste(raw_counts, tpm, sep = "/"), fpkm, sep = "/") %>%
        matrix(ncol = ncol(raw_counts)) %>%
        as.data.frame %>%
        setNames(colnames(raw_counts)) %>%
        dplyr::mutate(id = rownames(raw_counts))

    # Add anno
    dplyr::left_join(res, txi$anno, by = "id") %>%
        dplyr::select(id:transcript_type, dplyr::everything())
}

#' Prepare formated DE table in odysse format.
#'
#' The table contains annotation and the DE results. The \code{de_res} input
#' object must be annotated and correctly formatted (i.e.: from
#' \code{produce_deliverables}
#'
#' @param de_res The results from DESeq2::results converted to data.frame.
#' @param txi The txi object returned by the import_kallisto function.
#' @param contrast The contrast for the comparison (see ?DESeq2::results).
#' @param design_group The column from the design where the contrasts are
#' located. Default: "group"
#' @param ignoreTxVersion Should the transcript version be ignored for anno
#' mapping. Default: \code{FALSE}.
#' @param digits Integer indicating the number of decimal places
#'
#' @return A data.frame with the anno and the merged counts values.
#'
#' @examples
#' txi <- get_demo_txi()
#' design <- get_demo_design()
#' dds <- deseq2_analysis(txi, design, ~ group)
#' de_res <- format_de_results(dds, txi, c("group", "A", "B"))
#' samples_grp1 <- design[design$group == "A",]$sample
#' samples_grp2 <- design[design$group == "B",]$sample
#' de <- format_de_odysse(de_res, txi, samples_grp1, samples_grp2)
#'
#' @importFrom magrittr %>%
#' @importFrom DESeq2 results
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom stringr str_replace
#'
#' @export
format_de_odysse <- function(de_res, txi, samples_grp1, samples_grp2,
                             digits = 4) {
    validate_txi(txi)
    validate_de(de_res, txi, keep_stats = TRUE)
    stopifnot(is(samples_grp1, "character"))
    stopifnot(length(samples_grp1) > 1)
    stopifnot(is(samples_grp2, "character"))
    stopifnot(length(samples_grp2) > 1)
    stopifnot(is(digits, "numeric"))
    stopifnot(identical(digits, round(digits, 0)))
    stopifnot(digits >= 0)
    if (!is.null(txi$dummy)) {
        stopifnot(is(txi$dummy, "character"))
        stopifnot("abundance" %in% txi$dummy)
    }

    de_res <- dplyr::mutate(de_res,
               mean_TPM_grp1 = rowMeans(txi$abundance[,samples_grp1]),
               mean_TPM_grp2 = rowMeans(txi$abundance[,samples_grp2]),
               ratio = 2^log2FoldChange,
               fold_change = if_else(ratio < 1, -1 * (1/ratio), ratio)) %>%
               splicing_analysis(txi)

    de_res <- dplyr::select(de_res, id, ensembl_gene, symbol, entrez_id,
                            transcript_type, mean_TPM_grp1, mean_TPM_grp2,
                            pV = pvalue, qV = padj, percent_grp1, percent_grp2,
                            main_isoform_grp1, main_isoform_grp2, baseMean,
                            lfcSE, log2FoldChange, fold_change, ratio, stat)

    de_res <- dplyr::mutate(de_res,
           mean_TPM_grp1 = round_values(as.numeric(mean_TPM_grp1), digits),
           mean_TPM_grp2 = round_values(as.numeric(mean_TPM_grp2), digits),
           pV = round_values(as.numeric(pV), digits),
           qV = round_values(as.numeric(qV), digits),
           percent_grp1 = round_values(as.numeric(percent_grp1), digits),
           percent_grp2 = round_values(as.numeric(percent_grp2), digits),
           baseMean = round_values(as.numeric(baseMean), digits),
           lfcSE = round_values(as.numeric(lfcSE), digits),
           fold_change = round_values(as.numeric(fold_change), digits),
           log2FoldChange = round_values(as.numeric(log2FoldChange), digits),
           stat = round_values(as.numeric(stat), digits))
    as.data.frame(de_res)
}

splicing_analysis <- function(res, txi) {
    is.max <- function(x) seq_along(x) == which.max(x)
    if (txi$txOut) {
        sums <- dplyr::group_by(res, ensembl_gene) %>%
            dplyr::summarize(sum_g1 = sum(mean_TPM_grp1),
            sum_g2 = sum(mean_TPM_grp2))
        dplyr::left_join(res, sums, by = "ensembl_gene") %>%
            dplyr::group_by(ensembl_gene) %>%
            dplyr::mutate(percent_grp1 = mean_TPM_grp1 / sum_g1,
                   percent_grp2 = mean_TPM_grp2 / sum_g2) %>%
            dplyr::mutate(percent_grp1 = if_else(is.nan(percent_grp1),
                                          0, percent_grp1),
                   percent_grp2 = if_else(is.nan(percent_grp2),
                                          0, percent_grp2)) %>%
            dplyr::mutate(main_isoform_grp1 = is.max(mean_TPM_grp1),
                   main_isoform_grp2 = is.max(mean_TPM_grp2)) %>%
            dplyr::select(-sum_g1, -sum_g2)
    } else {
        dplyr::mutate(res,
               percent_grp1 = NA,
               percent_grp2 = NA,
               main_isoform_grp1 = NA,
               main_isoform_grp2 = NA)
    }
}

round_values <- function(values, digits = 4) {
    stopifnot(is.numeric(values))
    i <- is.na(values)
    new_values <- round(values, digits) %>% format(scientific = FALSE)
    new_values[i] <- NA
    as.numeric(new_values)
}
