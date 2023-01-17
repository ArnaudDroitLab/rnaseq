#' Extract raw counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' raw_counts <- get_raw_count_anno_df(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join 
#' @importFrom dplyr select 
#' @importFrom dplyr one_of 
#' @importFrom dplyr everything 
#'
#' @export
get_raw_count_anno_df <- function(txi) {
    as.data.frame(txi$counts) %>%
        tibble::rownames_to_column("id") %>%
        dplyr::left_join(txi$anno, by = "id") %>%
        dplyr::select(dplyr::one_of(colnames(txi$anno)),
                      dplyr::everything())
}

#' Extract TPM with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' tpm <- get_tpm_anno_df(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr select 
#' @importFrom dplyr one_of 
#' @importFrom dplyr everything 
#'
#' @export
get_tpm_anno_df <- function(txi) {
    as.data.frame(txi$abundance) %>%
        tibble::rownames_to_column("id") %>%
        dplyr::left_join(txi$anno, by = "id") %>%
        dplyr::select(dplyr::one_of(colnames(txi$anno)),
                      dplyr::everything())
}

#' Extract RUVg counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function and
#'             the `ruvg_normalisation` function
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' ruvg <- get_ruvg_anno_df(txi)
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr select 
#' @importFrom dplyr one_of 
#' @importFrom dplyr everything 
#'
#' @export
get_ruvg_anno_df <- function(txi) {
    as.data.frame(txi$ruvg_counts) %>%
        tibble::rownames_to_column("id") %>%
        dplyr::left_join(txi$anno, by = "id") %>%
        dplyr::select(dplyr::one_of(colnames(txi$anno)),
                      dplyr::everything())
}

#' Extract counts with annotation data.frame from txi
#'
#' @param txi: The txi object returned from the `import_kallisto` function
#' @param table_name: The name of the count table to return. Tolerated values:
#'            * For raw counts: counts, raw_counts, count and raw_count
#'            * For TPM: abundance, TPM and tpm
#'            * For RUVg: ruvg_counts, ruvg, RUVg, RUV and ruv
#'            * For FPKM: fpkm and FPKM
#'
#' @return A `data.frame` object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#' raw_counts <- get_anno_df(txi, "raw_count")
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr everything
#'
#' @export
get_anno_df <- function(txi, colname) {
    valid_colname <- c("count", "raw_count", "counts", "raw_counts", "tpm",
                       "abundance", "TPM", "ruvg", "RUVg", "RUV", "ruv",
                       "ruvg_counts", "fpkm", "FPKM", "combat", "Combat",
                        "combat_counts", "extra_count_matrix")
    stopifnot(colname %in% valid_colname)
    if (colname %in% c("count", "raw_count", "counts", "raw_counts")) {
        colname <- "counts"
    }
    if (colname %in% c("tpm", "abundance", "TPM")) {
        colname <- "abundance"
    }
    if (colname %in% c("ruvg_counts", "ruvg", "RUVg", "RUV", "ruv")) {
        colname <- "ruvg_counts"
    }
    if (colname %in% c("fpkm", "FPKM")) {
        colname <- "fpkm"
    }
    if (colname %in% c("combat", "Combat", "combat_counts")) {
        colname <- "combat_counts"
    }
    as.data.frame(txi[[colname]]) %>%
        tibble::rownames_to_column("id") %>%
        dplyr::left_join(txi$anno, by = "id") %>%
        dplyr::select(dplyr::one_of(colnames(txi$anno)),
                      dplyr::everything())
}
