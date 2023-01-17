#' Import quantifications from Kallisto
#'
#' @param filenames Paths to the abundance files.
#' @param anno The filename for the annotation in csv (see anno package)
#'  It must have the following columns: id, ensembl_gene, symbol, entrez_id and
#'  transcript_type. Values can be NA except for id and ensembl_gene.
#' @param txOut Return counts and abundance at the transcript level. Default:
#'              FALSE
#' @param ignoreTxVersion Ignore version of tx. Default = FALSE
#'
#' @return A txi object.
#'
#' @examples
#' abundances <- get_demo_abundance_files()
#' file_anno <- get_demo_anno_file()
#' txi <- import_kallisto(abundances, anno = file_anno)
#'
#' @importFrom magrittr %>%
#' @importFrom stringr str_extract
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom readr read_csv
#' @importFrom tximport tximport
#' @importFrom tximport summarizeToGene
#'
#' @export
import_kallisto <- function(filenames, anno, txOut = FALSE,
                            ignoreTxVersion = FALSE) {
    stopifnot(all(file.exists(filenames)))
    stopifnot(all(file.exists(anno)))
    stopifnot(txOut %in% c(TRUE, FALSE))
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    tx2gene <- readr::read_csv(anno, col_types = "ccccc") %>%
        as.data.frame %>%
        dplyr::select(TXNAME = id, GENEID = ensembl_gene)
    if (txOut == TRUE) {
        txi <- tximport::tximport(filenames, type = "kallisto",
                                  tx2gene = tx2gene, txOut = TRUE,
                                  ignoreTxVersion = ignoreTxVersion)
    } else {
        txi <- tximport::tximport(filenames, type = "kallisto",
                                  tx2gene = tx2gene,
                                  ignoreTxVersion = ignoreTxVersion)
    }
    txi$fpkm <- get_fpkm(txi)
    txi$anno <- get_anno(anno, txOut)
    expected_colnames <- c("id" , "ensembl_gene", "symbol", "entrez_id", "transcript_type")

    stopifnot(identical(colnames(txi$anno), expected_colnames))
    stopifnot(sum(is.na(txi$anno$id)) == 0)
    stopifnot(sum(is.na(txi$anno$ensembl_gene)) == 0)
    stopifnot(length(unique(txi$anno$id)) == length(txi$anno$id))

    txi$txOut <- txOut
    if (!ignoreTxVersion) {
        stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    } else {
        rownames(txi$counts) <- stringr::str_extract(rownames(txi$counts), "^[^\\.]*")
        rownames(txi$abundance) <- stringr::str_extract(rownames(txi$abundance), "^[^\\.]*")
        rownames(txi$fpkm) <- stringr::str_extract(rownames(txi$fpkm), "^[^\\.]*")
        txi$anno$id <- stringr::str_extract(txi$anno$id, "^[^\\.]*")
        txi$anno$ensembl_gene <- stringr::str_extract(txi$anno$ensembl_gene, "^[^\\.]*")
        stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    }
    stopifnot(identical(nrow(txi$counts), nrow(txi$anno)))
    stopifnot(all(rownames(txi$counts) %in% txi$anno$id))
    stopifnot(all(txi$anno$id %in% rownames(txi$counts)))

    txi
}

summarize_to_gene <- function(txi_tx, anno, ignoreTxVersion = FALSE) {
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    if (!file.exists(anno)) {
        tx2gene <- get(anno) %>%
            dplyr::select(TXNAME = id, GENEID = ensembl_gene)
    } else {
        tx2gene <- readr::read_csv(anno, col_types = "ccccc") %>%
            as.data.frame %>%
            dplyr::select(TXNAME = id, GENEID = ensembl_gene)
    }

    txi <- tximport::summarizeToGene(txi_tx, tx2gene = tx2gene,
                                     ignoreTxVersion = ignoreTxVersion)
    txi$fpkm <- get_fpkm(txi)
    txi$anno <- get_anno(anno, txOut = FALSE)
    txi$txOut <- FALSE
    stopifnot(all(rownames(txi$fpkm) %in% txi$anno$id))
    txi
}

arrange_anno <- function(anno) {
    lvls <- c("protein_coding", "miRNA", "lincRNA", "antisense", "snRNA",
              "lncRNA", "scRNA", "vault_RNA", "Spike-in", "sense_intronic",
              "misc_RNA", "snoRNA", "scaRNA", "rRNA",
              "3prime_overlapping_ncrna", "sense_overlapping", "sRNA",
              "ribozyme", "vaultRNA", "non_coding", "macro_lncRNA", "Mt_tRNA",
              "Mt_rRNA", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
              "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "IG_LV_gene",
              "IG_D_pseudogene", "processed_transcript",
              "translated_processed_pseudogene",
              "translated_unprocessed_pseudogene",
              "transcribed_processed_pseudogene",
              "transcribed_unitary_pseudogene",
              "transcribed_unprocessed_pseudogene", "pseudogene",
              "unitary_pseudogene", "unprocessed_pseudogene",
              "polymorphic_pseudogene", "processed_pseudogene",
              "IG_J_pseudogene", "IG_V_pseudogene", "IG_pseudogene",
              "IG_C_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
              "rRNA_pseudogene", "retained_intron", "non_stop_decay",
              "nonsense_mediated_decay", "TEC",
              "protein_coding_CDS_not_defined", "protein_coding_LoF")
    stopifnot(all(anno$transcript_type %in% lvls))

    dplyr::mutate(anno, transcript_type = factor(transcript_type,
                                                 levels = lvls)) %>%
            dplyr::arrange(id, transcript_type) %>%
            dplyr::mutate(transcript_type = as.character(transcript_type))
}

get_anno <- function(anno, txOut = TRUE) {
#    validate_anno(anno)
    anno <- readr::read_csv(anno, col_types = "ccccc")
    if (!txOut) {
        anno <- dplyr::mutate(anno, id = ensembl_gene) %>%
            arrange_anno %>%
            dplyr::filter(!duplicated(ensembl_gene))
    }
    as.data.frame(anno)
}

validate_anno <- function(anno) {
    valid_anno <- c("Hs.Gencode19", "Hs.Gencode27", "Hs.Gencode32",
                    "Hs.Gencode35", "Hs.Gencode37", "Hs.Ensembl79",
                    "Hs.Ensembl91", "Hs.Ensembl95", "Hs.Ensembl97",
                    "Hs.Ensembl98", "Hs.Ensembl100", "Hs.Ensembl101",
                    "Mm.Ensembl91", "Mm.Ensembl92", "Mm.Ensembl94",
                    "Mm.Ensembl99", "Mm.Ensembl100", "Mm.Ensembl101",
                    "Rn.Ensembl76", "Rn.Ensembl79", "Rn.Ensembl92",
                    "Rn.Ensembl98", "Bt.Ensembl99", "Mmu.Ensembl101",
                    "Mmu.Ensembl103", "peaux_colonisees")
    stopifnot(anno %in% valid_anno)
}

get_fpkm <- function(txi) {
    (txi$counts * 10^6) / (colSums(txi$counts) * (txi$length/1000))
}

#' Parse kallisto quant directory to retrieve filenames
#'
#' @param dir_kallisto Kallisto quantification base directory
#' @param file_extension The extension to parse (h5 or tsv). Default: h5
#'
#' @return A named list of the abundance files.
#'
#' @examples
#' dir_kallisto <- get_demo_kallisto_dir()
#' filenames <- get_filenames(dir_kallisto, file_extension = "tsv")
#'
#' @export
get_filenames <- function(dir_kallisto, file_extension = "h5") {
    stopifnot(dir.exists(dir_kallisto))
    stopifnot(is.character(file_extension))
    stopifnot(length(file_extension) == 1)
    stopifnot(file_extension %in% c("h5", "tsv"))
    filenames <- dir(dir_kallisto,
                     pattern = file_extension,
                     full.names = TRUE, recursive = TRUE)
    stopifnot(length(filenames) >= 1)
    names(filenames) <- basename(dirname(filenames))
    filenames
}
