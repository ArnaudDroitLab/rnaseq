#' Produce deliverables 
#'
#' The goal of this function is to parse the variaus infos tables and to
#' produce and save the results in a standardized object and directories
#' format.
#' 
#' The results will be returned as a list with the following slots:
#' 
#' * txi$gene: The txi object at the gene level.
#' * txi$tx: The txi object at the transcript level.
#' * counts$gene: A list of count matrices at the gene level matching the count
#'                matrices found in the txi object.
#' * counts$tx: A list of count matrices at the transcript level matching the
#'              count matrices found in the txi object.
#' * pca$gene: The list of PCA ggplot objects at the gene level.
#' * pca$tx: The list of PCA ggplot objects at the transcript level.
#' * de$gene: The list of DE results in the gene format.
#' * de$tx: The list of DE results in the transcript format.
#' * de$volcano$gene: The volcano plots produced at the gene level.
#' * de$volcano$tx: The volcano plots produced at the transcript level.
#' * de$report: The lines of the report.
#' * de$infos_tables: The infos tables used to produce the results.
#'
#' It's the content of the infos_tables that will determine which deliverables
#' will be produced. Some tables are mandatory for other steps to be performed.
#' DE analysis requires de_infos and design_infos tables. Volcano analysis
#' requires the DE analysis to also be performed beforehand.
#'
#' @param dir_kallisto Directory with Kallisto quantifications
#' @param anno Path to the anno file matching the reference used for the
#' Kallisto analysis
#' @param infos_tables A \code{list} of infos tables, or path the csv files
#' that contains these infos. The list may contains the following tables:
#' pca_infos, de_infos, design_infos, volcano_infos. None of the tables are
#' mandatory, but only the analysis corresponding to the available tables will
#' be performed. See the \code{batch_pca}, \code{batch_de},
#' \code{batch_volcano} and the \code{produce_report} documentation to see the
#' description of each tables.
#' @param analysis_level Should the analysis be done at the "gene", "tx" or
#' "both" level. Default: "both".
#' @param file_type Abundance file format to use (h5 or tsv).
#' @param digits Integer indicating the number of decimal places
#' @param ignoreTxVersion Should the version number in the IDs be ignored
#'        during import. Default: \code{TRUE}.
#' @param counts Save count tables in the <outdir>/counts folder? Default:
#' \code{TRUE}.
#' @param add_volcano_labels Add labels to volcano plots. Must be either
#' \code{NULL} or a vector of gene symbols. Default: \code{NULL}.
#' @param add_volcanos_labels A vector of the symbols to show on the volcano
#' plots. If \code{NULL}, no symbols will be shown. Default: \code{NULL}.
#' @param report_filename The name of the report that will be saved in the
#' \code{outdir/reports/<report_filename>.Rmd}. If \code{NULL}, the report
#' won't be saved, but the lines will be available in the returned list in the
#' \code{report} slot. Default: \code{NULL}
#' @param outdir A \code{character} string corresponding to the directory to
#' store the output files (pdf, csv, etc.). The outputs will be in
#' <outdir>/counts, <outdir>/PCA, <outdir>/DE, <outdir/Volcanos,
#' <outdir>/report, if their respective tables are present in the
#' \code{infos_tables}. If \code{NULL}, the results won't be saved on disk.
#' Default: \code{NULL}.
#' @param r_objects A \code{character} string corresponding to the directory
#' where the rds files will be saved. A standardized data structure will be
#' used. The expected directories are <r_objects>/txi <r_objects>/counts
#' <r_objects>/de and <r_objects>/volcanos, each with a "gene/" and "tx/"
#' folder, except for the txi folder. The folders will only be present if the
#' respective table is present in the \code{infos_tables}. If \code{NULL}, the
#' rds files won't be saved on disk. It is recommended to use the
#' \code{r_objects} param to save files as every steps are completed in case a
#' problem cause the function to stop. That way, you won't have to reprocess
#' the steps that completed correctly.  Default: "r_objects".
#' @param force Should the files be re-created if they already exists? Default:
#' \code{FALSE}.
#' @param ncores Number of cores to use for de analysis. Default \code{1}.
#' @param verbose Print progression of the analysis? Default: \code{FALSE}.
#'
#' @return Invisibly returns a \code{list} that contains all the results.
#'
#' @importFrom readr read_csv
#' @importFrom dplyr if_else
#'
#' @export
produce_deliverables <- function(dir_kallisto,
                                 anno,
                                 infos_tables,
                                 analysis_level = "both",
                                 file_type = "h5",
                                 digits = 4,
                                 ignoreTxVersion = TRUE,
                                 counts = TRUE,
                                 add_volcano_labels = NULL,
                                 report_filename = NULL,
                                 outdir = "deliverables",
                                 r_objects = "r_objects",
                                 force = FALSE,
                                 ncores = 1,
                                 verbose = FALSE) {

    # Note: I removed the normalizations (i.e.: ruvg, combat_seq) from this
    # function. If the user want to produce the normalizations, the txi will
    # have to be created before launching the produce_deliverables function and
    # saved in the <r_objects>/txi/txi.rds file.

    print_verbose("Validating parameters...", verbose)
    # Parameters validation
    stopifnot(dir.exists(dir_kallisto))
    stopifnot(is(analysis_level, "character"))
    stopifnot(analysis_level %in% c("gene", "tx", "both"))
    stopifnot(file_type %in% c("tsv", "h5"))
    stopifnot(file.exists(anno))
    stopifnot(is.logical(ignoreTxVersion))
    stopifnot(!is.na(ignoreTxVersion))
    stopifnot(is.logical(counts))
    stopifnot(is.numeric(ncores))
    stopifnot(round(ncores) == ncores)
    stopifnot(ncores > 0)

    if (!is.null(add_volcano_labels)) {
        stopifnot(is(add_volcano_labels, "character"))
        stopifnot(all(nchar(add_volcano_labels) > 0))
    }

    stopifnot(is(infos_tables, "list"))

    expected_names <- c("pca_infos", "de_infos", "design_infos",
                        "volcano_infos", "report_infos", "metadata")
    for (n in expected_names) {
        if (!is.null(infos_tables[[n]])) {
            if (is(infos_tables[[n]], "character")) {
                stopifnot(file.exists(infos_tables[[n]]))
                infos_tables[[n]] <- readr::read_csv(infos_tables[[n]])
            }
            stopifnot(is(infos_tables[[n]], "data.frame"))
            stopifnot(nrow(infos_tables[[n]]) > 0)
        }
    }
    stopifnot(any(expected_names %in% names(infos_tables)))
    if ("de_infos" %in% names(infos_tables)) {
        stopifnot("design_infos" %in% names(infos_tables))
    }
    if ("volcano_infos" %in% names(infos_tables)) {
        stopifnot("de_infos" %in% names(infos_tables))
    }

#    infos_tables <- complete_and_validate_tables(infos_tables)
    print_verbose("    Done!", verbose)

    ## Directories
    print_verbose("Creating directory structure...", verbose)
    if (!is.null(outdir)) {
        stopifnot(is(outdir, "character"))
        if (!dir.exists(outdir)) {
            dir.create(outdir)
        }
    }
    if (!is.null(r_objects)) {
        stopifnot(is(r_objects, "character"))
        if (!dir.exists(r_objects)) {
            dir.create(r_objects)
        }
    }

    res <- list()

    ## Import quantifications
    print_verbose("Importing quantification results...", verbose)
    if (!is.null(r_objects)) {
        msg <- "Checking if r_objects files already present..."
        print_verbose(paste0("    ", msg), verbose)

        r_objects_txi <- paste0(r_objects, "/txi")
        output_txi_rds <- paste0(r_objects_txi, "/txi.rds")
        if (!dir.exists(r_objects_txi)) {
            dir.create(r_objects_txi)
        }
        if (!force & file.exists(output_txi_rds)) {
            msg <- "Files are present, importing rds files..."
            print_verbose(paste0("        ", msg), verbose)

            res$txi <- readRDS(output_txi_rds)
            validate_txi(res$txi$gene)
            validate_txi(res$txi$tx)
        } else {
            if (file.exists(output_txi_rds) & force) {
                msg <- "Files are present, but force is TRUE.\n"
                msg <- paste0(msg, "Re-creating txi...")
                print_verbose(paste0("        ", msg), verbose)
            } else {
                msg <- "Files not found, creating txi objects..."
                print_verbose(paste0("        ", msg), verbose)
            }

            files <- get_filenames(dir_kallisto, file_type)
            res$txi <- produce_txi(files = files, anno = anno,
                               ignoreTxVersion = ignoreTxVersion)
            saveRDS(res$txi, output_txi_rds)
        }
    } else {
        msg <- "r_objects is set to NULL, creating txi objects "
        msg <- paste0(msg, "without saving to disk...")
        print_verbose(paste0("    ", msg), verbose)
        files <- get_filenames(dir_kallisto, file_type)
        res$txi <- produce_txi(files = files, anno = anno,
                           ignoreTxVersion = ignoreTxVersion)
    }

    res$counts <- list()
    ## Counts
    print_verbose("Counts extraction...", verbose)
    save_counts <- function(lvl) {
        msg <- paste0("    Current level: ", lvl)
        print_verbose(msg, verbose)
        current_txi <- res$txi[[lvl]]
        current_res <- list()
        if (counts) {
            if (!is.null(outdir)) {
                outdir_counts <- paste0(outdir, "/counts")
                if (!dir.exists(outdir_counts)) {
                    dir.create(outdir_counts)
                }
            }
            expected_matrices <- c("counts", "abundance", "fpkm", "ruvg_counts",
                                   "combat_counts", "extra_count_matrix")
            for (m in expected_matrices) {
                msg <- paste0("        Current matrix: ", m)
                print_verbose(msg, verbose)
                if (!m %in% current_txi$dummy) {
                    if (!is.null(current_txi[[m]])) {
                        df <- get_anno_df(current_txi, m)
                        if (!is.null(outdir)) {
                            current_count_file <- paste0(outdir_counts, "/", m,
                                                         "_", lvl, ".csv")
                            if (force | !file.exists(current_count_file)) {
                                msg <- "            Saving to disk."
                                print_verbose(msg, verbose)
                                readr::write_csv(df, current_count_file)
                            } else {
                                msg <- "            File is already present"
                                msg <- paste0(msg, " and force is FALSE.")
                                msg <- paste0(msg, " Not saving to disk.")
                                print_verbose(msg, verbose)
                            }
                        } else {
                            msg <- "            outdir is FALSE.\n"
                            msg <- paste0(msg,
                                          "            Not saving to disk.")
                            print_verbose(msg, verbose)
                        }
                        current_res[[m]] <- df
                    } else {
                        msg <- "            Matrix not found in txi. Skipping."
                        print_verbose(msg, verbose)
                    }
                } else {
                    msg <- paste0("        Matrix ", m,
                                  " is a dummy matrix. Skipping")
                    print_verbose(msg, verbose)
                }
            }
        }
        current_res
    }
    if (analysis_level %in% c("both", "gene")) {
        res$counts[["gene"]] <- save_counts("gene")
    }
    if (analysis_level %in% c("both", "tx")) {
        res$counts[["tx"]] <- save_counts("tx")
    }

    ## PCA
    print_verbose("PCA analysis...", verbose)
    if (!is.null(infos_tables$pca_infos)) { 
        msg <- "pca_infos table found."
        print_verbose(paste0("    ", msg), verbose)
        msg <- "Preparing directory structure..."
        print_verbose(paste0("    ", msg), verbose)

        outdir_pca <- paste0(outdir, "/pca")
        r_objects_pca <- paste0(r_objects, "/pca")
        if (!dir.exists(outdir_pca)) {
            dir.create(outdir_pca)
        }
        if (!dir.exists(r_objects_pca)) {
            dir.create(r_objects_pca)
        }
        res$pca <- list()
        if (analysis_level %in% c("both", "gene")) {
            msg <- "Launching batch_pca at gene levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$pca$gene <- batch_pca(pca_infos = infos_tables$pca_infos,
                                      txi = res$txi$gene,
                                      metadata = infos_tables$metadata,
                                      outdir = paste0(outdir_pca, "/gene"),
                                      r_objects = paste0(r_objects_pca,
                                                         "/gene"),
                                      force = force,
                                      cores = ncores)
        }
        if (analysis_level %in% c("both", "tx")) {
            msg <- "Launching batch_pca at transcripts levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$pca$tx <- batch_pca(pca_infos = infos_tables$pca_infos,
                                    txi = res$txi$tx,
                                    metadata = infos_tables$metadata,
                                    outdir = paste0(outdir_pca, "/tx"),
                                    r_objects = paste0(r_objects_pca, "/tx"),
                                    force = force,
                                    cores = ncores)
        }
    } else {
        msg <- "pca_infos table not found, skipping PCA analysis."
        print_verbose(paste0("    ", msg), verbose)
    }


    # DE
    print_verbose("Differential expression (DE) analysis...", verbose)
    if (!is.null(infos_tables$de_infos)) {
        msg <- "de_infos table found."
        print_verbose(paste0("    ", msg), verbose)
        msg <- "Preparing directory structure..."
        print_verbose(paste0("    ", msg), verbose)

        outdir_de <- paste0(outdir, "/de")
        r_objects_de <- paste0(r_objects, "/de")
        if (!dir.exists(outdir_de)) {
            dir.create(outdir_de)
        }
        if (!dir.exists(r_objects_de)) {
            dir.create(r_objects_de)
        }
        res$de <- list()
        if (analysis_level %in% c("both", "gene")) {
            msg <- "Launching batch_de at gene levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$de$gene <- batch_de(de_infos = infos_tables$de_infos,
                                    txi = res$txi$gene,
                                    design = infos_tables$design_infos,
                                    outdir = paste0(outdir_de, "/gene"),
                                    r_objects = paste0(r_objects_de, "/gene"),
                                    force = force,
                                    cores = ncores)
        }
        if (analysis_level %in% c("both", "tx")) {
            msg <- "Launching batch_de at transcript levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$de$tx <- batch_de(de_infos = infos_tables$de_infos,
                                  txi = res$txi$tx,
                                  design = infos_tables$design_infos,
                                  outdir = paste0(outdir_de, "/tx"),
                                  r_objects = paste0(r_objects_de, "/tx"),
                                  force = force,
                                  cores = ncores)
        }
    } else {
        msg <- "de_infos table not found, skipping DE analysis."
        print_verbose(paste0("    ", msg), verbose)
    }

    # Volcanos
    print_verbose("Volcano analysis...", verbose)
    if (!is.null(infos_tables$volcano_infos)) {
        msg <- "volcano_infos table found."
        print_verbose(paste0("    ", msg), verbose)
        msg <- "Preparing directory structure..."
        print_verbose(paste0("    ", msg), verbose)

        outdir_volcano <- paste0(outdir, "/volcano")
        r_objects_volcano <- paste0(r_objects, "/volcano")
        if (!dir.exists(outdir_volcano)) {
            dir.create(outdir_volcano)
        }
        if (!dir.exists(r_objects_volcano)) {
            dir.create(r_objects_volcano)
        }
        res$volcano <- list()
        if (analysis_level %in% c("both", "gene")) {
            msg <- "Launching batch_volcano at gene levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$volcano$gene <- batch_volcano(volcano_infos =
                                              infos_tables$volcano_infos,
                                              de_results = res$de$gene,
                                              add_labels = add_volcano_labels,
                                              outdir = paste0(outdir_volcano,
                                                              "/gene"),
                                              r_objects =
                                                  paste0(r_objects_volcano,
                                                         "/gene"),
                                              force = force,
                                              cores = ncores)
        }
        if (analysis_level %in% c("both", "tx")) {
            msg <- "Launching batch_volcano at transcript levels..."
            print_verbose(paste0("    ", msg), verbose)
            res$volcano$tx <- batch_volcano(volcano_infos =
                                            infos_tables$volcano_infos,
                                            de_results = res$de$tx,
                                            add_labels = add_volcano_labels,
                                            outdir = paste0(outdir_volcano,
                                                            "/tx"),
                                            r_objects =
                                                paste0(r_objects_volcano,
                                                       "/tx"),
                                            force = force,
                                            cores = ncores)
        }
    } else {
        msg <- "volcano_infos table not found, skipping volcano analysis."
        print_verbose(paste0("    ", msg), verbose)
    }

    # reports
    print_verbose("Report analysis...", verbose)
    if (!is.null(infos_tables$report_infos)) {
        msg <- "report_infos table found."
        print_verbose(paste0("    ", msg), verbose)
        msg <- "Starting report creation..."
        print_verbose(paste0("    ", msg), verbose)

        if (!is.null(report_filename)) {
            msg <- "report_filename is available"
            print_verbose(paste0("    ", msg), verbose)
            msg <- paste0("Report will be saved as: ", report_filename)
            print_verbose(paste0("    ", msg), verbose)
        } else {
            msg <- "report_filename is NULL"
            print_verbose(paste0("    ", msg), verbose)
            msg <- "Report won't be saved to file."
            print_verbose(paste0("    ", msg), verbose)
        }

        res$report <- produce_report(report_infos = infos_tables$report_infos,
                                     report_filename = report_filename)
    } else {
        msg <- "report_infos table not found."
        print_verbose(paste0("    ", msg), verbose)
        msg <- "Skipping report production."
        print_verbose(paste0("    ", msg), verbose)
    }

    print_verbose("produce_deliverables completed successfully!", verbose)
    res$infos_tables <- infos_tables
    invisible(res)
}

print_verbose <- function(msg, verbose) {
    if (verbose) {
        message(msg)
    }
}

complete_and_validate_tables <- function(infos_tables) {
    it <- infos_tables

    if (!is.null(it$pca_infos)) {
        it$pca_infos <- complete_pca_infos(it$pca_infos)
        stopifnot(length(validate_pca_infos(it$pca_infos, it$metadata, res$txi)) == 0)
    }

    if (!is.null(it$de_infos)) {
        it$de_infos <- complete_de_infos(it$de_infos)
        stopifnot(length(validate_de_infos(it$de_infos, it$design_infos, res$txi)) == 0)
    }
    
    if (!is.null(it$volcano_infos)) {
        it$volcano_infos <- complete_volcano_infos(it$volcano_infos)
        stopifnot(length(validate_volcano_infos(it$volcano_infos, de_results)) == 0)
    }
    
    it
}
