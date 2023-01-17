#' Produce report in Rmd format
#'
#' The function will build a rmarkdown file based on the the
#' \code{report_infos} table. The table must contain at least the following
#' columns: \code{add} and \code{value}. The \code{add} must be one of the
#' following: text, file or plot. The value will vary depending on the type of
#' add.
#'
#' You can also add the \code{id_report} column that will need to contains a
#' unique identifiers for each line of the report_infos table.
#'
#' For text, it must be a character string corresponding to the text that
#' will be directly added to the file specified with \code{report_filename}
#' file.
#'
#' For file, it must be a valid path to a filename
#'
#' For plot, it must be a valid path to a \code{rds} file that contains the
#' plot to show. The plot must be into a format that can be shown using the
#' \code{print} function
#'
#' The \code{extra} column will be used only when the \code{add} value is
#' "plot" and will be added as is to the code chunk. Be aware that not validity
#' checks are performed. Make sure you use valid chunk option syntax.
#'
#' The content of the rmarkdown will be created by parsing each line of the
#' \code{report_infos} table and adding the requested elements in the same
#' order as they are found in the table.
#'
#' There will be no validation that the created rmarkdown file is in the
#' correct format!
#'
#' @param report_infos A \code{data.frame} or the path to a csv file that
#' describes the report to produce.
#' @param report_filename The name of the rmarkdown file to create. If
#' \code{NULL}, the report won't be saved to file. Default: report.Rmd
#' @param verbose Print progression? Default: \code{FALSE}.
#'
#' @return Invisibly returns the lines saved to the Rmd file.
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_detect
#' @importFrom R.utils getRelativePath
#'
#' @export
produce_report <- function(report_infos, report_filename = "report.Rmd", verbose = FALSE) {

    # 1. Data validation
    stopifnot(is(verbose, "logical"))
    msg <- "Starting data validation..."
    print_verbose(msg, verbose)
    stopifnot(is(report_infos, "data.frame") | is(report_infos, "character"))
    if (is.character(report_infos)) {
        stopifnot(file.exists(report_infos))
        report_infos <- readr::read_csv(report_infos, show_col_types = FALSE)
    }
    stopifnot(nrow(report_infos) > 0)
    stopifnot(all(c("add", "value") %in% colnames(report_infos)))
    report_infos <- complete_report_infos(report_infos)
    stopifnot(is(report_infos$id_report, "character"))
    stopifnot(!any(duplicated(report_infos$id_report)))
    stopifnot(is(report_infos$extra, "character"))
    stopifnot(all(report_infos$add %in% c("text", "file", "plot")))
    all_ids <- dplyr::filter(report_infos, add == "plot") %>%
            dplyr::pull(value) %>%
            basename %>%
            tools::file_path_sans_ext()
    stopifnot(!any(duplicated(all_ids)))
    for (i in 1:nrow(report_infos)) {
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]
        if (current_add == "text") {
            stopifnot(is(current_value, "character"))
        } else if (current_add == "file") {
            stopifnot(file.exists(current_value))
        } else if (current_add == "plot") {
            stopifnot(file.exists(current_value))
            stopifnot(stringr::str_detect(current_value, "\\.rds$"))
        }
    }
    if (!is.null(report_filename)) {
        stopifnot(is(report_filename, "character"))
        if (dirname(report_filename) != ".") {
            stopifnot(dir.exists(dirname(report_filename)))
        }
    }
    msg <- "Done!\n"
    print_verbose(msg, verbose)

    # 2. Parse report_infos
    msg <- "Parsing report_infos line by line..."
    print_verbose(msg, verbose)
    produce_lines <- function(i) {
        if (i %% 50 == 0) {
            msg <- paste0("    Parsed: ", i, " lines")
            print_verbose(msg, verbose)
        }
        current_add <- report_infos$add[i]
        current_value <- report_infos$value[i]
        current_extra <- report_infos$extra[i]

        if (current_add == "text") {
            lines <- c(current_value, "\n")
        } else if (current_add == "file") {
            lines <- c(readLines(current_value), "\n")
        } else if (current_add == "plot") {
            correct_path_value <- R.utils::getRelativePath(current_value, dirname(report_filename))
            current_id <- basename(current_value) %>%
                tools::file_path_sans_ext()
            current_line <- paste0("```{r ", current_id, " ", current_extra, "}\n")
            current_line <- paste0(current_line, "print(readRDS('", correct_path_value, "'))", "\n")
            current_line <- paste0(current_line, "```")
            lines <- c(current_line, "\n")
        }
        lines
    }
    all_lines <- purrr::map(1:nrow(report_infos), produce_lines) %>% unlist
    msg <- "Parsing report_infos line by line... Done!"
    print_verbose(msg, verbose)

    msg <- "Printing Rmd file..."
    print_verbose(msg, verbose)

    if (!is.null(report_filename)) {
        file_conn<-file(report_filename)
        writeLines(all_lines, file_conn)
        close(file_conn)
    }
    msg <- "Done!"
    print_verbose(msg, verbose)

    invisible(lines)
}

complete_report_infos <- function(report_infos) {
    # Fill missing
    if (!"id_report" %in% colnames(report_infos))
        report_infos$id_report <- paste("report_", 1:nrow(report_infos))
    if (!"extra" %in% colnames(report_infos))
        report_infos$extra <- ""
    report_infos
}
