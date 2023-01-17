#' Get the list of default human housekeeping genes.
#'
#' Data from de Jonge *et al.*, 2007:
#' Evidence Based Selection of Housekeeping Genes
#' 
#' @return A \code{vector} with human housekeeping genes symbols
#' 
#' @examples
#' human_hsk <- get_human_hsk()
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_csv
#' @importFrom dplyr pull
#'
#' @export
get_human_hsk <- function() {
    file_hsk <- system.file("extdata/human_hsk.csv", package = "rnaseq")
    readr::read_csv(file_hsk, col_type = "c") %>%
        dplyr::pull(Symbol)
}

#' Normalize data with RUVSeq RUVg algorithm
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#'            function.
#' @param housekeeping_genes A \code{vector} of gene symbols
#' @param ignoreTxVersion Ignore version of tx. Default = FALSE
#' @param txOut txi was produced at tx level?. Default = FALSE
#'
#' @return The original txi object with a \code{$ruvg_counts} element added.
#'
#' @examples
#' txi <- get_demo_txi()
#' # For this demo, we use a subset of housekeeping genes. In a real example,
#' # It is recommended to use the default value or a complete set of
#' # housekeeping genes.
#' txi_ruv <- ruvg_normalization(txi, housekeeping_genes = c("RPL9", "RPL24"))
#'
#' @importFrom EDASeq newSeqExpressionSet
#' @importFrom EDASeq betweenLaneNormalization
#' @importFrom EDASeq normCounts
#' @importFrom RUVSeq RUVg
#' @importFrom magrittr %>%
#' @importFrom stringr str_replace
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#'
#' @export
ruvg_normalization <- function(txi, housekeeping_genes = get_human_hsk(),
                               ignoreTxVersion = FALSE, txOut = FALSE) {
    stopifnot(is(housekeeping_genes, "character"))
    stopifnot(ignoreTxVersion %in% c(TRUE, FALSE))

    MyData <- txi$counts %>% round
    MyData <- MyData + 1
    if (ignoreTxVersion) {
        rownames(MyData) <- stringr::str_replace(rownames(MyData), "\\..*$", "")
    }

    stopifnot(all(housekeeping_genes %in% txi$anno$symbol))

    hskList <- data.frame(symbol = housekeeping_genes) %>%
            dplyr::left_join(txi$anno, by = "symbol")
    if (txOut) {
            hskList <- dplyr::pull(hskList, id)
    } else {
            hskList <- dplyr::pull(hskList, ensembl_gene) %>% unique
    }

    filter <- apply(MyData, 1, function(x) length(x[x>1])>=2)
    MyData <-  MyData[filter,]
    genes <- rownames(MyData[!(rownames(MyData) %in% hskList), ])
    spikes <- rownames(MyData[rownames(MyData) %in% hskList,])
    instances <- as.factor(colnames(MyData))

    phenoData <- data.frame(instances, row.names=colnames(MyData))
    txi$ruv_sets <- list()
    txi$ruv_sets$set <- EDASeq::newSeqExpressionSet(as.matrix(MyData), phenoData = phenoData)
    txi$ruv_sets$set_UQ <- EDASeq::betweenLaneNormalization(txi$ruv_sets$set, which="upper")
    txi$ruv_sets$set_RUVg <- RUVSeq::RUVg(txi$ruv_sets$set_UQ, spikes, k=1)

    txi$ruvg_counts <- EDASeq::normCounts(txi$ruv_sets$set_RUVg)
    txi
}

#' Produce RUVSeq recommended graph to evaluate normalization
#'
#' @param txi The \code{txi} object returned by the \code{ruvg_normalization}
#'            function.
#' @param output The name of the pdf that will be produced.
#'               Default: RUVg_graphs.pdf
#' 
#' @return Silently return the \code{txi} object used as input.
#' 
#' @examples
#' txi <- get_demo_txi()
#' # For this demo, we use a subset of housekeeping genes. In a real example,
#' # It is recommended to use the default value or a complete set of
#' # housekeeping genes.
#' txi_ruv <- ruvg_normalization(txi, housekeeping_genes = c("RPL9", "RPL24"))
#' \dontrun{
#'    produce_RUVSeq_graphs(txi_ruv)
#' }
#'
#' @importFrom EDASeq plotRLE
#' @importFrom EDASeq plotPCA
#'
#' @export
produce_RUVSeq_graphs <- function(txi, output = "RUVg_graphs.pdf") {
    stopifnot("ruv_sets" %in% names(txi))
    stopifnot("set" %in% names(txi$ruv_sets))
    stopifnot("set_UQ" %in% names(txi$ruv_sets))
    stopifnot("set_RUVg" %in% names(txi$ruv_sets))

    plotRelativeLogExpression = function(set, title){
        EDASeq::plotRLE(set, outline=FALSE, ylim=c(-2.5, 2.5), main=title, xlab="Patients", las=2)
    }
    plotPrincipalComponent = function(set, title){
        EDASeq::plotPCA(set, cex=1.2)
    }

    pdf(output)
    #No normalization
    plotRelativeLogExpression(txi$ruv_sets$set,"No_normalization")
    plotPrincipalComponent(txi$ruv_sets$set,"No_normalization")

    #betweenLaneNormalization
    plotRelativeLogExpression(txi$ruv_sets$set_UQ, "Upper-quartile_normalization")
    plotPrincipalComponent(txi$ruv_sets$set_UQ,"Upper-quartile_normalization")

    #Remove Unwanted Variations
    plotRelativeLogExpression(txi$ruv_sets$set_RUVg, "RUVg_normalization")
    plotPrincipalComponent(txi$ruv_sets$set_RUVg,"RUVg_normalization")
    dev.off()

    invisible(txi)
}
