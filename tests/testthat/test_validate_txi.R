library(rnaseq)

valid_txi <- suppressMessages(get_demo_txi(large = TRUE))
valid_txi_ercc <- suppressMessages(get_demo_txi(large = TRUE, ercc = TRUE))
valid_txi_tx <- suppressMessages(get_demo_txi(large = TRUE, txOut = TRUE))
valid_txi_tx_ercc <- suppressMessages(get_demo_txi(large = TRUE, txOut = TRUE, ercc = TRUE))

test_that("Demo data works correctly", {
    txi <- valid_txi
    expect_true(validate_txi(txi))
})

test_that("Valid txOut works correctly ", {
     txi <- valid_txi_ercc
     expect_true(validate_txi(txi))
     txi <- valid_txi_tx
     expect_true(validate_txi(txi))
     txi <- valid_txi_tx_ercc
     expect_true(validate_txi(txi))
})

test_that("Invalid txi class throws error", {
    txi <- "a"
    msg <- 'is(txi, "list") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Missing matrices/df throws error", {
    msg <- 'all(c("counts", "abundance", "length", "anno") %in% names(txi)) is not TRUE'
    # counts
    txi <- valid_txi
    txi$counts <- NULL
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # abundance
    txi <- valid_txi
    txi$abundance <- NULL
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # length
    txi <- valid_txi
    txi$length <- NULL
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # anno
    txi <- valid_txi
    txi$anno <- NULL
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # empty
    txi <- list()
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Invalid matrices class throws error", {
    produce_invalid_matrix <- function(fill, m) {
        res <- matrix(rep(fill, ncol(m) * nrow(m)), ncol = ncol(m))
        rownames(res) <- rownames(m)
        colnames(res) <- colnames(m)
        res
    }
    # counts
    txi <- valid_txi
    txi$counts <- produce_invalid_matrix("a", txi$counts)
    msg <- 'is.numeric(txi$counts) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # abundance
    txi <- valid_txi
    txi$abundance <- produce_invalid_matrix("a", txi$counts)
    msg <- 'is.numeric(txi$abundance) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # length
    txi <- valid_txi
    txi$length <- produce_invalid_matrix("a", txi$counts)
    msg <- 'is.numeric(txi$length) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # ruvg
    txi <- valid_txi
    txi$ruvg_counts <- produce_invalid_matrix("a", txi$counts)
    msg <- 'is.numeric(txi$ruvg_counts) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # combat
    txi <- valid_txi
    txi$combat_counts <- produce_invalid_matrix("a", txi$counts)
    msg <- 'is.numeric(txi$combat_counts) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Invalid counts class throws error", {
    # counts
    txi <- valid_txi
    txi$counts <- "a"
    msg <- 'is(txi$counts, "matrix") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # abundance
    txi <- valid_txi
    txi$abundance <- "a"
    msg <- 'is(txi$abundance, "matrix") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # length
    txi <- valid_txi
    txi$length <- "a"
    msg <- 'is(txi$length, "matrix") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # ruvg_counts
    txi <- valid_txi
    txi$ruvg_counts <- "a"
    msg <- 'is(txi$ruvg_counts, "matrix") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # combat_counts
    txi <- valid_txi
    txi$combat_counts <- "a"
    msg <- 'is(txi$combat_counts, "matrix") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Different matrices colnames throw error", {
    # counts
    txi <- valid_txi
    colnames(txi$counts)[1] <- "x"
    msg <- 'identical(colnames(txi$counts), colnames(txi$abundance)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # abundance
    txi <- valid_txi
    colnames(txi$abundance)[1] <- "x"
    msg <- 'identical(colnames(txi$counts), colnames(txi$abundance)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # length
    txi <- valid_txi
    colnames(txi$length)[1] <- "x"
    msg <- 'identical(colnames(txi$counts), colnames(txi$length)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # ruvg_counts
    txi <- valid_txi
    txi$ruvg_counts <- txi$counts
    colnames(txi$ruvg_counts)[1] <- "x"
    msg <- 'identical(colnames(txi$counts), colnames(txi$ruvg_counts)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # combat_counts
    txi <- valid_txi
    txi$combat_counts <- txi$counts
    colnames(txi$combat_counts)[1] <- "x"
    msg <- 'identical(colnames(txi$counts), colnames(txi$combat_counts)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Invalid anno throws error", {
    # Invalid class
    txi <- valid_txi
    txi$anno <- "a"
    msg <- 'is(txi$anno, "data.frame") is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # Missing columns
    ## id
    txi <- valid_txi
    txi$anno$id <- NULL
    msg <- 'expected_col %in% colnames(txi$anno) are not all TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    ## ensembl_gene
    txi <- valid_txi
    txi$anno$ensembl_gene <- NULL
    msg <- 'expected_col %in% colnames(txi$anno) are not all TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    ## symbol
    txi <- valid_txi
    txi$anno$symbol <- NULL
    msg <- 'expected_col %in% colnames(txi$anno) are not all TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    ## entrez_id
    txi <- valid_txi
    txi$anno$entrez_id <- NULL
    msg <- 'expected_col %in% colnames(txi$anno) are not all TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    ## transcript_type
    txi <- valid_txi
    txi$anno$transcript_type <- NULL
    msg <- 'expected_col %in% colnames(txi$anno) are not all TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
    # Different id
    txi <- valid_txi
    txi$anno$id <- paste0(txi$anno$id, ".1")
    msg <- 'identical(txi$anno$id, rownames(txi$counts)) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})

test_that("Invalid txOut throws error", {
    # txOut is FALSE
    txi <- valid_txi
    txi$anno$ensembl_gene <- paste0(txi$anno$id, "a")
    msg <- 'all(anno_not_ercc$id == anno_not_ercc$ensembl_gene) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)

    txi <- valid_txi_ercc
    txi$anno$ensembl_gene <- paste0(txi$anno$id, "a")
    msg <- 'all(anno_not_ercc$id == anno_not_ercc$ensembl_gene) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)

    # txOut is TRUE
    txi <- valid_txi_tx
    txi$anno$ensembl_gene <- txi$anno$id
    msg <- 'all(anno_not_ercc$id != anno_not_ercc$ensembl_gene) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)

    txi <- valid_txi_tx_ercc
    txi$anno$ensembl_gene <- txi$anno$id
    msg <- 'all(anno_not_ercc$id != anno_not_ercc$ensembl_gene) is not TRUE'
    expect_error(validate_txi(txi), msg, fixed = TRUE)
})
