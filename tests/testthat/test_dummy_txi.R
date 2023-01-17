library(rnaseq)

txi <- get_demo_txi()
full_tables <- list()
full_tables$counts <- txi$counts
full_tables$abundance <- txi$abundance
full_tables$length <- txi$length
full_tables$anno <- txi$anno

txi_tx <- get_demo_txi(txOut = TRUE)
full_tables_tx <- list()
full_tables_tx$counts <- txi_tx$counts
full_tables_tx$abundance <- txi_tx$abundance
full_tables_tx$length <- txi_tx$length
full_tables_tx$anno <- txi_tx$anno

check_matrix_presence <- function(dummy_txi) {
    expect_true(!is.null(dummy_txi$counts))
    expect_true(!is.null(dummy_txi$abundance))
    expect_true(!is.null(dummy_txi$length))
    expect_true(!is.null(dummy_txi$anno))
}

test_that("Valid params works correctly with all matrices", {
    current_tables <- full_tables
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(!identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(!identical(dummy_txi$counts, dummy_txi$length))
    expect_true(!identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Valid params works correctly with all matrices except counts", {
    current_tables <- full_tables
    current_tables$counts <- NULL
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(!identical(dummy_txi$counts, dummy_txi$length))
    expect_true(!identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Valid params works correctly with all matrices except abundance", {
    current_tables <- full_tables
    current_tables$abundance <- NULL
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(!identical(dummy_txi$counts, dummy_txi$length))
    expect_true(!identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Valid params works correctly with all matrices except length", {
    current_tables <- full_tables
    current_tables$length <- NULL
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(!identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(identical(dummy_txi$counts, dummy_txi$length))
    expect_true(!identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Valid params works correctly with only counts matrix", {
    current_tables <- full_tables
    current_tables$abundance <- NULL
    current_tables$length <- NULL
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(identical(dummy_txi$counts, dummy_txi$length))
    expect_true(identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Valid params works correctly with only abundance matrix", {
    current_tables <- full_tables
    current_tables$abundance <- NULL
    current_tables$length <- NULL
    dummy_txi <- create_dummy_txi(current_tables, txOut = FALSE)
    check_matrix_presence(dummy_txi)
    expect_true(identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(identical(dummy_txi$counts, dummy_txi$length))
    expect_true(identical(dummy_txi$abundance, dummy_txi$length))
})


test_that("Valid params works correctly with all matrices and txOut is TRUE", {
    current_tables <- full_tables_tx
    dummy_txi <- create_dummy_txi(current_tables, txOut = TRUE)
    check_matrix_presence(dummy_txi)
    expect_true(!identical(dummy_txi$counts, dummy_txi$abundance))
    expect_true(!identical(dummy_txi$counts, dummy_txi$length))
    expect_true(!identical(dummy_txi$abundance, dummy_txi$length))
})

test_that("Invalid tables format throws correct error", {
    msg <- 'is(tables, "list") is not TRUE'
    expect_error(create_dummy_txi("a"), msg, fixed = TRUE)
    expect_error(create_dummy_txi(1), msg, fixed = TRUE)
})

test_that("Missing both counts and abundance tables throws correct error", {
    current_tables <- full_tables
    current_tables$counts <- NULL
    current_tables$abundance <- NULL
    msg <- 'any(c("counts", "abundance") %in% names(tables)) is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Missing anno throws correct error", {
    current_tables <- full_tables
    current_tables$anno <- NULL
    msg <- '"anno" %in% names(tables) is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Invalid anno format throws correct error", {
    current_tables <- full_tables
    current_tables$anno <- "a"
    msg <- 'is(tables$anno, "data.frame") is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
    current_tables$anno <- 1
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Invalid counts format throws correct error", {
    current_tables <- full_tables
    current_tables$counts <- "a"
    msg <- 'is(tables$counts, "matrix") is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
    current_tables$counts <- 1
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Invalid abundance format throws correct error", {
    current_tables <- full_tables
    current_tables$abundance <- "a"
    msg <- 'is(tables$abundance, "matrix") is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
    current_tables$abundance <- 1
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Invalid txOut format throws correct error", {
    current_tables <- full_tables
    msg <- 'is(txOut, "logical") is not TRUE'
    expect_error(create_dummy_txi(current_tables, txOut = "a"), msg, fixed = TRUE)
    expect_error(create_dummy_txi(current_tables, txOut = 1), msg, fixed = TRUE)

})
test_that("Identical counts and abundance tables throws correct error", {
    current_tables <- full_tables
    current_tables$abundance <- current_tables$counts
    msg <- '!identical(tables$counts, tables$abundance) is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Identical counts and length tables throws correct error", {
    current_tables <- full_tables
    current_tables$length <- current_tables$counts
    msg <- '!identical(tables$counts, tables$length) is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})

test_that("Identical abundance and length tables throws correct error", {
    current_tables <- full_tables
    current_tables$length <- current_tables$abundance
    msg <- '!identical(tables$abundance, tables$length) is not TRUE'
    expect_error(create_dummy_txi(current_tables), msg, fixed = TRUE)
})
