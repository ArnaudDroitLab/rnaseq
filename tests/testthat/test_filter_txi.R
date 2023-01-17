library(rnaseq)

valid_txi <- suppressMessages(get_demo_txi(large = TRUE))

test_that("Valid filter_txi works", {
    # Some
    txi <- valid_txi
    txi_filter <- filter_txi(txi, c("a", "b"))
    expect_true(validate_txi(txi_filter))
    expect_equal(colnames(txi_filter$counts), c("a", "b"))

    # All
    txi <- valid_txi
    txi_filter <- filter_txi(txi, letters[1:8])
    expect_true(validate_txi(txi_filter))
    expect_equal(colnames(txi_filter$counts), letters[1:8])

    # None
    txi <- valid_txi
    txi_filter <- filter_txi(txi, character())
    expect_true(validate_txi(txi_filter))
    expect_equal(colnames(txi_filter$counts), NULL)
})

test_that("Invalid txi throws error", {
    # Invalid class
    txi <- "a"
    msg <- 'is(txi, "list") is not TRUE'
    expect_error(filter_txi(txi, c("a", "b")), msg, fixed = TRUE)
    # Other tests are done through validate_txi unit test suite
})

test_that("Invalid samples throws error", {
    txi <- valid_txi
    # One invalid sample
    msg <- 'all(samples %in% colnames(txi[[name]])) is not TRUE'
    expect_error(filter_txi(txi, c("a", "z")), msg, fixed = TRUE)

    # All invalid samples
    msg <- 'all(samples %in% colnames(txi[[name]])) is not TRUE'
    expect_error(filter_txi(txi, c("x", "z")), msg, fixed = TRUE)
})
