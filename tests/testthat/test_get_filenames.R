library(rnaseq)

test_that("Correct params values works correctly", {
    # Correct dir_kallisto
    dir_quant <- get_demo_kallisto_dir()
    filenames <- get_filenames(dir_quant, file_extension = "tsv")
    expect_true(all(file.exists(filenames)))

    # Correct file_extension
    filenames <- get_filenames(dir_quant, file_extension = "tsv")
    expect_true(all(file.exists(filenames)))
})

test_that("Invalid dir_kallisto throw error", {
    msg <- 'dir.exists(dir_kallisto) is not TRUE'
    expect_error(get_filenames("a"), msg, fixed = TRUE)
})

test_that("Invalid file_extension throw error", {
    msg <- 'is.character(file_extension) is not TRUE'
    expect_error(get_filenames(get_demo_kallisto_dir(), 1), msg, fixed = TRUE)
    expect_error(get_filenames(get_demo_kallisto_dir(), NULL), msg, fixed = TRUE)
    expect_error(get_filenames(get_demo_kallisto_dir(), NA), msg, fixed = TRUE)

    msg <- 'length(file_extension) == 1 is not TRUE'
    expect_error(get_filenames(get_demo_kallisto_dir(), c("tsv", "a")), msg, fixed = TRUE)

    msg <- 'file_extension %in% c("h5", "tsv") is not TRUE'
    expect_error(get_filenames(get_demo_kallisto_dir(), "a"), msg, fixed = TRUE)
})
