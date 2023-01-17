library(rnaseq)

valid_txi <- suppressMessages(get_demo_txi(large = TRUE))

validate_gg <- function(gg, id_plot, meta, group) {
    if (meta & group) {
        expect_equal(nrow(gg$data), 4)
        expected_cols <- c("Dim1", "Dim2", "sample", "shape", "color", "group")
        expect_equal(colnames(gg$data), expected_cols)
        if (id_plot == "test_grp1") {
            expect_equal(gg$data$sample, c("a", "c", "e", "g"))
            expect_equal(gg$data$shape, paste0("shape", c(1,3,3,1)))
            expect_equal(gg$data$color, paste0("color", c(4,2,2,4)))
        }
        if (id_plot == "test_grp2") {
            expect_equal(gg$data$sample, c("b", "d", "f", "h"))
            expect_equal(gg$data$shape, paste0("shape", c(2,4,2,4)))
            expect_equal(gg$data$color, paste0("color", c(3,1,3,1)))
        }
    }
    if (meta & !group) {
        expect_equal(nrow(gg$data), 8)
        expected_cols <- c("Dim1", "Dim2", "sample", "shape", "color", "group")
        expect_equal(colnames(gg$data), expected_cols)
        expect_equal(gg$data$sample, letters[1:8])
        expect_equal(gg$data$shape, paste0("shape", c(1:4,3,2,1,4)))
        expect_equal(gg$data$color, paste0("color", c(4:1,2,3,4,1)))
        expect_equal(gg$data$group, paste0("group", rep(1:2, 4)))
    }
    if (!meta & !group) {
        expect_equal(nrow(gg$data), 8)
        expected_cols <- c("Dim1", "Dim2", "sample")
        expect_equal(colnames(gg$data), expected_cols)
        expect_equal(gg$data$sample, letters[1:8])
    }
}

test_that("Demo data works correctly: file pca infos", {
    pca_infos <- get_demo_pca_infos_file()
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    gg_list <- batch_pca(pca_infos, txi, metadata)
    expect_true(is(gg_list, "list"))
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = TRUE)
    validate_gg(gg_list[[2]], "test_grp2", meta = TRUE, group = TRUE)
})

test_that("Demo data works correctly: data.frame pca infos", {
    pca_infos <- get_demo_pca_infos_file() %>% read_csv(show_col_types = FALSE)
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = TRUE)
    validate_gg(gg_list[[2]], "test_grp2", meta = TRUE, group = TRUE)
})

test_that("No group in pca_infos works correctly", {
    pca_infos <- get_demo_pca_infos_file() %>%
        read_csv(show_col_types = FALSE) %>%
        mutate(group = NA, group_val = NA)
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
    validate_gg(gg_list[[2]], "test_grp2", meta = TRUE, group = FALSE)
})

test_that("Invalid pca_anno values throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'file.exists(pca_infos) is not TRUE'
    expect_error(batch_pca("a", txi, metadata), msg, fixed = TRUE)
    msg <- 'is(pca_infos, "data.frame") | is(pca_infos, "character") is not TRUE'
    expect_error(batch_pca(1, txi, metadata), msg, fixed = TRUE)
    msg <- 'nrow(pca_infos) > 0 is not TRUE'
    expect_error(batch_pca(data.frame(), txi, metadata), msg, fixed = TRUE)
})

test_that("Invalid txi values throws correct errors", {
    pca_infos <- get_demo_pca_infos_file() %>% read_csv(show_col_types = FALSE)
    metadata <- get_demo_metadata_file()
    msg <- 'is(txi, "list") is not TRUE'
    expect_error(batch_pca(pca_infos, "a", metadata), msg, fixed = TRUE)
    msg <- 'all(c("counts", "abundance", "length", "anno") %in% names(txi)) is not TRUE'
    expect_error(batch_pca(pca_infos, list(), metadata), msg, fixed = TRUE)
})

test_that("Invalid metadata values throws correct errors", {
    pca_infos <- get_demo_pca_infos_file() %>% read_csv(show_col_types = FALSE)
    txi <- valid_txi
    msg <- 'file.exists(metadata) is not TRUE'
    expect_error(batch_pca(pca_infos, txi, "a"), msg, fixed = TRUE)
    msg <- 'is(metadata, "data.frame") | is(metadata, "character") is not TRUE'
    expect_error(batch_pca(pca_infos, txi, 1), msg, fixed = TRUE)
    msg <- 'nrow(metadata) > 0 is not TRUE'
    expect_error(batch_pca(pca_infos, txi, data.frame()), msg, fixed = TRUE)
})

test_that("Invalid directories values throws correct errors", {
    pca_infos <- get_demo_pca_infos_file() %>% read_csv(show_col_types = FALSE)
    metadata <- get_demo_metadata_file()
    txi <- valid_txi
    msg <- 'is(outdir, "character") is not TRUE'
    expect_error(batch_pca(pca_infos, txi, metadata, outdir = 1), msg, fixed = TRUE)
    expect_error(batch_pca(pca_infos, txi, metadata, outdir = NA), msg, fixed = TRUE)
    msg <- 'is(r_objects, "character") is not TRUE'
    expect_error(batch_pca(pca_infos, txi, metadata, r_objects = 1), msg, fixed = TRUE)
    expect_error(batch_pca(pca_infos, txi, metadata, r_objects = NA), msg, fixed = TRUE)
})

test_that("Invalid force values throws correct errors", {
    pca_infos <- get_demo_pca_infos_file() %>% read_csv(show_col_types = FALSE)
    metadata <- get_demo_metadata_file()
    txi <- valid_txi
    msg <- 'is(force, "logical") is not TRUE'
    expect_error(batch_pca(pca_infos, txi, metadata, force = 1), msg, fixed = TRUE)
    msg <- '!is.na(force) is not TRUE'
    expect_error(batch_pca(pca_infos, txi, metadata, force = NA), msg, fixed = TRUE)
})

test_that("Valid pca_anno id_plot column works", {
    txi <- valid_txi
    pca_infos <- data.frame(id_plot = "a")
    gg_list <- batch_pca(pca_infos, txi)
    validate_gg(gg_list[[1]], "test_grp1", meta = FALSE, group = FALSE)
})

test_that("Valid pca_anno group, group_val and id_metadata columns work", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", group = "group",
                            group_val = "group1", id_metadata = "id")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = TRUE)
})

test_that("Valid pca_anno shape column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", shape = "color")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno color column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", color = "shape")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno size column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", size = 2)
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno title column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", title = "abc")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno legend.position column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.position = "right")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.position = "left")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.position = "top")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.position = "bottom")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno legend.box column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.box = "vertical")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "id", legend.box = "horizontal")
    gg_list <- batch_pca(pca_infos, txi, metadata)
    validate_gg(gg_list[[1]], "test_grp1", meta = TRUE, group = FALSE)
})

test_that("Valid pca_anno show_names column works", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    pca_infos <- data.frame(id_plot = "a", show_names = TRUE)
    gg_list <- batch_pca(pca_infos, txi)
    validate_gg(gg_list[[1]], "test_grp1", meta = FALSE, group = FALSE)
    pca_infos <- data.frame(id_plot = "a", show_names = FALSE)
    gg_list <- batch_pca(pca_infos, txi)
    validate_gg(gg_list[[1]], "test_grp1", meta = FALSE, group = FALSE)
})

test_that("Invalid pca_anno id_plot column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- '"id_plot" %in% colnames(pca_infos) is not TRUE'
    expect_error(batch_pca(data.frame(a = 1), txi, metadata), msg, fixed = TRUE)
    msg <- '!any(duplicated(pca_infos$id_plot)) is not TRUE'
    expect_error(batch_pca(data.frame(id_plot = c("a", "a")), txi, metadata), msg, fixed = TRUE)
    msg <- '!any(is.na(pca_infos$id_plot)) is not TRUE'
    expect_error(batch_pca(data.frame(id_plot = NA), txi, metadata), msg, fixed = TRUE)
    expect_error(batch_pca(data.frame(id_plot = c("a", NA)), txi, metadata), msg, fixed = TRUE)
})

test_that("Invalid pca_anno group columns content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", group = NA, group_val = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", group = "a", group_val = NA)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", group = "a", group_val = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    expect_error(suppressMessages(batch_pca(pca_infos, txi, metadata)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", group = "group", group_val = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi, metadata)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno use_normalisation column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", use_normalisation = 1)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", use_normalisation = NA)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", use_normalisation = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno min_counts column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", min_counts = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", min_counts = NA)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", min_counts = -1)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno id_metadata column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", id_metadata = 1, metadata)
    expect_error(suppressMessages(batch_pca(pca_infos, txi, metadata)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = NA)
    expect_error(suppressMessages(batch_pca(pca_infos, txi, metadata)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi, metadata)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", id_metadata = "group")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno size column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", size = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", size = NA)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", size = -1)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", size = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno shape column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", shape = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", shape = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno color column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", color = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", color = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno title column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", title = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno legend.position column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", legend.position = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", legend.position = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})

test_that("Invalid pca_anno legend.box column content throws correct errors", {
    txi <- valid_txi
    metadata <- get_demo_metadata_file()
    msg <- 'length(validate_pca_infos(pca_infos, metadata, txi)) == 0 is not TRUE'
    pca_infos <- data.frame(id_plot = "a", legend.box = 0)
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
    pca_infos <- data.frame(id_plot = "a", legend.box = "a")
    expect_error(suppressMessages(batch_pca(pca_infos, txi)), msg, fixed = TRUE)
})
