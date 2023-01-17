library(rnaseq)

suppressMessages(txi <- get_demo_txi())
txi_ruv <- ruvg_normalization(txi, housekeeping_genes = c("RPL9", "RPL24"))

validate_table <- function(df, table_name) {
    stopifnot(table_name %in% c("counts", "abundance", "fpkm", "ruvg_counts"))

    m <- as.matrix(df[-1:-5])
    rownames(m) <- df$id
    if (table_name == "ruvg_counts") {
        current_txi <- txi_ruv
        expect_identical(rownames(m), rownames(current_txi[[table_name]]))
        expect_identical(m[,1], current_txi[[table_name]][,1])
        expect_identical(m[,2], current_txi[[table_name]][,2])
        expect_identical(m[,3], current_txi[[table_name]][,3])
        expect_identical(m[,4], current_txi[[table_name]][,4])
    } else {
        current_txi <- txi
        expect_identical(m, current_txi[[table_name]])
        expect_true(all(current_txi$anno$id %in% df$id))
        expect_true(all(current_txi$anno$ensembl_gene %in% df$ensembl_gene))
        expect_true(all(current_txi$anno$symbol %in% df$symbol))
        expect_true(all(current_txi$anno$entrez_id %in% df$entrez_id))
        expect_true(all(current_txi$anno$transcript_type %in% df$transcript_type))
    }
    expect_true(all(df$id %in% current_txi$anno$id))
    expect_true(all(df$ensembl_gene %in% current_txi$anno$ensembl_gene))
    expect_true(all(df$symbol %in% current_txi$anno$symbol))
    expect_true(all(df$entrez_id %in% current_txi$anno$entrez_id))
    expect_true(all(df$transcript_type %in% current_txi$anno$transcript_type))
}

test_that("get_raw_count_anno_df works correctly", {
    validate_table(get_raw_count_anno_df(txi), "counts")
})

test_that("get_tpm_anno_df works correctly", {
    validate_table(get_tpm_anno_df(txi), "abundance")
})

test_that("get_ruvg_anno_df works correctly", {
    validate_table(get_ruvg_anno_df(txi_ruv), "ruvg_counts")
})

test_that("get_anno_df works correctly with raw counts", {
    validate_table(get_anno_df(txi, "raw_counts"), "counts")
    validate_table(get_anno_df(txi, "counts"), "counts")
    validate_table(get_anno_df(txi, "raw_count"), "counts")
    validate_table(get_anno_df(txi, "count"), "counts")
})

test_that("get_anno_df works correctly with abundance", {
    validate_table(get_anno_df(txi, "abundance"), "abundance")
    validate_table(get_anno_df(txi, "tpm"), "abundance")
    validate_table(get_anno_df(txi, "TPM"), "abundance")
})

test_that("get_anno_df works correctly with ruvg", {
    validate_table(get_anno_df(txi_ruv, "ruvg"), "ruvg_counts")
    validate_table(get_anno_df(txi_ruv, "RUVg"), "ruvg_counts")
    validate_table(get_anno_df(txi_ruv, "RUV"), "ruvg_counts")
    validate_table(get_anno_df(txi_ruv, "ruv"), "ruvg_counts")
})

test_that("get_anno_df works correctly with fpkm", {
    validate_table(get_anno_df(txi, "fpkm"), "fpkm")
    validate_table(get_anno_df(txi, "FPKM"), "fpkm")
})

test_that("get_anno_df throws error with invalid table", {
    msg <- "colname %in% valid_colname is not TRUE"
    expect_error(get_anno_df(txi, "raw count"), msg)
    expect_error(get_anno_df(txi, "TPMs"), msg)
    expect_error(get_anno_df(txi, "Fpkm"), msg)
    expect_error(get_anno_df(txi, "ruvG"), msg)
})
