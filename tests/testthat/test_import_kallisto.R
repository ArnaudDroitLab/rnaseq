library(rnaseq)

validate_txi <- function(txi) {
    expected_names <- c("abundance", "counts", "length", "countsFromAbundance",
                        "fpkm", "anno", "txOut")
    expect_equal(names(txi), expected_names)
    expect_equal(class(txi$counts), c("matrix", "array"))
    expect_equal(class(txi$abundance), c("matrix", "array"))
    expect_equal(class(txi$fpkm), c("matrix", "array"))
    expect_equal(class(txi$anno), "data.frame")
    expect_equal(class(txi$txOut), "logical")
    expect_equal(nrow(txi$counts), nrow(txi$abundance))
    expect_equal(ncol(txi$counts), ncol(txi$abundance))
    expect_equal(nrow(txi$counts), nrow(txi$fpkm))
    expect_equal(ncol(txi$counts), ncol(txi$fpkm))
    expect_equal(nrow(txi$counts), nrow(txi$anno))
    expect_identical(rownames(txi$counts), rownames(txi$abundance))
    expect_identical(rownames(txi$counts), rownames(txi$fpkm))
    expect_identical(rownames(txi$counts), txi$anno$id)
    expect_equal(colnames(txi$anno), c("id", "ensembl_gene", "symbol",
                                       "entrez_id", "transcript_type"))
}

test_that("Demo data works correctly", {
    abundances <- get_demo_abundance_files()
    file_anno <- get_demo_anno_file()
    suppressMessages(txi <- import_kallisto(abundances, file_anno))
    validate_txi(txi)
})

test_that("Custom anno expected transcript_type order", {
    anno_filename <- "extdata/valid_anno_multiple_transcript_type_expected_order.csv"
    anno_filename <- system.file(anno_filename, package = "rnaseq")
    dir_kallisto <- system.file("extdata/quant_test_valid", package = "rnaseq")
    filenames <- get_filenames(dir_kallisto, "tsv")
    suppressMessages(txi <- import_kallisto(filenames, anno = anno_filename))
    validate_txi(txi)
    expect_equal(txi$anno$id, c("gene1", "gene2", "gene3"))
    expect_equal(txi$anno$ensembl_gene, c("gene1", "gene2", "gene3"))
    expect_equal(txi$anno$symbol, c("g1", "g2", "g3"))
    expect_equal(txi$anno$entrez_id, c("123", "456", "789"))
    expect_equal(txi$anno$transcript_type, rep("protein_coding", 3))
})

test_that("Custom anno unexpected transcript_type order", {
    anno_filename <- "extdata/valid_anno_multiple_transcript_type_unexpected_order.csv"
    anno_filename <- system.file(anno_filename, package = "rnaseq")
    dir_kallisto <- system.file("extdata/quant_test_valid", package = "rnaseq")
    filenames <- get_filenames(dir_kallisto, "tsv")
    suppressMessages(txi <- import_kallisto(filenames, anno = anno_filename))
    validate_txi(txi)
    expect_equal(txi$anno$id, c("gene1", "gene2", "gene3"))
    expect_equal(txi$anno$ensembl_gene, c("gene1", "gene2", "gene3"))
    expect_equal(txi$anno$symbol, c("g1", "g2", "g3"))
    expect_equal(txi$anno$entrez_id, c("123", "456", "789"))
    expect_equal(txi$anno$transcript_type, rep("protein_coding", 3))
})

## TODO: test invalid cases (NA, etc...)
