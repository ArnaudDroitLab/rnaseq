# rnaseq

A package to facilitate RNA-Seq analysis using pre-computed annotations. Fresh version, old commits can still be viewed [here](https://github.com/CharlesJB/rnaseq)

## Install

In R:
 
```r
devtools::install_github("CharlesJB/rnaseq")
```

## Main functions:

* `import_kallisto`
* `get_raw_count_anno_df`
* `get_tpm_anno_df`
* `produce_pca`
* `deseq2_analysis`

Please see documentation using the `?` (i.e.: `?import_kallisto`)

## RUVg normalization

* `ruvg_normalization`

See also the `use_ruv` parameter of the `deseq2_analysis` function.
