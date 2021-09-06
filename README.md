# rmyknife

<!-- TOC depthFrom:2 -->

- [💾 Installation](#💾-installation)
- [🧠 Memoise for BiomaRt](#🧠-memoise-for-biomart)
- [♻️ Using Ensembl mirror](#♻️-using-ensembl-mirror)
- [⏳ History](#⏳-history)

<!-- /TOC -->

The goal of rmyknife is to provide a centralized place for R functions that I often use but that do not belong to a more specialized package yet.
This package will probably be pretty messy.

## 💾 Installation

You can install the github version of rmyknife with:

```r
library(devtools)
devtools::install_github("paulklemm/rmyknife")
```

## 🧠 Memoise for BiomaRt

For interfacing Biomart, we are using the [memoise](https://github.com/r-lib/memoise) package.
It will store all queries to a folder names `rmyknife_memoise` in your r temp folder.
You can set the following options to customize this behavior.

```r
options(rmyknife.use_memoise = FALSE)
options(rmyknife.memoise_path = "<some/other/path>")
```

## ♻️ Using Ensembl mirror

Sometimes the Ensembl biomart server are unstable.
To still be able to conduct queries, we support alternative hosts which can be set as follows.

```r
# rmyknife.use_biomart_mirror = TRUE will use "useast.ensembl.org" as default
# You can specify a custom mirror using rmyknife.biomart_mirror_host
options(
  rmyknife.verbose = TRUE,
  rmyknife.use_biomart_mirror = TRUE,
  rmyknife.biomart_mirror_host = "useast.ensembl.org"
  # rmyknife.biomart_mirror_host = "uswest.ensembl.org"
  # rmyknife.biomart_mirror_host = "asia.ensembl.org"
)

```

## ⏳ History

- *2021-09-06*
  - Add `rnorvegicus` to `get_ensembl_dataset_from_version`. Bump to 0.2.3
- *2021-04-19*
  - Add `set_ggplot_defaults`. Bump to 0.2.2
- *2021-02-22*
  - Add option to get children of GO-terms
  - Bump version to 0.2.1
- *2021-02-19*
  - Get GO-name based on biomaRt, not GO.db
  - Allow custom biomart to be used in `attach_biomart`
  - Bump version to 0.2.0
- *2020-11-23*
  - Fix rare attach_biomart bug
  - Bump version to 0.1.10
- *2020-10-26*
  - Remove tidylog dependency
  - Bump version to 0.1.9
- *2020-08-05*
  - Add `print_params` function
  - Bump version to 0.1.7
- *2020-08-04*
  - Add `get_tempfile_from_url`
  - Bump version to 0.1.6
- *2020-07-27*
  - Add `attach_gene_symbol_from_entrez` and `ensembl_to_entrez` from mygo package
  - Bump version to 0.1.5
- *2020-07-02*
  - Fix minor deprecated calls to avoid warning messages from tibble and dplyr
  - Bump version to 0.1.4
- *2020-06-24*
  - Add `deseq2_result_to_tibble` function for converting deseq2 result objects to a dataframe
  - Bump version to 0.1.3
- *2020-05-25*
  - Add `write_xls` function
  - Add options `options(rmyknife.use_biomart_mirror = TRUE)`, `options(rmyknife.biomart_mirror_host = "useast.ensembl.org")` and `options(rmyknife.verbose = TRUE)`
  - Bump version to 0.1.2
  - Close [Remove hack while Ensembl is being transferred #1](https://github.com/paulklemm/rmyknife/issues/1)
- *2020-04-29*
  - Add `dt_datatable` function
  - Bump version to 0.1.0
- *2020-04-22*
  - Add function `get_gene_counts_per_cell`
  - Bump version to 0.0.5
- *2020-03-21*
  - Add [Remove hack while Ensembl is being transferred #1](https://github.com/paulklemm/rmyknife/issues/1)
  - Bump version to `0.0.4`
- *2020-02-04*
  - Added `get_genes_of_goterm_godb` function
- *2019-11-14*
  - Added `get_seurat_clustering` function
- *2019-11-13*
  - Added `read_cellranger_feature_bc_matrix` and `read_cellranger_feature_bc_matrix_h5` function for reading 10x genomics cellranger data
