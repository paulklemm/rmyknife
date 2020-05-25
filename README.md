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
  rmyknife.use_biomart_mirror = TRUE,
  rmyknife.biomart_mirror_host = "useast.ensembl.org"
)

```

## ⏳ History

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
