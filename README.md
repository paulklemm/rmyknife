# rmyknife

<!-- TOC depthFrom:2 -->

- [💾 Installation](#💾-installation)
- [📦 Reproducible project environments](#📦-reproducible-project-environments)
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

## 📦 Reproducible project environments

Three layers pin an analysis project: the singularity image (OS, R, system libraries), the R package library (`renv.lock` plus dated repositories), and the code (git).
`renv` already covers the middle layer.
These functions add the container layer as a checksummed, in-repo fact recorded in `environment.lock`, and make the whole thing verifiable and archivable.

Run everything **inside the project's image**, because `renv` builds the library with the running R.

```r
# Set the project up. Safe on existing projects: an existing renv.lock and
# .Rprofile are left alone, only the repository pins are guaranteed.
rmyknife::project_init()

# A site that needs a bind mount, and a project that uses a tool container:
rmyknife::project_init(bind = "/data:/data", aux_images = "/path/to/tool.sif")

# Check that everything is present, consistent and restorable
rmyknife::project_verify()

# Archive the compute environment: the image and the built package library
rmyknife::project_backup()

# Bring one back
rmyknife::project_restore("backup/myproject_2026-07-24_a1b2c3d.tar.zst", destination = "restored")
```

`project_backup()` archives the **compute environment only** — the singularity images and the built renv library, with `renv.lock` to describe it.
The project code, its history and its `.Rprofile` are not included, because they live in git.
The archive records the commit the environment served, so the two halves can be paired back up.

`project_verify()` reports two independent kinds of restorability.
**From backup** is offline and exact: the archived image plus the binary library, so no compilation and no network.
**From lockfile** is a rebuild from scratch and needs every package to resolve, which is commonly partial for a project converted from a pre-renv state — a locally installed package will never restore over the network.
A project can be perfectly restorable from its backup while its lockfile is still messy, so the useful order is convert, back up immediately, then clean the lockfile up at leisure.

For a project that has no `renv` yet, the CRAN snapshot defaults to the **image build date** rather than today, since its packages are frozen at image build time.

Launching through a `latest/` symlink is fine.
`environment.lock` always records the versioned image the symlink resolves to, never the moving pointer, and `project_verify()` resolves symlinks before comparing.
So once `latest/` moves on to a newer image, verify tells you that you are no longer running the container the project was pinned to.

Note that `project_backup()` copies the images, so archives are large — expect a few GB for a typical R image plus its library.
Use `include = "library"` for a quick snapshot without them.

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

- _2026-07-24_
  - Add `project_init`, `project_verify`, `project_backup` and `project_restore` for reproducible project environments
  - Bump to `0.4.0`
- _2026-04-20_
  - Code-review package
  - Bump to `0.3.6`
- _2026-04-02_
  - Fix error when pruning using `make`
  - Bump to `0.3.5`
- _2026-02-10_
  - `ensembl_to_entrez` will now use biomaRt if supplied
  - Bump to 0.3.4
- _2026-01-21_
  - `make` lists the number of targets to be made as well as the total number
  - Bump to 0.3.3
- _2025-08-26_
  - `make` now loads successfully built targets even if one of them fails
- _2025-08-20_
  - Expose extensions in datatable function
  - Bump to 0.3.1
- _2024-09-26_
  - Allow `make` to check if all variables in local environment are up to date
  - Bump to 0.3.0
- *2023-08-04*
  - Add workers to `make`
  - Bump to 0.2.9
- *2022-10-07*
  - Enhance `make`
  - Bump to 0.2.8
- *2022-10-06*
  - Add `make` for making and loading targets
  - Bump to 0.2.7
- *2022-09-28*
  - Add `get_gzipped_stream`, `get_uniprot` and `get_uniprot_with_ensembl`
  - Bump to 0.2.6
- *2022-08-15*
  - Add `highlight` parameter to `plot_volcano`. Bump to 0.2.5
- *2022-07-22*
  - Add `plot_volcano`. Bump to 0.2.4
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
