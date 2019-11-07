# rmyknife

The goal of rmyknife is to provide a centralized place for R functions that I often use but that do not belong to a more specialized package yet.
This package will probably be pretty messy.

## Installation

You can install the github version of rmyknife with:

```r
library(devtools)
devtools::install_github("paulklemm/rmyknife")
```

## Memoise for BiomaRt

For interfacing Biomart, we are using the [memoise](https://github.com/r-lib/memoise) package.
It will store all queries to a folder names `rmyknife_memoise` in your r temp folder.
You can set the following options to customize this behavior.

```r
options(rmyknife.use_memoise = FALSE)
options(rmyknife.memoise_path = "<some/other/path>")
```
