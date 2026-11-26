# `taxPPro` — Taxonomic and Phylogenetic Propagation

This package contains functions and data, e.g. taxonomic trees, used for
performing filtering and propagation of bugphyzz annotations.

The functions in this package are used in the bugphyzzExports script.
Check the [`waldronlab/bugphyzzExports`](https://github.com/waldronlab/bugphyzzExports) repository on GitHub for details.

This repository is one step of the bugphyzz pipeline; the table under [How the data are produced](https://github.com/waldronlab/bugphyzz#how-the-data-are-produced) in the bugphyzz README shows how the repositories fit together.

## Installation

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install('waldronlab/taxPPro')
```

