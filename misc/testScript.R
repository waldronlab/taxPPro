
library(bugphyzz)
library(purrr)

## Import bugphyzz
phys <- physiologies()
bugphyzz_taxids <- map(phys, ~ .x$NCBI_ID) |>
    unlist() |>
    unique() |>
    {\(y) y[!is.na(y)]}() |>
    as.character()

## Get NCBI table and ids
ncbi_tax <- get_ncbi_taxonomy(force_download = FALSE)
ncbi_taxids <- ncbi_tax$NCBI_ID

## All ids
proc_all_ids_fname <-
    system.file('extdata/proc_all_ids.txt', package = 'taxPPro')
proc_all_ids <- read.table(proc_all_ids_fname, header = FALSE)[[1]] |>
    as.character()

## Exclude all options form the NCBI
proc_exclude_all_fname <-
    system.file('extdata/proc_exclude_all.txt', package = 'taxPPro')
proc_exclude_all <- read.table(proc_exclude_all_fname, header = FALSE)[[1]] |>
    as.character()
x <- suppressWarnings(taxizedb::taxid2rank(proc_exclude_all))
table(x)[c('genus', 'species', 'strain')]

## Exclude informal names
proc_exclude_inf_fname <-
    system.file('extdata/proc_exclude_inf.txt', package = 'taxPPro')
proc_exclude_inf <- read.table(proc_exclude_inf_fname, header = FALSE)[[1]] |>
    as.character()

## Exclude unclassified
proc_exclude_uncla_fname <-
    system.file('extdata/proc_exclude_unclassified.txt', package = 'taxPPro')
proc_exclude_uncla <- read.table(proc_exclude_uncla_fname, header = FALSE)[[1]] |>
    as.character()




sum(bugphyzz_taxids %in% ncbi_taxids)
sum(bugphyzz_taxids %in% proc_all_ids)
sum(bugphyzz_taxids %in% proc_exclude_all)
sum(bugphyzz_taxids %in% proc_exclude_inf)
sum(bugphyzz_taxids %in% proc_exclude_uncla)


missing <- bugphyzz_taxids[!bugphyzz_taxids %in% proc_exclude_uncla]






phys |>
    map(~ dplyr::filter(.x,  NCBI_ID == '997346')) |>
    purrr::discard(~!nrow(.x))







