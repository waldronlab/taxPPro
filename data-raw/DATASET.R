## code to prepare `ncbi_direct_download` dataset goes here

usethis::use_data(ncbi_direct_download, overwrite = TRUE)

## Include all ids
proc_all_ids_fname <-
    system.file('extdata/proc_all_ids.txt', package = 'taxPPro')
proc_all_ids <- read.table(proc_all_ids_fname, header = FALSE)[[1]] |>
    as.character()

## Exclude all
proc_exclude_all_fname <-
    system.file('extdata/proc_exclude_all.txt', package = 'taxPPro')
proc_exclude_all <- read.table(proc_exclude_all_fname, header = FALSE)[[1]] |>
    as.character()

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

## Exclude uncultured


## Exclude Unclassified and Informal names


## Exclude Unclassifed and uncultured


## Exclude Informal names and uncultured








