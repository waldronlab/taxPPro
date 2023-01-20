library(dplyr)

## proc all is 'base' + 'informal' + 'unclassified' + 'uncultured'
## exclude all is 'base'
## exclude informal is 'base' + 'unclassified' + uncultured'
## exclude unclassified informal is 'base' + 'uncultured'
## exclude unclassified is 'base' + 'uncultured' +  'informal'
## exclude

options <- c(
    b = 'base',
    i = 'informal',
    l = 'unclassified',
    t = 'uncultured'
)

x <- 'blit'

stringr::str_split(x, pattern = '') |>
    unlist() |>
    sort() |>
    paste0(collapse = '')


fname <- system.file('extdata/proc_all_ids.txt', package = 'taxPPro')
taxids <- read.table(fname, header = FALSE)[[1]]
ncbi <- get_ncbi_taxonomy()

sum(taxids %in% ncbi$NCBI_ID)

ncbi_total <- ncbi[ncbi$NCBI_ID %in% taxids,]

## Number of genera per phylum
x <- ncbi_total |>
    filter(Rank == 'genus') |>
    count(kingdom, phylum) |>
    filter(kingdom == 'Archaea') |>
    arrange(-n)
head(x)

## Number of species per phylum
y <- ncbi_total |>
    filter(Rank == 'species') |>
    count(kingdom, phylum) |>
    filter(kingdom == 'Bacteria') |>
    arrange(-n)
head(y)


zzz <- c(
    bilt = 'proc_all_ids.txt',
    b = 'proc_exclude_all.txt',
    blt = 'proc_exclude_informal.txt',
    bt = 'proc_exclude_unclassified_informal.txt',
    bit = 'proc_exclude_unclassified.txt',
    bi = 'proc_exclude_unclassified_uncultured.txt',
    bl = 'proc_exclude_uncultured_informal.txt',
    bil = 'proc_exclude_uncultured.txt'
)

the_names <- names(.files())
output <- vector('list', length(the_names))
for (i in seq_along(output)) {
    message('Getting data for', the_names[i])
    output[[i]] <- get_ncbi_taxids(keyword = the_names[i])
    names(output)[i] <- the_names[i]
}

