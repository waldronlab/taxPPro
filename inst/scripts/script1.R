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

data <- get_ncbi_taxids(keyword = 'all')


