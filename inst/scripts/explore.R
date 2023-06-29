library(taxPPro)
library(purrr)
library(dplyr)
library(bugphyzz)

# ncbi_taxonomy <- get_ncbi_taxonomy(force_download = TRUE)
fname <- system.file(
    'extdata/proc_include_all_06292023.txt', package = 'taxPPro',
    mustWork = TRUE
)
taxids <- read.table(fname, header = FALSE)[[1]]

phys <- physiologies(remove_false = TRUE, full_source = FALSE)

bp_ids <- map(phys, ~ .x$NCBI_ID) |>
    unlist(recursive = TRUE, use.names = FALSE) |>
    unique() |>
    {\(y) y[!is.na(y)] }() |>
    {\(y) y[y != 'unknown']}()


not_ids <- bp_ids[which(!bp_ids %in% taxids)]
class_not_ids <- taxizedb::classification(not_ids, db = 'ncbi')
class_not_ids_na <- keep(class_not_ids, ~ any(is.na(.x)))
tz_output <- taxize::classification(names(class_not_ids_na), db = 'ncbi')
tz_na <- keep(tz_output, ~ any(is.na(.x)))


x <- taxizedb::taxid2rank(bp_ids, db = 'ncbi')

pos <- which(!is.na(x))

df <- data.frame(
    ids = bp_ids[pos],
    rank = x[pos]
)

count(df, rank)


ex_all_fname <- name <- system.file(
    'extdata/proc_exclude_unclassified_uncultured_06292023.txt', package = 'taxPPro',
    mustWork = TRUE
)
exclude_all <- read.table(ex_all_fname, header = FALSE)[[1]]


mean(bp_ids %in% exclude_all)


ids <- split(df, factor(df$rank))
ids <- map(ids, ~ .x$ids)


map(ids, ~ mean(.x %in% taxids))
map_int(ids, length)




sp <- ids$species
output <- taxizedb::children(sp, db = 'ncbi')

strain_ids <- map(output, ~ {
    if (!nrow(.x)) {
        return(NA)
    } else {
        res <- filter(.x, rank == 'strain') |>
            pull(id)
        return(res)
    }
})
strain_ids <- unlist(unique(strain_ids))
strain_ids <- strain_ids[!is.na(strain_ids)]


mean(strain_ids %in% taxids)






ranks2 <- taxizedb::taxid2rank(taxids, db = 'ncbi')

df2 <- data.frame(
    id = taxids,
    rank = ranks2
)


ids2 <- split(df2, factor(df2$rank))
new_sp_ids <- ids2$species$id

children <- taxizedb::children(new_sp_ids, db = 'ncbi')


strains <- map(children, ~ {
    if (!nrow(.x)) {
        return(NA)
    } else {
        res <- filter(.x, rank == 'strain') |>
            pull(id)
        return(res)
    }
})
strains <- unlist(unique(strains))
strains <- strain_ids[!is.na(strains)]



strains %in% taxids


ids2$strain$id %in% strains



ranks <- c(
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
)

# taxizedb::db_download_ncbi(verbose = TRUE, overwrite = TRUE)

output <- taxizedb::classification(taxids, db = 'ncbi')
output_na <- keep(output, ~any(is.na(.x)))

id_ranks <- taxizedb::taxid2rank(taxids, db = 'ncbi')

df <- data.frame(ids = as.character(taxids), rank = id_ranks)

myFun <- function(x) {
    df <- dplyr::filter(x, rank %in% ranks)
}

sp_ids <- filter(df, rank == 'species')
sp <- output[sp_ids$ids]
filtered <- map(sp, myFun)

gn_ids <- map_chr(filtered, ~ {
    filter(.x, rank == 'genus') |>
        pull(id)
})

