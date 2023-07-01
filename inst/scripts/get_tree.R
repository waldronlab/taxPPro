library(taxPPro)
library(data.tree)
## The option 'bit' allows having more taxa, only excluding unclassified bacteria
tbl <- getNCBI(keyword = 'bit')
tree <- as.Node(tbl[, 'pathString', drop = FALSE], pathDelimiter = '|||')
tree$Do(addStrains)
tree_list <- as.list(tree)
usethis::use_data(tree_list, overwrite = TRUE)
## starts at 10:05 finished 11:12 (a little more than an hour)


library(purrr)
library(dplyr)
library(taxPPro)
library(data.tree)

ranks <- c(
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
    'strain'
)
ranks2 <- c('genus', 'species', 'strain')
fname <- system.file(
    'extdata/proc_exclude_unclassified_uncultured_06292023.txt',
    package = 'taxPPro', mustWork = TRUE
)
taxids <- read.table(fname, header = FALSE)[[1]]
taxranks <- taxizedb::taxid2rank(taxids, db = 'ncbi')
taxdf <- data.frame(id = taxids, rank = taxranks)
tax_separated <- map(split(taxdf, factor(taxdf$rank)), ~ .x$id)
gss <- tax_separated[ranks2]
classification <- map(gss, ~ taxizedb::classification(.x, db = 'ncbi'))

sp <- classification$species

sp_gn <- unique(unlist(map(sp, ~ filter(.x, rank == 'genus')$id)))
mean(sp_gn %in% gss$genus)

st <- classification$strain
st_sp <- unique(unlist(map(st, ~ filter(.x, rank == 'species')$id)))
mean(st_sp %in% gss$species)

st_gn <- unique(unlist(map(st, ~ filter(.x, rank == 'genus')$id)))
mean(st_gn %in% gss$genus)
