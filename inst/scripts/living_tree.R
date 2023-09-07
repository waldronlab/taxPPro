## This script is for the living tree project (SILVA)
##
## This is the url of the file (active as of Aug 30, 2023)
## I also backed up the file in this package (to make things faster).
## url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'

library(ape)
library(taxonomizr)
library(dplyr)
library(purrr)


## I need to get taxid information from the accession IDs with taxononmizr.
## This is not enough for all tips, so I also need to get taxids from the organisms names with taxizedb.
## Still, I need to get some missing taxids with taxize.

tree <- read.tree('https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree')
tip_labels <- tree$tip.label
sql <- '~/accessionTaxa.sql' # this takes a while - date Aug 29, 2023 between 70 and 80 GB - Got it with taxonomizr

acc <- sub("^'([^,]+).*", "\\1", tip_labels)
tax_names <- sub('^.+, (.+), .+, .+, .+$'  ,'\\1', tip_labels) |>
    {\(y) gsub('"', '', y)}()

taxids <- accessionToTaxa(accessions = acc, sqlFile = sql, version = 'base')
pos <- which(is.na(taxids))
missing_taxids <- taxizedb::name2taxid(tax_names[pos], db = 'ncbi')
taxids[pos] <- missing_taxids

taxids_ranks <- taxizedb::taxid2rank(taxids, db = 'ncbi')
pos2 <- which(is.na(taxids_ranks))
taxids[pos2] <- as.character(taxize::get_uid(tax_names[pos2], db = 'ncbi')) # luckily, tax names were unique
taxids_ranks[pos2] <- flatten_chr(taxize::tax_rank(taxids[pos2], db = 'ncbi'))

tree_data <- data.frame(
    tip_label = tip_labels,
    accesion = acc,
    taxid = unname(taxids),
    taxname = tax_names,
    taxrank = taxids_ranks

)

new_tree_data <- tree_data |>
    filter(taxrank == 'species') |> ## only species in the tree
    arrange(taxid, tip_label) |>
    slice_head(n = 1, by = 'taxid') ## Remove duplicated taxids

new_tree <- ape:::keep.tip(tree, tip = new_tree_data$tip_label)
new_tree_data <- new_tree_data[match(new_tree$tip.label, new_tree_data$tip_label),]

# all(new_tree_data$tip_label == new_tree$tip.label)

new_tree$tip.label <- new_tree_data$taxid
new_tree_data <- new_tree_data |>
    mutate(
        old_tip_label = tip_label,
        tip_label = taxid
    )

taxonomy <- taxizedb::classification(new_tree_data$tip_label, db = 'ncbi')
new_taxonomy <- map(taxonomy, ~ {
    ranks_ <- c(
        'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
    )
    .x |>
        filter(rank %in% ranks_) |>
        select(rank, id) |>
        tidyr::pivot_wider(
            names_from = 'rank', values_from = 'id'
        ) |>
        {\(y) set_names(y, sub("$", "_taxid", colnames(y)))}()
}) |>
    bind_rows(.id = 'original_taxid') |>
    relocate(
        original_taxid, kingdom_taxid = superkingdom_taxid, phylum_taxid, class_taxid,
        order_taxid, family_taxid, genus_taxid, species_taxid
    )

new_tree_data <- left_join(
    new_tree_data, new_taxonomy, by = c('tip_label' = 'original_taxid')
)


# export data -------------------------------------------------------------
tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
write.table(
    x = new_tree_data, file = tree_data_fname, sep = '\t',
    quote = TRUE, row.names = FALSE
)

tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')
write.tree(phy = new_tree, file = tree_fname)



# Add labels to internal nodes --------------------------------------------
tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')

ltree <- read.tree(tree_fname)
ldata <- read.table(tree_data_fname, sep = '\t', header = TRUE, row.names = NULL)

getMRCA <- function(t, df) {
    res <- phytools::findMRCA(t, tips = df[['tip_label']])
    if (is.null(res))
        res <- NA
    res
}

















