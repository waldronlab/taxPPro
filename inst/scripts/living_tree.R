## This script is for the living tree project
##
## This is the url of the file (active as of Aug 30, 2023)
## I also backed up the file in this package (to make things faster).
## url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'

library(ape)
library(taxonomizr)
library(dplyr)

tree <- read.tree('https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree')
tip_labels <- tree$tip.label
sql <- '~/accessionTaxa.sql' # this takes a while - date Aug 29, 2023
acc <- sub("^'([^,]+).*", "\\1", tip_labels)
taxids <- accessionToTaxa(accessions = acc, sqlFile = sql, version = 'base')
tax_names <- sub('^.+, (.+), .+, .+, .+$'  ,'\\1', tip_labels) |>
    {\(y) gsub('"', '', y)}()
pos <- which(is.na(taxids))
missing_taxids <- taxizedb::name2taxid(tax_names[pos], db = 'ncbi')
taxids[pos] <- missing_taxids

tree_data <- data.frame(
    tip_label = tip_labels,
    accesion = acc,
    taxid = unname(taxids),
    taxname = tax_names

)

new_tree_data <- tree_data |>
    arrange(taxid, tip_label) |>
    slice_head(n = 1, by = 'taxid')

new_tree <- ape:::keep.tip(tree, tip = new_tree_data$tip_label)
new_tree_data <- new_tree_data[match(new_tree$tip.label, new_tree_data$tip_label),]

# all(new_tree_data$tip_label == new_tree$tip.label)

new_tree$tip.label <- new_tree_data$taxid
new_tree_data <- new_tree_data |>
    mutate(
        old_tip_label = tip_label,
        tip_label = taxid
    )

tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
write.table(
    x = new_tree_data, file = tree_data_fname, sep = '\t',
    quote = TRUE, row.names = FALSE
)

tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')
write.tree(phy = new_tree, file = tree_fname)
