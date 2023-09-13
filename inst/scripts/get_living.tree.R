## This script is for the living tree project (SILVA)
##
## This is the url of the file (active as of Aug 30, 2023)
## I also backed up the file in this package (to make things faster).
## url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'


# coreurl <- 'https://imedea.uib-csic.es/mmg/ltp/?smd_process_download=1&download_id=459'
# allurl <- 'https://imedea.uib-csic.es/mmg/ltp/?smd_process_download=1&download_id=458'


library(ape)
library(taxonomizr)
library(dplyr)
library(purrr)
library(tidyr)


sql <- '~/accessionTaxa.sql'
# tree <- read.tree('https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree')

# coreurl <- 'https://imedea.uib-csic.es/mmg/ltp/?smd_process_download=1&download_id=459'
allurl <- 'https://imedea.uib-csic.es/mmg/ltp/?smd_process_download=1&download_id=458'
tree <- read.tree(allurl)


tip_labels <- tree$tip.label
accessions <- sub("^'([^,]+).*", "\\1", tip_labels)
# taxnames <- sub('^.+, (.+), .+, .+, .+$'  ,'\\1', tip_labels) |>
#     {\(y) gsub('"', '', y)}()

taxnames <- sub("^[^,]*,([^,]*),.*$",'\\1', tip_labels) |>
    {\(y) gsub('"', '', y)}() |>
    stringr::str_squish()

taxids <- accessionToTaxa(
    accessions = accessions, sqlFile = sql, version = 'base'
)

missing_taxa <- taxnames[which(is.na(taxids))]
missing_taxa[which(missing_taxa == 'Micromonospora okii')] <- 'Micromonospora sp. TP-A0468'
missing_taxa[which(missing_taxa == 'Sala cibi')] <- 'Salella cibi'
not_missing_anymore_taxa <- taxizedb::name2taxid(missing_taxa, db = 'ncbi')

taxids[which(is.na(taxids))] <- not_missing_anymore_taxa

taxonomy <- taxizedb::classification(unique(taxids), db = 'ncbi')
missing_taxonomy_positions <- which(map_lgl(taxonomy, ~ all(is.na(.x))))
not_missing_anymore_taxonomy <- taxize::classification(
    names(missing_taxonomy_positions), db =  'ncbi'
)
chr_vct <- map_chr(not_missing_anymore_taxonomy, ~ tail(.x$id, 1))

for (i in seq_along(chr_vct)) {
    pos <- which(taxids == names(chr_vct[i]))
    message(
        'Replacing ', unique(taxids[pos]), ' with ', chr_vct[i]
    )
    taxids[pos] <- chr_vct[i]
}

chr_vct2 <- chr_vct[!chr_vct %in% names(taxonomy)]
not_missing_anymore_taxonomy <- not_missing_anymore_taxonomy[names(chr_vct2)]
taxonomy[names(chr_vct2)] <- not_missing_anymore_taxonomy
names(taxonomy[names(chr_vct2)]) <- chr_vct2
taxonomy <- taxonomy[which(!map_lgl(taxonomy, ~ all(is.na(.x))))]

taxonomic_ranks <- c(
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'
)

new_taxonomy <- map(taxonomy, ~ {
    x <- .x
    x <- x[which(x$rank %in% taxonomic_ranks),]
    m <- matrix(x$id, byrow = TRUE, nrow = 1)
    colnames(m) <- x$rank
    d <- as.data.frame(m)
    colnames(d) <- sub("$", "_taxid", colnames(d))
    d
}) |>
    bind_rows(.id = 'taxid') |>
    relocate(
        taxid, kingdom_taxid = superkingdom_taxid, phylum_taxid, class_taxid,
        order_taxid, family_taxid, genus_taxid, species_taxid
    )

tree_data <- data.frame(
    tip_label = tip_labels,
    accession = accessions,
    taxid = taxids,
    taxname = taxnames
) |>
    left_join(new_taxonomy, by = 'taxid')

getMRCA <- function(tree, tips) {
    res <- phytools::findMRCA(tree = tree, tips = tips)
    if (is.null(res))
        res <- NA
    res
}

tx <- paste0(taxonomic_ranks, '_taxid')
tx[which(tx == 'superkingdom_taxid')] <- 'kingdom_taxid'

mrcas <- flatten(map(tx, ~ split(tree_data, factor(tree_data[[.x]]))))
mrcas <-map(mrcas, ~ .x[['tip_label']])
mrcas <- map_int(mrcas, ~ getMRCA(tree, .x))
mrcas <- mrcas[!is.na(mrcas)]
mrcas_df <- data.frame(node = unname(mrcas), node_label = names(mrcas))
mrcas_df <- mrcas_df |>
    group_by(node) |>
    mutate(n_labels = length(unique(node_label))) |>
    mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    ungroup() |>
    distinct()

nodes <- data.frame(node = length(tree$tip.label) + 1:tree$Nnode) |>
    left_join(mrcas_df, by = 'node') |>
    mutate(
        node_label = ifelse(is.na(node_label), paste0('n', node), node_label),
        n_labels = ifelse(is.na(n_labels), 0, n_labels),
    )

tree$node.label <- nodes$node_label

nodes_with_taxid <- nodes |>
    filter(!grepl('^n', node_label)) |>
    separate_longer_delim(node_label, delim = '+') |>
    rename(taxid = node_label) |>
    group_by(node) |>
    mutate(node_label = paste0(taxid, collapse = '+')) |>
    ungroup()

all(nodes_with_taxid$node_label %in% tree$node.label)

nodes_taxonomy <- taxizedb::classification(
    x = unique(nodes_with_taxid$taxid), db = 'ncbi'
)
nodes_new_taxonomy <- map(nodes_taxonomy, ~ {
    x <- .x
    x <- x[which(x$rank %in% taxonomic_ranks),]
    m <- matrix(x$id, byrow = TRUE, nrow = 1)
    colnames(m) <- x$rank
    d <- as.data.frame(m)
    colnames(d) <- sub("$", "_taxid", colnames(d))
    d
}) |>
    bind_rows(.id = 'taxid') |>
    relocate(
        taxid, kingdom_taxid = superkingdom_taxid, phylum_taxid, class_taxid,
        order_taxid, family_taxid, genus_taxid, species_taxid
    ) |>
    discard(~ all(is.na(.x)))

node_data <- left_join(nodes_with_taxid, nodes_new_taxonomy, by= 'taxid')



# Export data -------------------------------------------------------------
tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')
ape::write.tree(tree, tree_fname)

tree_data_fname <- file.path('inst', 'extdata', 'livingTree_tips.tsv')
write.table(
    tree_data, tree_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

node_data_fname <- file.path('inst', 'extdata', 'livingTree_nodes.tsv')
write.table(
    node_data, node_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

# taxids_ranks <- taxizedb::taxid2rank(taxids, db = 'ncbi')
# pos2 <- which(is.na(taxids_ranks))
# taxids[pos2] <- as.character(taxize::get_uid(tax_names[pos2], db = 'ncbi')) # luckily, tax names were unique
# taxids_ranks[pos2] <- flatten_chr(taxize::tax_rank(taxids[pos2], db = 'ncbi'))
#
# tree_data <- data.frame(
#     tip_label = tip_labels,
#     accesion = acc,
#     taxid = unname(taxids),
#     taxname = tax_names,
#     taxrank = taxids_ranks
#
# )

# new_tree_data <- tree_data |>
#     filter(taxrank == 'species') |> ## only species in the tree
#     arrange(taxid, tip_label) |>
#     slice_head(n = 1, by = 'taxid') ## Remove duplicated taxids
#
# new_tree <- ape:::keep.tip(tree, tip = new_tree_data$tip_label)
# new_tree_data <- new_tree_data[match(new_tree$tip.label, new_tree_data$tip_label),]
#
# # all(new_tree_data$tip_label == new_tree$tip.label)
#
# new_tree$tip.label <- new_tree_data$taxid
# new_tree_data <- new_tree_data |>
#     mutate(
#         old_tip_label = tip_label,
#         tip_label = taxid
#     )
#
# taxonomy <- taxizedb::classification(new_tree_data$tip_label, db = 'ncbi')
# new_taxonomy <- map(taxonomy, ~ {
#     ranks_ <- c(
#         'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
#     )
#     .x |>
#         filter(rank %in% ranks_) |>
#         select(rank, id) |>
#         tidyr::pivot_wider(
#             names_from = 'rank', values_from = 'id'
#         ) |>
#         {\(y) set_names(y, sub("$", "_taxid", colnames(y)))}()
# }) |>
#     bind_rows(.id = 'original_taxid') |>
#     relocate(
#         original_taxid, kingdom_taxid = superkingdom_taxid, phylum_taxid, class_taxid,
#         order_taxid, family_taxid, genus_taxid, species_taxid
#     )
#
# new_tree_data <- left_join(
#     new_tree_data, new_taxonomy, by = c('tip_label' = 'original_taxid')
# )


# export data -------------------------------------------------------------
# tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
# write.table(
#     x = new_tree_data, file = tree_data_fname, sep = '\t',
#     quote = TRUE, row.names = FALSE
# )

# tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')
# write.tree(phy = new_tree, file = tree_fname)



# Add labels to internal nodes --------------------------------------------
# tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
# tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')

# ltree <- read.tree(tree_fname)
# ldata <- read.table(tree_data_fname, sep = '\t', header = TRUE, row.names = NULL)
# ldata <- modify(ldata, as.character)

# getMRCA <- function(t, df) {
#     res <- phytools::findMRCA(t, tips = df[['tip_label']])
#     if (is.null(res))
#         res <- NA
#     res
# }

# taxRanks <- c(
#     'kingdom', 'phylum', 'class', 'order', 'family', 'genus'
# ) |>
#     {\(y) sub("$", "_taxid", y)}()

# mrcas <- map(taxRanks, ~ {
#     splitted_df <- split(ldata, factor(ldata[[.x]]))
#     output <- map_chr(splitted_df, \(y) as.character(getMRCA(t = ltree, df = y)))
#     return(output)
# })
# names(mrcas) <- taxRanks
# mrcas <- map(mrcas, ~ {
#     v <- .x[!is.na(.x)]
#     df <- data.frame(names = names(v), nodes = unname(v))
#     df |>
#         group_by(nodes) |>
#         mutate(names = paste0(names, collapse = '+')) |>
#         ungroup() |>
#         distinct()
#
# })
#
# mrcas_df <- bind_rows(mrcas) |>
#         group_by(nodes) |>
#         mutate(names = paste0(names, collapse = '+')) |>
#         ungroup() |>
#         distinct()
# original_labels_df <- data.frame(
#     nodes = as.character(length(ltree$tip.label) + 1:ltree$Nnode),
#     names_original = ltree$node.label
# )
#
# new_labels_df <- left_join(original_labels_df, mrcas_df, by = 'nodes')
# new_node_labels <- ifelse(is.na(new_labels_df$names), "", new_labels_df$names)
#
# ltree$node.label <- new_node_labels



# Export data again -------------------------------------------------------
#
# tree_data_fname <- file.path('inst', 'extdata', 'livingTree.tsv')
# write.table(
#     x = ldata, file = tree_data_fname, sep = '\t',
#     quote = TRUE, row.names = FALSE
# )
#
# tree_fname <- file.path('inst', 'extdata', 'livingTree.newick')
# write.tree(phy = ltree, file = tree_fname)




## just some code to check that I get the same using MRCA
# ltree$node.label
# original_labels_df <- data.frame(names_original = ltree$node.label, nodes = as.character(length(ltree$tip.label) + 1:ltree$Nnode))
# testdf <- left_join(mrcas_df, original_labels_df, by = 'nodes')
#
# vec <- as.integer(testdf$names)
# vec <- vec[!is.na(vec)]
# vec_names <- taxizedb::taxid2name(vec, db = 'ncbi')
# test2 <- data.frame(names = as.character(vec), new_names = vec_names)
#
# test3 <- left_join(testdf, test2, by = 'names')





