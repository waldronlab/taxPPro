
## This code is for converting the tips of the metaphlan tree v3.1 from
## genome ids (GCA) to taxids

library(ape)
library(purrr)
library(dplyr)
library(tidyr)

column_names <- c(
    'genome_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
    'species'
)

taxonomic_ranks <- c(
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'
)

getMRCA <- function(tree, tips) {
    res <- phytools::findMRCA(tree = tree, tips = tips)
    if (is.null(res))
        res <- NA
    res
}

## It's faster to import the tree if it's stored in extdata instead of
## downloading from the url: http://cmprod1.cibio.unitn.it/biobakery3/mpa_v31_CHOCOPhlAn_201901_species_tree.nwk
## See also this post: https://forum.biobakery.org/t/announcing-the-metaphlan-3-1-phylogenetic-tree/5237
mpa_tree_fname <- system.file(
    'extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.nwk',
    package = 'taxPPro', mustWork = TRUE
)
mpa_tree <- ape::read.tree(mpa_tree_fname)

## I downloaded these data (assembly) from https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ on 08292023
## I should add some code here to save the file to a temporary location rather than using paths
assembly_data_current_fname <- '~/assembly_summary_genbank.txt'
assembly_data_historical_fname <- '~/assembly_summary_genbank_historical.txt'
assembly_data_current <- readr::read_tsv(
    assembly_data_current_fname, show_col_types = FALSE, skip = 1,
) |>
    {\(y) set_names(y, sub('#', '', colnames(y)))}() |>
    select(assembly_accession, taxid)
    # select(assembly_accession, taxid, species_taxid, organism_name)
assembly_data_historical <- readr::read_tsv(
    assembly_data_historical_fname, show_col_types = FALSE, skip = 1
) |>
    {\(y) set_names(y, sub('#', '', colnames(y)))}() |>
    select(assembly_accession, taxid)
    # select(assembly_accession, taxid, species_taxid, organism_name)
assembly_data <- bind_rows(assembly_data_current, assembly_data_historical) |>
    mutate(genome_id = sub('\\.\\d+', '', assembly_accession))
    # mutate(version = as.integer(sub('^.*\\.(\\d+)$', '\\1', assembly_accession)))

mpa_data <- data.frame(tip_label = mpa_tree$tip.label) |>
    separate(
        col = 'tip_label', into = column_names, sep = '\\|', remove = FALSE
    ) |>
    # modify(~ sub('^\\w__', '', .x)) |>
    left_join(assembly_data, by = 'genome_id') |>
    ## the next line of code filters only the most recent version of genome ID
    ## I could have done this in the code chunk above, but the dataset was
    ## too big (I think), so better do it here after the left join.
    slice_max(order_by = assembly_accession, n = 1, by = genome_id) |>
    select(tip_label, genome_id, taxid)
    # select(
    #     tip_label, genome_id, assembly_accession, taxid,
    #     original_species_taxid = species_taxid
    # ) |>
    # mutate(original_taxid = as.character(original_taxid))

## Update taxonomy information with the most recent one
taxonomy <- taxizedb::classification(x = unique(mpa_data$taxid), db = 'ncbi')
missing_taxids <- which(map_lgl(taxonomy, ~ all(is.na(.x)))) # need this positon info in a few lines below
not_missing_anymore_taxids <- taxize::classification(names(missing_taxids), db = 'ncbi')
for (i in seq_along(missing_taxids)) {
    taxon_info <- not_missing_anymore_taxids[[i]]
    new_taxon_name <- taxon_info |>
        filter(rank %in% taxonomic_ranks) |>
        pull(id) |>
        tail(1)
    taxonomy[[missing_taxids[i]]] <- taxon_info ## I need position here
    names(taxonomy)[[missing_taxids[i]]] <- new_taxon_name ## update taxid name to most recent
    mpa_data[mpa_data$taxid == names(missing_taxids)[i], 'taxid'] <- new_taxon_name

}
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
mpa_data <- left_join(mpa_data, new_taxonomy, by = 'taxid')
new_mpa_tree <- mpa_tree

tx <- paste0(taxonomic_ranks, '_taxid')
tx[which(tx == 'superkingdom_taxid')] <- 'kingdom_taxid'

mrcas <- flatten(map(tx, ~ split(mpa_data, factor(mpa_data[[.x]]))))
mrcas <-map(mrcas, ~ .x[['tip_label']])
mrcas <- map_int(mrcas, ~ getMRCA(mpa_tree, .x))
mrcas <- mrcas[!is.na(mrcas)]
mrcas_df <- data.frame(node = unname(mrcas), node_label = names(mrcas))
mrcas_df <- mrcas_df |>
    group_by(node) |>
    mutate(n_labels = length(unique(node_label))) |>
    mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    ungroup() |>
    distinct()

nodes <- data.frame(node = length(mpa_tree$tip.label) + 1:mpa_tree$Nnode) |>
    left_join(mrcas_df, by = 'node') |>
    mutate(
        node_label = ifelse(is.na(node_label), paste0('n', node), node_label),
        n_labels = ifelse(is.na(n_labels), 0, n_labels)
    )

mpa_tree$node.label <- nodes$node_label

nodes_with_taxid <- nodes |>
    filter(!grepl('^n', node_label)) |>
    separate_longer_delim(node_label, delim = '+')

nodes_taxonomy <- taxizedb::classification(
    x = unique(nodes_with_taxid$node_label), db = 'ncbi'
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
nodes_new_taxonomy <- nodes_new_taxonomy |>
    rename(node_label = taxid)

# sp_taxid_dups <- mpa_data$species_taxid[which(duplicated(mpa_data$species_taxid))]
#
# new_mpa_data <- mpa_data |>
#     # mutate(dup = ifelse(species_taxid %in% sp_taxid_dups, 'dup', NA)) |>
#     mutate(no_gca = as.double(sub('^GCA_', '', assembly_accession))) |>
#     slice_max(no_gca, n = 1, by = species_taxid) |>
#     # select(-no_gca, -dup) |>
#     select(-no_gca) |>
#     mutate(
#         old_tip_label = tip_label,
#         tip_label = as.character(species_taxid)
#     ) |>
#     select(
#         tip_label, genome_id, assembly_accession, old_tip_label,
#         kingdom_taxid, phylum_taxid, class_taxid, order_taxid, family_taxid,
#         genus_taxid, species_taxid
#     )
#
# tip_labels <- new_mpa_data$tip_label
# names(tip_labels) <- new_mpa_data$old_tip_label
# new_mpa_tree <- keep.tip(phy = mpa_tree, tip = names(tip_labels))
# new_tip_labels <- unname(tip_labels[new_mpa_tree$tip.label])
# new_mpa_tree$tip.label <- new_tip_labels
# new_mpa_data <- new_mpa_data[match(new_mpa_tree$tip.label, new_mpa_data$tip_label),]

# all(new_mpa_data$tip_label == new_mpa_tree$tip.label)




# new_mpa_data2 <- new_mpa_data |>
#     drop_na()

# taxRanks <- c(
#     'kingdom', 'phylum', 'class', 'order', 'family', 'genus'
# ) |>
#     {\(y) sub("$", "_taxid", y)}()
#
# mrcas <- map(taxRanks, ~ {
#     splitted_df <- split(new_mpa_data2, factor(new_mpa_data2[[.x]]))
#     output <- map_chr(splitted_df, \(y) as.character(getMRCA(t = new_mpa_tree, df = y)))
#     return(output)
# })
# names(mrcas) <- taxRanks

## values are the node numbers
## names are the taxids

# Export data -------------------------------------------------------------
new_mpa_tree_fname <- file.path('inst', 'extdata', 'mpav31.newick')
ape::write.tree(mpa_tree, new_mpa_tree_fname)

new_mpa_data_fname <- file.path('inst', 'extdata', 'mpav31_tips.tsv')
write.table(
    mpa_data, new_mpa_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

node_mpa_data_fname <- file.path('inst', 'extdata', 'mpav31_nodes.tsv')
write.table(
    nodes_new_taxonomy, node_mpa_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

