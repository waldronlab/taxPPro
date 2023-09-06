
## This code is for converting the tips of the metaphlan tree v3.1 from
## genome ids (GCA) to taxids

library(ape)
library(tidyr)
library(purrr)
library(dplyr)

column_names <- c(
    'genome_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
    'species'
)

getMRCA <- function(t, df) {
    res <- phytools::findMRCA(t, tips = df[['tip_label']])
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
    select(assembly_accession, taxid, species_taxid, organism_name)
assembly_data_historical <- readr::read_tsv(
    assembly_data_historical_fname, show_col_types = FALSE, skip = 1
) |>
    {\(y) set_names(y, sub('#', '', colnames(y)))}() |>
    select(assembly_accession, taxid, species_taxid, organism_name)
assembly_data <- bind_rows(assembly_data_current, assembly_data_historical) |>
    mutate(genome_id = sub('\\.\\d+', '', assembly_accession)) |>
    mutate(version = as.integer(sub('^.*\\.(\\d+)$', '\\1', assembly_accession)))

mpa_data <- data.frame(tip_label = mpa_tree$tip.label) |>
    separate(
        col = 'tip_label', into = column_names, sep = '\\|', remove = FALSE
    ) |>
    modify(~ sub('^\\w__', '', .x)) |>
    left_join(assembly_data, by = 'genome_id') |>
    slice_max(order_by = version, n = 1, by = genome_id) |>
    select(
        tip_label, genome_id, assembly_accession, original_taxid = taxid,
        original_species_taxid = species_taxid
    ) |>
    mutate(original_taxid = as.character(original_taxid))

taxonomy <- taxizedb::classification(x = unique(mpa_data$original_taxid), db = 'ncbi')
missing <- which(map_lgl(taxonomy, ~ all(is.na(.x))))
no_missing_anymore <- taxize::classification(names(missing), db = 'ncbi')
for (i in seq_along(missing)) {
    taxonomy[[missing[i]]] <- no_missing_anymore[[i]]
}

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

mpa_data <- left_join(mpa_data, new_taxonomy, by = 'original_taxid')

sp_taxid_dups <- mpa_data$species_taxid[which(duplicated(mpa_data$species_taxid))]

new_mpa_data <- mpa_data |>
    mutate(dup = ifelse(species_taxid %in% sp_taxid_dups, 'dup', NA)) |>
    mutate(no_gca = as.double(sub('^GCA_', '', assembly_accession))) |>
    slice_max(no_gca, n = 1, by = species_taxid) |>
    select(-no_gca, -dup) |>
    mutate(
        old_tip_label = tip_label,
        tip_label = as.character(species_taxid)
    ) |>
    select(
        tip_label, genome_id, assembly_accession, old_tip_label,
        kingdom_taxid, phylum_taxid, class_taxid, order_taxid, family_taxid,
        genus_taxid, species_taxid
    )

tip_labels <- new_mpa_data$tip_label
names(tip_labels) <- new_mpa_data$old_tip_label
new_mpa_tree <- keep.tip(phy = mpa_tree, tip = names(tip_labels))
new_tip_labels <- unname(tip_labels[new_mpa_tree$tip.label])
new_mpa_tree$tip.label <- new_tip_labels
new_mpa_data <- new_mpa_data |>
    arrange(match(tip_label, new_tip_labels))

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
ape::write.tree(new_mpa_tree, new_mpa_tree_fname)

new_mpa_data_fname <- file.path('inst', 'extdata', 'mpav31.tsv')
write.table(
    new_mpa_data, new_mpa_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)
