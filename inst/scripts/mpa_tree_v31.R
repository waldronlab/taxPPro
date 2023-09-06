
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

mpa_tree_fname <- system.file(
    'extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.nwk',
    package = 'taxPPro', mustWork = TRUE
)
mpa_tree <- ape::read.tree(mpa_tree_fname)
mpa_data <- data.frame(tip_label = mpa_tree$tip.label) |>
    separate(
        col = 'tip_label', into = column_names, sep = '\\|', remove = FALSE
    ) |>
    modify(~ sub('^\\w__', '', .x))

## I downloaded these data from https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ on 08292023
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

mpa_data <- left_join(mpa_data, assembly_data, by = 'genome_id') |>
    slice_max(order_by = version, n = 1, by = genome_id)

sp_taxid_dups <- mpa_data$species_taxid[which(duplicated(mpa_data$species_taxid))]

new_mpa_data <- mpa_data |>
    mutate(dup = ifelse(species_taxid %in% sp_taxid_dups, 'dup', NA)) |>
    mutate(no_gca = as.double(sub('^GCA_', '', assembly_accession))) |>
    slice_max(no_gca, n = 1, by = species_taxid) |>
    select(-no_gca, -dup) |>
    mutate(
        old_tip_label = tip_label,
        tip_label = species_taxid
    ) |>
    relocate(tip_label)

tip_labels <- new_mpa_data$tip_label
names(tip_labels) <- new_mpa_data$old_tip_label
new_mpa_tree <- keep.tip(phy = mpa_tree, tip = names(tip_labels))
new_tip_labels <- unname(tip_labels[new_mpa_tree$tip.label])
new_mpa_tree$tip.label <- new_tip_labels

new_mpa_tree_fname <- file.path('inst', 'extdata', 'mpav31.newick')
ape::write.tree(new_mpa_tree, new_mpa_tree_fname)

new_mpa_data <- new_mpa_data |>
    arrange(match(tip_label, new_tip_labels))

new_mpa_data_fname <- file.path('inst', 'extdata', 'mpvav31.tsv')
write.table(
    new_mpa_data, new_mpa_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)
