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




mpa_data[which(mpa_data$taxid != mpa_data$species_taxid),] |> View()



sp_dups <- mpa_data$species[which(duplicated(mpa_data$species))]

mpa_data |>
    filter(species %in% sp_dups) |>
    View()



mpa_data[which(duplicated(mpa_data$species_taxids)),] |> View()
mpa_data[which(duplicated(mpa_data$species_taxid)),] |> View()



output_fname <- file.path('inst', 'extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.tsv')

write.table(
    x = df, file = output_fname,
    sep = '\t', row.names = FALSE, quote = FALSE
)

