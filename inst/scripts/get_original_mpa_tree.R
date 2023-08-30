library(ape)
library(tidyr)
library(purrr)

col_names <- c(
    'genome_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
    'species'
)

t_fname <- system.file('extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.nwk', package = 'taxPPro', mustWork = TRUE)
t <- ape::read.tree(t_fname)
df <- data.frame(tip_label = t$tip.label)
df <- df |>
    separate(col = 'tip_label', into = col_names, sep = '\\|', remove = FALSE) |>
    modify(~ sub('^\\w__', '', .x))

assembly_fname <- '~/assembly_summary_genbank.txt' # I downloaded this from https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ on 08292023
assembly <- readr::read_tsv(assembly_fname, show_col_types = FALSE)

acc <- assembly$assembly_accession
acc <- sub('\\.\\d+', '', acc)


missing <- df$genome_id[which(!df$genome_id %in% acc)]




output_fname <- file.path('inst', 'extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.tsv')

write.table(
    x = df, file = output_fname,
    sep = '\t', row.names = FALSE, quote = FALSE
)

