library(ape)
library(tidyr)
fname <- system.file('extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.nwk', package = 'taxPPro', mustWork = TRUE)
t <- ape::read.tree(fname)
df <- data.frame(tip_label = t$tip.label)
df <- df |>
    separate(col = 'tip_label', into = col_names, sep = '\\|', remove = FALSE) |>
    modify(~ sub('^\\w__', '', .x))
output_fname <- file.path('inst', 'extdata', 'mpa_v31_CHOCOPhlAn_201901_species_tree.tsv')
write.table(
    x = df, file = output_fname,
    sep = '\t', row.names = FALSE, quote = FALSE
)

