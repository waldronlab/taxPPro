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


url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'
sql <- '~/accessionTaxa.sql'
tree <- read.tree(url)


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


## ape write.tree won't allow commas and spaces when writin tree name
tree_data$tip_label <- gsub(' ', '_', tree_data$tip_label)
tree_data$tip_label <- gsub('[,;]', '-', tree_data$tip_label)
tree_data$tip_label <- gsub('"', '', tree_data$tip_label)
tree_data$tip_label <- gsub("'", '', tree_data$tip_label)
tree_data$tip_label <- gsub("\\[T\\]", '', tree_data$tip_label)
tree_data$tip_label <- gsub("\\(", '|', tree_data$tip_label)
tree_data$tip_label <- gsub("\\)", '|', tree_data$tip_label)
tree_data$tip_label <- gsub("\\[", '|', tree_data$tip_label)
tree_data$tip_label <- gsub("\\]", '|', tree_data$tip_label)
tree_data$tip_label <- gsub(":", '-', tree_data$tip_label)

tree$tip.label<- gsub(' ', '_', tree$tip.label)
tree$tip.label <- gsub('[,;]', '-', tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)
tree$tip.label <- gsub("'", '', tree$tip.label)
tree$tip.label <- gsub("\\[T\\]", '', tree$tip.label)
tree$tip.label <- gsub("\\(", '|', tree$tip.label)
tree$tip.label <- gsub("\\)", '|', tree$tip.label)
tree$tip.label <- gsub("\\]", '|', tree$tip.label)
tree$tip.label <- gsub("\\[", '|', tree$tip.label)
tree$tip.label <- gsub(":", '-', tree$tip.label)

all(tree$tip.label %in% tree_data$tip_label)

# Export data -------------------------------------------------------------
tree_fname <- file.path('inst', 'extdata', 'livingTree_2022.newick')
ape::write.tree(tree, tree_fname)

tree_data_fname <- file.path('inst', 'extdata', 'livingTree_tips_2022.tsv')
write.table(
    tree_data, tree_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

node_data_fname <- file.path('inst', 'extdata', 'livingTree_nodes_2022.tsv')
write.table(
    node_data, node_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)
