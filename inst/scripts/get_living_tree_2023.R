
library(ape)
library(taxonomizr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(phytools)

sql <- '~/accessionTaxa.sql'
treeUrl <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_08_2023.ntree'
treeFilePath <- file.path('inst', 'extdata', 'LTP_all_08_2023.ntree')
if (file.exists(treeFilePath)) {
    message('File already cached. Imported from local machine.')
    tree <- read.tree(treeFilePath)
} else {
    message("File doesn't exist. Downloading from the LTP website.")
    tree <- read.tree(treeUrl)
    tree$tip.label <- gsub('"', '', tree$tip.label)
    tree$tip.label <- gsub("'", '', tree$tip.label)
    tree$tip.label <- gsub(' ', '_', tree$tip.label)
    tree$tip.label <- gsub('[,;:\\(\\)]', '-', tree$tip.label)
    message("Caching file.")
    write.tree(tree, treeFilePath)
}


# Accession to taxids -----------------------------------------------------
accession_rgx <- '_-.*--_'
accessions <- str_extract(tree$tip.label, accession_rgx) |>
    str_remove_all('[_-]')
taxnames <- str_extract(tree$tip.label, paste0("^.*", accession_rgx)) |>
    str_remove(accession_rgx) |>
    str_remove('-$') |>
    str_replace_all('_', ' ') |>
    str_squish()

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
new_taxonomy <- purrr::map(taxonomy, ~ {
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

# tip_data data.frame -----------------------------------------------------
tip_data <- data.frame(
    tip_label = tree$tip.label,
    accession = accessions,
    taxid = taxids,
    taxname = taxnames,
    rank = taxizedb::taxid2rank(taxids, db = 'ncbi')

) |>
    left_join(new_taxonomy, by = 'taxid')
rownames(tip_data) <- tip_data$tip_label

# node_data ---------------------------------------------------------------
getMRCA <- function(tree, tips) {
    res <- phytools::findMRCA(tree = tree, tips = tips)
    if (is.null(res))
        res <- NA
    res
}

tx <- paste0(taxonomic_ranks, '_taxid')
tx[which(tx == 'superkingdom_taxid')] <- 'kingdom_taxid'
mrcas <- flatten(purrr::map(tx, ~ split(tip_data, factor(tip_data[[.x]]))))
mrcas <- purrr::map(mrcas, ~ .x[['tip_label']])
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
nodes_new_taxonomy <- purrr::map(nodes_taxonomy, ~ {
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
node_data$rank <- taxizedb::taxid2rank(node_data$taxid, db = 'ncbi')


node_data$NCBI_ID <- taxPPro::addRankPrefix(node_data$taxid, node_data$rank)
tip_data$NCBI_ID <- taxPPro::addRankPrefix(tip_data$taxid, tip_data$rank)

## ape write.tree won't allow commas and spaces when writin tree name
# tree_data$tip_label <- gsub(' ', '_', tree_data$tip_label)
# tree_data$tip_label <- gsub('[,;]', '-', tree_data$tip_label)
# tree_data$tip_label <- gsub('"', '', tree_data$tip_label)
# tree_data$tip_label <- gsub("'", '', tree_data$tip_label)
# tree_data$tip_label <- gsub("\\[T\\]", '', tree_data$tip_label)
#
# tree$tip.label<- gsub(' ', '_', tree$tip.label)
# tree$tip.label <- gsub('[,;]', '-', tree$tip.label)
# tree$tip.label <- gsub('"', '', tree$tip.label)
# tree$tip.label <- gsub("'", '', tree$tip.label)
# tree$tip.label <- gsub("\\[T\\]", '', tree$tip.label)


all(tree$tip.label %in% tip_data$tip_label)


# tip_data$taxid
tip_data$Taxon_name <- taxizedb::taxid2name(tip_data$taxid, db = 'ncbi')
tip_data$Rank <- taxizedb::taxid2rank(tip_data$taxid, db = 'ncbi')
tip_data$rank <- NULL


# node_data$taxid
node_data$Taxon_name <- taxizedb::taxid2name(node_data$taxid, db = 'ncbi')
node_data$Rank <- taxizedb::taxid2rank(node_data$taxid, db = 'ncbi')
node_data$rank <- NULL


# Adjust tips with zero length --------------------------------------------
## All tips are in the second column in the matrix tree$edge

## I was thinking about adding the minimun branch length to the
## ones with zero, but this might change the result
# tips_node_positions <- which(tree$edge[,2] %in% 1:Ntip(tree))
# tip_branch_lenghts <- tree$edge.length[tips_node_positions]
# min_length_for_tip <- min(tip_branch_lenghts[tip_branch_lenghts > 0])

pos_zero <- which((tree$edge[,2] %in% 1:Ntip(tree)) & (tree$edge.length == 0))
tree$edge.length[pos_zero] <- tree$edge.length[pos_zero] + 1e-06

# Add genus information ---------------------------------------------------

# node_data_g <- node_data[which(node_data$Rank == 'genus'),]$node_label

rl <- tip_data |>
    filter(!is.na(genus_taxid)) |>
    {\(y) split(y, factor(y$genus_taxid))}() |>
    purrr::keep(~ nrow(.x) > 1)
names(rl) <- paste0('g__', names(rl))
tree_extended <- tree
system.time({
    for (i in seq_along(rl)) {
        ## Here, I need to find the MRCA again because the node numbers change when
        ## a new node is added to the tree
        message('Adding ', names(rl)[i], ' - ', i, '/', length(rl))
        myMRCA <- findMRCA(tree = tree_extended, tips = rl[[i]]$tip_label, type = 'node')
        tree_extended <- bind.tip(
            tree = tree_extended, edge.length = 1e-06, where = myMRCA,
            tip.label = names(rl)[i]
        )
    }
})


## This is just a small check that the tips are being mapped to the right
## internal node
all_labels <- c(tree_extended$tip.label, tree_extended$node.label)
myMRCA2 <- findMRCA(tree = tree_extended, tips = rl[[i]]$tip_label, type = 'node')
all_labels[myMRCA2] == sub('g__', '', names(rl)[i])


# node_data_g <- node_data[which(node_data$Rank == 'genus'),]$node_label
# names(node_data_g) <- node_data[which(node_data$Rank == 'genus'),]$taxid
# names(node_data_g) <- paste0('g__', names(node_data_g))
#
# tree_extended <- tree
# system.time({
#     for (i in seq_along(node_data_g)) {
#         message('Adding ', names(node_data_g)[i], ' - ', i, '/', length(node_data_g))
#         tree_extended <- bind.tip(
#             tree = tree_extended, edge.length = 1e-06, where = node_data_g[i],
#             tip.label = names(node_data_g)[i]
#         )
#     }
# })

extra_tip_data <- data.frame(tip_label = grep('g__', tree_extended$tip.label, value = TRUE))
extra_tip_data$accession <- NA
extra_tip_data$taxid <- sub('g__', '', extra_tip_data$tip_label)
extra_tip_data$NCBI_ID <- extra_tip_data$tip_label
extra_tip_data$Rank <- 'genus'
extra_tip_taxonomy <- taxizedb::classification(
    x = unique(extra_tip_data$taxid), db = 'ncbi'
)
extra_tip_new_taxonomy <- purrr::map(extra_tip_taxonomy, ~ {
    x <- .x
    x
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
        order_taxid, family_taxid, genus_taxid
    ) |>
    discard(~ all(is.na(.x)))

extra_tip_data <- left_join(extra_tip_data, extra_tip_new_taxonomy, by = 'taxid')
extra_tip_data$Taxon_name <- taxizedb::taxid2name(extra_tip_data$taxid, db = 'ncbi')
rownames(extra_tip_data) <- extra_tip_data$tip_label
tip_data_extended <- bind_rows(tip_data, extra_tip_data)
tip_data_extended <- tip_data_extended[tree_extended$tip.label,]

# Export data -------------------------------------------------------------
tree_fname <- file.path('inst', 'extdata', 'LTP_all_08_2023.newick')
# ape::write.tree(tree, tree_fname)
ape::write.tree(tree_extended, tree_fname)

tip_data_fname <- file.path('inst', 'extdata', 'LTP_all_08_2023.tip_data')
write.table(
    # tip_data, tip_data_fname, sep = '\t', quote = TRUE,
    tip_data_extended, tip_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)

node_data_fname <- file.path('inst', 'extdata', 'LTP_all_08_2023.node_data')
write.table(
    node_data, node_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)
