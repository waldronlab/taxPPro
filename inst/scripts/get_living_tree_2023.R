
# Setup -------------------------------------------------------------------
library(ape)
library(taxonomizr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(phytools)


# Import tree -------------------------------------------------------------

## The file accessionTaxa.sql was downloaded with taxonomizr
## The file is kind of large, so better to downloaded just one time
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

# Adjust branches with zero length ----------------------------------------
pos_zero <- which((tree$edge[,2] %in% 1:Ntip(tree)) & (tree$edge.length == 0))
tree$edge.length[pos_zero] <- tree$edge.length[pos_zero] + 1e-05

# Accession to taxids -----------------------------------------------------
accession_rgx <- '_-.*--_'
accessions <- tree$tip.label |>
    str_extract(accession_rgx) |>
    str_remove_all('[_-]')
taxnames <- tree$tip.label |>
    str_extract(paste0("^.*", accession_rgx)) |>
    str_remove(accession_rgx) |>
    str_remove('-$') |>
    str_replace_all('_', ' ') |>
    str_squish()
taxids <- accessionToTaxa(accessions, sql, version = 'base') # base doesn't include decimal version

# Complete missing taxa ---------------------------------------------------
missing_taxa <- taxnames[which(is.na(taxids))]
## Below, the two only cases that need to be updated manually
missing_taxa[which(missing_taxa == 'Micromonospora okii')] <- 'Micromonospora sp. TP-A0468'
missing_taxa[which(missing_taxa == 'Sala cibi')] <- 'Salella cibi'
not_missing_anymore_taxa <- taxizedb::name2taxid(missing_taxa, db = 'ncbi')
taxids[which(is.na(taxids))] <- not_missing_anymore_taxa

# Get full taxonomy -------------------------------------------------------
taxonomy <- taxizedb::classification(unique(taxids), db = 'ncbi')
missing_taxonomy_positions <- which(map_lgl(taxonomy, ~ all(is.na(.x))))
not_missing_anymore_taxonomy <- taxize::classification(
    names(missing_taxonomy_positions), db =  'ncbi'
)
chr_vct <- map_chr(not_missing_anymore_taxonomy, ~ tail(.x$id, 1))
for (i in seq_along(chr_vct)) {
    pos <- which(taxids == names(chr_vct[i]))
    message('Replacing ', unique(taxids[pos]), ' with ', chr_vct[i])
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
    taxid = taxids
) |>
    mutate(Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi')) |>
    mutate(Rank = taxizedb::taxid2rank(taxid, db = 'ncbi')) |>
    mutate(NCBI_ID = taxPPro::addRankPrefix(taxid, Rank)) |>
    left_join(new_taxonomy, by = 'taxid')
rownames(tip_data) <- tip_data$tip_label

# node_data ---------------------------------------------------------------
getMRCA <- function(tree, tips) {
    res <- phytools::findMRCA(tree = tree, tips = tips)
    if (is.null(res))
        res <- NA
    res
}

ranks_ids <- paste0(taxonomic_ranks, '_taxid')
ranks_ids[which(ranks_ids == 'superkingdom_taxid')] <- 'kingdom_taxid'
ranks_ids <- ranks_ids[!ranks_ids %in% c('species_taxid', 'strain_taxid')]

mrcas_df <- tip_data |>
    select(-species_taxid, -strain_taxid) |> # only do this for genus and above
    {\(y) purrr::map(ranks_ids, ~ split(y, factor(y[[.x]])))}() |>
    flatten() |>
    purrr::map(~ .x[['tip_label']]) |>
    map_int(~ getMRCA(tree, .x)) |>
    {\(y) y[!is.na(y)]}() |>
    {\(y) data.frame(node = unname(y), node_label = names(y))}() |>
    group_by(node) |>
    mutate(n_labels = length(unique(node_label))) |>
    mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    ungroup() |>
    distinct() |>
    arrange(node)
tree$node.label <- mrcas_df$node_label[match(Ntip(tree) + 1:Nnode(tree), mrcas_df$node)]

nodes_taxonomy <- taxizedb::classification(
    x = unique(unlist(strsplit(mrcas_df$node_label, '\\+'))),
    db = 'ncbi'
) |> purrr::map(~ {
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
        order_taxid, family_taxid, genus_taxid
    ) |>
    discard(~ all(is.na(.x)))

node_data <- mrcas_df |>
    separate_longer_delim(node_label, delim = '+') |>
    rename(taxid = node_label) |>
    group_by(node) |>
    mutate(node_label = paste0(taxid, collapse = '+')) |>
    ungroup() |>
    mutate(Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi')) |>
    mutate(Rank = taxizedb::taxid2rank(taxid, db = 'ncbi')) |>
    mutate(NCBI_ID = taxPPro::addRankPrefix(taxid, Rank)) |>
    left_join(nodes_taxonomy, by = 'taxid') |>
    as.data.frame()

all(node_data$node_label %in% tree$node.label)

# Add genus information ---------------------------------------------------
genus_list <- tip_data |>
    filter(!is.na(genus_taxid)) |>
    {\(y) split(y, factor(y$genus_taxid))}()
names(genus_list) <- paste0('g__', names(genus_list))
tree_extended <- tree

## The next step took about 25 minutes
system.time({
    for (i in seq_along(genus_list)) {
        myTips <- genus_list[[i]]$tip_label
        message('Adding ', names(genus_list)[i], ' - ', i, '/', length(genus_list), '. This one contains ', length(myTips), ' tip(s).')
        if (length(myTips) == 1) {
            ## It the tip is alone, just add another tip with the genus.
            ## This creates a new internal node with the same edge length, and two child tips with branch length of 0
            tip_number <- which(tree_extended$tip.label == myTips)
            tree_extended <- bind.tip(
                tree = tree_extended, tip.label = names(genus_list)[i],
                edge.length = 0, ## This really doesn't matter. edge.length is always 0 when adding a tip to a tip. That's why I need to add lengths of 1e-06 a few lines below.
                where = tip_number
            )
            tn1 <- which(tree_extended$tip.label == myTips) ## tip and node numbers change when binding a tip
            tn2 <- which(tree_extended$tip.label == names(genus_list)[i])
            tree_extended$edge.length[which(tree_extended$edge[,2] == tn1)] <- 1e-06
            tree_extended$edge.length[which(tree_extended$edge[,2] == tn2)] <- 1e-06
        } else {
            ## Here, I need to find the MRCA again because the node numbers change when
            ## a new node is added to the tree
            myMRCA <- findMRCA(tree = tree_extended, tips = myTips, type = 'node')
            tree_extended <- bind.tip(
                tree = tree_extended, edge.length = 1e-06, where = myMRCA,
                tip.label = names(genus_list)[i]
            )
        }
    }
})

## This is just a small check that the tips are being mapped to the right
## internal node
all_labels <- c(tree_extended$tip.label, tree_extended$node.label)
myMRCA2 <- findMRCA(tree = tree_extended, tips = genus_list[[3709]]$tip_label, type = 'node')
all_labels[myMRCA2] == sub('g__', '', names(genus_list)[3709])

extra_tip_data <- data.frame(
    tip_label = grep('g__', tree_extended$tip.label, value = TRUE)
)
# extra_tip_data$accession <- NA # just create an explicit column with NAs for combining later with tip_data
extra_tip_data$taxid <- sub('^g__', '', extra_tip_data$tip_label)
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


## Update node numbers in node data
node_data_extended <- node_data
node_data_extended$node <- match(node_data$node_label, tree_extended$node.label) + Ntip(tree_extended)

## Some checks (all should have name) no 'NA's (character string not an NA value)
sum(tree_extended$node.label[node_data_extended$node - Ntip(tree_extended)] == 'NA')
tree_extended$node.label[23504 - Ntip(tree_extended)] == '2157'

## note that node_data_extended, only contains data about genera with
## two or more tips, i.e. I could get a MRCA. The rest are only tips.
## It's probably better to filter out genera from node_data when
## importing with the ltp function and use only the tips for managing
## genus information.

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
    node_data_extended, node_data_fname, sep = '\t', quote = TRUE,
    row.names = FALSE
)
