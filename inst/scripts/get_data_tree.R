## Script used to generate tree_list (ncbi_tree)

library(purrr)
library(dplyr)
library(taxPPro)
library(data.tree)
library(tidyr)
library(bugphyzz)
library(treeio)

# taxizedb::db_download_ncbi(overwrite = TRUE)

fname1 <- system.file(
    'extdata/proc_species_07122023.txt',
    package = 'taxPPro', mustWork = TRUE
)
exclude_all_ids <- read.table(fname1, header = FALSE)[[1]]
exclude_all_ranks <- taxizedb::taxid2rank(exclude_all_ids, db = 'ncbi')

pos <- which(!is.na(exclude_all_ranks))
exclude_all_ids <- exclude_all_ids[pos]
exclude_all_ranks <- exclude_all_ranks[pos]
names(exclude_all_ranks) <- exclude_all_ids
exclude_all_names <- taxizedb::taxid2name(exclude_all_ids, db = 'ncbi')
names(exclude_all_names) <- exclude_all_ids

exclude_all_ids_sp <- exclude_all_ids[which(exclude_all_ranks == 'species')]
exclude_all_sp_class <- taxizedb::classification(exclude_all_ids_sp, db = 'ncbi')
exclude_all_sp_tbl <- exclude_all_sp_class |>
    map(classif2Table) |>
    bind_rows() |>
    drop_na() |>
    filter(superkingdom %in% c('2', '2157'))

fname2 <- system.file(
    'extdata/proc_strains_07122023.txt',
    package = 'taxPPro', mustWork = TRUE
)
exclude_uu_ids <- read.table(fname2, header = FALSE)[[1]]
exclude_uu_ranks <- taxizedb::taxid2rank(exclude_uu_ids, db = 'ncbi')


pos2 <- which(!is.na(exclude_uu_ranks))
exclude_uu_ids <- exclude_uu_ids[pos2]
exclude_uu_ranks <- exclude_uu_ranks[pos2]
names(exclude_uu_ranks) <- exclude_uu_ids
exclude_uu_names <- taxizedb::taxid2name(exclude_uu_ids, db = 'ncbi')
names(exclude_uu_names) <- exclude_uu_ids

exclude_uu_ids_st <- exclude_uu_ids[which(exclude_uu_ranks == 'strain')]
exclude_uu_st_class <- taxizedb::classification(exclude_uu_ids_st, db = 'ncbi')
exclude_uu_st_tbl <- exclude_uu_st_class |>
    map(classif2Table) |>
    bind_rows() |>
    drop_na() |>
    filter(superkingdom %in% c('2', '2157')) |>
    select(species, strain)
merged_sp_st <- left_join(
    exclude_all_sp_tbl, exclude_uu_st_tbl, by = 'species',
) |>
    distinct() |>
    mutate(strain = ifelse(is.na(strain), '', strain))
phys <- bugphyzz:::physiologies(full_source = FALSE)
bp_ids <- map(phys, ~ .x$NCBI_ID) |>
    unlist(use.names = FALSE, recursive = TRUE) |>
    unique() |>
    {\(y) y[y != 'unknown']}() |>
    {\(y) y[!is.na(y)]}()
bp_ranks <- taxizedb::taxid2rank(bp_ids, db = 'ncbi')


pos3 <- which(!is.na(bp_ranks))
bp_ids <- bp_ids[pos3]
bp_ranks <- bp_ranks[pos3]
bp_names <- taxizedb::taxid2name(bp_ids, db = 'ncbi')
names(bp_names) <- bp_ids

names(bp_ranks) <- bp_ids
bp_ids_sp <- names(bp_ranks)[which(bp_ranks == 'species')]
extra_ids_sp <- bp_ids_sp[!bp_ids_sp %in% merged_sp_st$species]
extra_class_sp <- taxizedb::classification(extra_ids_sp, db = 'ncbi')
extra_tbl_sp <- map(extra_class_sp, classif2Table) |>
    bind_rows() |>
    drop_na() |>
    relocate(
        superkingdom, phylum, class, order, family, genus, species
    ) |>
    filter(superkingdom %in% c('2', '2157')) |>
    distinct()

final_tbl <- bind_rows(merged_sp_st, extra_tbl_sp) |>
    distinct() |>
    mutate(strain = ifelse(is.na(strain), '', strain))
df <- data.frame(
    pathString = paste0(
        'ArcBac|',
        'k__', final_tbl$superkingdom,
        '|p__', final_tbl$phylum,
        '|c__', final_tbl$class,
        '|o__', final_tbl$order,
        '|f__', final_tbl$family,
        '|g__', final_tbl$genus,
        '|s__', final_tbl$species,
        '|t__', final_tbl$strain
    ))
df$pathString <- sub('\\|t__$', '', df$pathString)
tree <- as.Node(df[, 'pathString', drop = FALSE], pathDelimiter = '|')


x <- final_tbl |>
    rename(kingdom = superkingdom) |>
    map(unique)
ranks <- vector('list', length(x))
for (i in seq_along(x)) {
    ranks[[i]] <- rep(names(x)[i], length(x[[i]]))
    names(ranks[[i]]) <- x[[i]]

}
ranks <- flatten_chr(ranks)
ranks <- ranks[names(ranks) != '']
names <- taxizedb::taxid2name(names(ranks), db = 'ncbi')
names(names) <- names(ranks)

tree$Do(function(node) {
    if (!node$isRoot) {
        taxid <- sub('^\\w__', '', node$name)
        node$Taxon_name <- names[taxid]
        node$taxid <- taxid
        node$Rank <- ranks[taxid]
    }
})

tree_list <- as.list(tree)

tree_test <- as.Node(tree_list)


tbl_for_tree <- final_tbl |>
    rename(domain = superkingdom) |>
    mutate(
        root = 'ArcBac',
        domain = paste0('k__', domain),
        phylum = paste0('p__', phylum),
        class = paste0('c__', class),
        order = paste0('o__', order),
        family = paste0('f__', family),
        genus = paste0('g__', genus),
        species = paste0('s__', species)

    ) |>
    relocate(root) |>
    select(-strain) |>
    distinct()

tree_sp <- TreeSummarizedExperiment::toTree(tbl_for_tree)
tree_sp$node.label <- sub('^.*:', '', tree_sp$node.label)

usethis::use_data(tree_list, overwrite = TRUE)
## usethis::use_data(tree_sp, overwrite = TRUE)
## This whole procedure takes about 10 minutes (local machine).
