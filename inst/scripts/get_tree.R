library(purrr)
library(dplyr)
library(taxPPro)
library(data.tree)
library(tidyr)
library(bugphyzz)

fname1 <- system.file(
    'extdata/proc_exclude_all_06292023.txt',
    package = 'taxPPro', mustWork = TRUE
)
exclude_all_ids <- read.table(fname1, header = FALSE)[[1]]
exclude_all_ranks <- taxizedb::taxid2rank(exclude_all_ids, db = 'ncbi')
exclude_all_ids_sp <- exclude_all_ids[which(exclude_all_ranks == 'species')]
exclude_all_sp_class <- taxizedb::classification(exclude_all_ids_sp, db = 'ncbi')
exclude_all_sp_tbl <- exclude_all_sp_class |>
    map(classif2Table) |>
    bind_rows() |>
    drop_na() |>
    filter(superkingdom %in% c('2', '2157'))

fname2 <- system.file(
    'extdata/proc_exclude_unclassified_uncultured_06292023.txt',
    package = 'taxPPro', mustWork = TRUE
)
exclude_uu_ids <- read.table(fname2, header = FALSE)[[1]]
exclude_uu_ranks <- taxizedb::taxid2rank(exclude_uu_ids, db = 'ncbi')
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

phys <- bugphyzz:::physiologies(remove_false = TRUE, full_source = FALSE)
bp_ids <- map(phys, ~ .x$NCBI_ID) |>
    unlist(use.names = FALSE, recursive = TRUE) |>
    unique() |>
    {\(y) y[y != 'unknown']}() |>
    {\(y) y[!is.na(y)]}()
bp_ranks <- taxizedb::taxid2rank(bp_ids, db = 'ncbi')
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
    distinct()

df <- data.frame(

    pathString = paste0(
        'ArcBac|',
        'd__', final_tbl$superkingdom,
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
tree_list <- as.list(tree)
usethis::use_data(tree_list, overwrite = TRUE)
