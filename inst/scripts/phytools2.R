## ----include = FALSE-----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, message=FALSE------------------------------------------------------------------------------------------------
library(bugphyzz)
library(taxPPro)
library(data.tree)
library(phytools)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)


## ----import physiology, message=FALSE------------------------------------------------------------------------------------
phys_name <- 'acetate producing'
phys <- physiologies(phys_name)


## ----warning=FALSE-------------------------------------------------------------------------------------------------------
phys_data_ready <- phys[[1]] |>
    filterData() |>
    getDataReady()
phys_data_list <- split(phys_data_ready, factor(phys_data_ready$NCBI_ID))


## ------------------------------------------------------------------------------------------------------------------------
data('tree_list')
ncbi_tree <- as.Node(tree_list)


## ------------------------------------------------------------------------------------------------------------------------
ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data


## ------------------------------------------------------------------------------------------------------------------------
ncbi_tree$Do(function(node) {
    if (node$name %in% names(phys_data_list)) {
        node$attribute_tbl <- phys_data_list[[node$name]]
    } else {
        node$attribute_tbl <- NULL
    }
})


## ------------------------------------------------------------------------------------------------------------------------
Attribute_group_var <- unique(phys_data_ready$Attribute_group)
Attribute_group_var <- Attribute_group_var[!is.na(Attribute_group_var)]
Attribute_type_var <- unique(phys_data_ready$Attribute_type)
Attribute_type_var <- Attribute_type_var[!is.na(Attribute_type_var)]
ncbi_tree$Do(
    function(node) taxPool(node = node, grp = Attribute_group_var, typ = Attribute_type_var),
    traversal = 'post-order'
)
ncbi_tree$Do(inh1, traversal = 'pre-order')


## ------------------------------------------------------------------------------------------------------------------------
new_phys_data <- ncbi_tree$Get(
    'attribute_tbl', filterFun = function(node) grepl('^[gst]__', node$name)
) |>
    discard(~ all(is.na(.x))) |>
    bind_rows() |>
    arrange(NCBI_ID, Attribute) |>
    filter(!NCBI_ID %in% phys_data_ready$NCBI_ID) |>
    mutate(taxid = sub('^\\w__', '', NCBI_ID)) |>
    bind_rows(phys_data_ready)

tip_data_annotated <- left_join(
    tip_data,
    select(new_phys_data, taxid, Attribute, Score),
    by = 'taxid'
    )

annotated_tips <- tip_data_annotated |>
    select(tip_label, Attribute, Score) |>
    filter(!is.na(Attribute)) |>
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score', values_fill = 0
    ) |>
    tibble::column_to_rownames(var = 'tip_label') |>
    as.matrix()
no_annotated_tips_chr <- tip_data_annotated |>
    filter(!tip_label %in% rownames(annotated_tips)) |>
    pull(tip_label) |>
    unique()
no_annotated_tips <- matrix(
    data = rep(rep(1/ncol(annotated_tips), ncol(annotated_tips)), length(no_annotated_tips_chr)),
    nrow = length((no_annotated_tips_chr)),
    byrow = TRUE,
    dimnames = list(
        rownames = no_annotated_tips_chr, colnames = colnames(annotated_tips)
    )
)

input_matrix <- rbind(annotated_tips, no_annotated_tips)
input_matrix <- input_matrix[tree$tip.label,]


## ------------------------------------------------------------------------------------------------------------------------
fit <- fitMk(
    tree = tree, x = input_matrix, model = 'ER', pi = 'estimated',
    lik.func = 'pruning', logscale = TRUE
)
asr <- ancr(object = fit, tips = TRUE)
res <- asr$ace
node_rows <- length(tree$tip.label) + 1:tree$Nnode
rownames(res)[node_rows] <- tree$node.label


## ------------------------------------------------------------------------------------------------------------------------
new_taxa_from_tips <- res[rownames(no_annotated_tips),] |>
    as.data.frame() |>
    tibble::rownames_to_column(var = 'tip_label') |>
    left_join(tip_data, by = 'tip_label') |>
    mutate(
        Rank = taxizedb::taxid2rank(taxid, db = 'ncbi')
    ) |>
    filter(Rank %in% c('genus', 'species', 'strain')) |>
    mutate(
        NCBI_ID = case_when(
            Rank == 'genus' ~ paste0('g__', taxid),
            Rank == 'species' ~ paste0('s__', taxid),
            Rank == 'strain' ~ paste0('t__', taxid)
        )
    ) |>
    select(-ends_with('_taxid'), -tip_label, -taxid, -accession) |>
    relocate(NCBI_ID, Taxon_name, Rank) |>
    pivot_longer(
        names_to = 'Attribute', values_to = 'Score', cols = 4:last_col()
    ) |>
    mutate(
        Evidence = 'asr',
        Attribute_source = NA,
        Confidence_in_curation = NA,
        taxid = sub('\\w__', '', NCBI_ID),
        Attribute_type = Attribute_type_var,
        Attribute_group = Attribute_group_var,
        Frequency = case_when(
            Score == 1 ~ 'always',
            Score > 0.9 ~ 'usually',
            Score >= 0.5 ~ 'sometimes',
            Score > 0 & Score < 0.5 ~ 'rarely',
            Score == 0 ~ 'never'
        )
    )

nodes_annotated <- res[which(grepl('^\\d+(\\+\\d+)*', rownames(res))),]
new_taxa_from_nodes <- nodes_annotated |>
    as.data.frame() |>
    tibble::rownames_to_column(var = 'NCBI_ID') |>
    filter(grepl('^\\d+(\\+\\d+)*', NCBI_ID)) |>
    mutate(NCBI_ID = strsplit(NCBI_ID, '\\+')) |>
    tidyr::unnest(NCBI_ID) |>
    mutate(Rank = taxizedb::taxid2rank(NCBI_ID)) |>
    mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
    mutate(
        NCBI_ID = case_when(
            Rank == 'kingdom' ~ paste0('k__', NCBI_ID),
            Rank == 'phylum' ~ paste0('p__', NCBI_ID),
            Rank == 'class' ~ paste0('c__', NCBI_ID),
            Rank == 'order' ~ paste0('o__', NCBI_ID),
            Rank == 'family' ~ paste0('f__', NCBI_ID),
            Rank == 'genus' ~ paste0('g__', NCBI_ID),
            Rank == 'species' ~ paste0('s__', NCBI_ID),
            Rank == 'strain' ~ paste0('t__', NCBI_ID)
        )
    ) |>
    filter(
        Rank %in% c(
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
            'species', 'strain'
        )
    ) |>
    mutate(Evidence = 'asr') |>
    relocate(NCBI_ID, Rank, Evidence) |>
    pivot_longer(
        cols = 4:last_col(), names_to = 'Attribute', values_to = 'Score'
    ) |>
    mutate(
        Attribute_source = NA,
        Confidence_in_curation = NA,
        Attribute_group = Attribute_group_var,
        Attribute_type = Attribute_type_var,
        taxid = sub('\\w__', '', NCBI_ID),
        Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi'),
        Frequency = case_when(
            Score == 1 ~ 'always',
            Score > 0.9 ~ 'usually',
            Score >= 0.5 ~ 'sometimes',
            Score > 0 & Score < 0.5 ~ 'rarely',
            Score == 0 ~ 'never'
        )
    )

new_taxa_for_ncbi_tree <- new_taxa_from_tips |>
    relocate(NCBI_ID, Rank, Attribute, Score, Evidence) |>
    bind_rows(new_taxa_from_nodes)
new_taxa_for_ncbi_tree_list <- split(
    new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
)



## ------------------------------------------------------------------------------------------------------------------------
ncbi_tree$Do(function(node) {
    cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
    cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
    if (cond1 && cond2) {
        node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
    }
})


## ------------------------------------------------------------------------------------------------------------------------
ncbi_tree$Do(
    function(node) taxPool(node = node, grp = Attribute_group_var, typ = Attribute_type_var ),
    traversal = 'post-order'
)
ncbi_tree$Do(inh2, traversal = 'pre-order')


## ------------------------------------------------------------------------------------------------------------------------
output <- ncbi_tree$Get(
    attribute = 'attribute_tbl', simplify = FALSE,
    filterFun = function(node) node$name != 'ArcBac'
    ) |>
    bind_rows()
min_thr <- 1 / length(unique(phys_data_ready$Attribute))
add_taxa_1 <- phys_data_ready |>
    filter(!NCBI_ID %in% unique(output$NCBI_ID))
add_taxa_2 <- new_taxa_for_ncbi_tree |>
    filter(!NCBI_ID %in% unique(output$NCBI_ID))
final_output <- bind_rows(list(output, add_taxa_1, add_taxa_2)) |>
    filter(Score > min_thr)

