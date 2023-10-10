
## Load packages ####
library(bugphyzz)
library(taxPPro)
library(data.tree)
library(phytools)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ape)

phys_names <- c(
    ## multistate-intersection
    'aerophilicity',

    ## multistate-union
    'antimicrobial resistance',

    ## binary
    'acetate producing'
)

phys <- physiologies(phys_names)

phys_data_ready <- vector('list', length(phys))
taxidWarnings <- vector('list', length(phys))
for (i in seq_along(phys_data_ready)) {
    name <- names(phys)[i]
    message('Preparing ', name, '.')
    names(phys_data_ready)[i] <- name
    names(taxidWarnings)[i] <- name
    wngs <- list()
    suppressWarnings({
        withCallingHandlers({
            dat <- getDataReady(filterData(phys[[i]]))
            if (length(dat) > 0)
                phys_data_ready[[i]] <- dat
        },
        warning = function(w) {
            if (grepl('taxizedb', w$message)) {
                msg <- sub('.*unrank.*: (\\d+.*)$', '\\1', w$message)
                wngs <<- c(wngs, list(msg))
            }
        })
    })
    if (length(wngs) > 0)
        taxidWarnings[[i]] <- wngs
}
taxidWarnings <- discard(taxidWarnings, is.null)
if (!is.null(taxidWarnings)) {
    message('Some NCBI IDs (taxids) need to be updated:')
    print(taxidWarnings)
}
phys_data_ready <- list_flatten(phys_data_ready)

## Load NCBI taxonomy tree ####
data('tree_list')
ncbi_tree <- as.Node(tree_list)

## Load the living tree project (LTP) tree and tip data ####
ltp <- ltp()
tree <- reorder(ltp$tree, 'postorder')
tip_data <- ltp$tip_data

tx <- grep('_taxid$', colnames(tip_data), value = TRUE)
nodes <- flatten(map(tx, ~ split(tip_data, factor(tip_data[[.x]]))))
nodes <- map(nodes, ~ .x[['tip_label']])
## integers are nodes in tree and names is the taxid for that node
node_names <- map_int(nodes, ~ getMRCATaxPPro(tree, .x))
node_names <- node_names[!is.na(node_names)]
nodes_df <- data.frame(
    node = unname(node_names),
    node_label = names(node_names)
    ) |>
    group_by(node) |>
    # mutate(n_labels = length(unique(node_label))) |>
    mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    ungroup() |>
    distinct()

start_time <- Sys.time()
output <- vector('list', length(phys_data_ready))
for (i in seq_along(phys_data_ready)) {

    ## Define variables
    current_phys <- names(phys_data_ready)[i]
    names(output)[i] <- current_phys
    dat <- phys_data_ready[[i]]
    Attribute_group_var <- unique(dat$Attribute_group)
    Attribute_group_var <- Attribute_group_var[!is.na(Attribute_group_var)]
    Attribute_type_var <- unique(dat$Attribute_type)
    Attribute_type_var <- Attribute_type_var[!is.na(Attribute_type_var)]

    ##  Mappint to NCBI tree
    message('Mapping source annotations to the NCBI tree for ', current_phys, '.')
    node_list <- split(
        x = dat, f = factor(dat$NCBI_ID)
    )

    tim <- system.time({
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
    })
    print(tim)

    ## Taxonomic pooling and inheritance (round 1)
    message('Performing round 1 of propagation for ', current_phys, '.')
    tim <- system.time({
        ncbi_tree$Do(
           function(node) {
                taxPool(
                    node = node,
                    grp = Attribute_group_var,
                    typ = Attribute_type_var)
            },
            traversal = 'post-order'
        )
        ncbi_tree$Do(inh1, traversal = 'pre-order')
    })
    print(tim)

    new_dat <- ncbi_tree$Get(
        'attribute_tbl', filterFun = function(node) {
            grepl('^[gst]__', node$name)
        }
    ) |>
        discard(~ all(is.na(.x))) |>
        bind_rows() |>
        arrange(NCBI_ID, Attribute) |>
        filter(!NCBI_ID %in% dat$NCBI_ID) |>
        # mutate(taxid = sub('^\\w__', '', NCBI_ID)) |>
        bind_rows(dat)

    if (all(!new_dat$taxid %in% tip_data$taxid)) {
        message('Not enough data for ASR. Skipping ', current_phys, '.')
        next
    }

    ## Perform ASR with phytools
    tip_data_annotated <- left_join(
        tip_data,
        select(new_dat, taxid, Attribute, Score),
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

    pruned_tree <- ape::keep.tip(tree, tip = rownames(annotated_tips))
    pruned_tree <- reorder(pruned_tree, 'postorder')

    pruned_tip_data <- tip_data |>
        filter(tip_label %in% pruned_tree$tip.label)

    ## Annotated the pruned tree

    pruned_node_data <- data.frame(
        node = length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
    )

    tx <- grep('_taxid$', colnames(pruned_tip_data), value = TRUE)
    nodes <- tx |>
        map(~ split(pruned_tip_data, factor(pruned_tip_data[[.x]]))) |>
        flatten() |>
        map(~ .x[['tip_label']])

    ## integers are nodes in tree and names is the taxid for that node
    node_names <- map_int(nodes, ~ getMRCATaxPPro(pruned_tree, .x))
    node_names <- node_names[!is.na(node_names)]
    nodes_df <- data.frame(
        node = unname(node_names),
        node_label = names(node_names)
    ) |>
        group_by(node) |>
        # mutate(n_labels = length(unique(node_label))) |>
        mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
        ungroup() |>
        distinct() |>
        arrange(node)
    node_data <- left_join(pruned_node_data, nodes_df, by = 'node') |>
        mutate(
            node_label = ifelse(
                is.na(node_label), paste0('n', as.character(node)), node_label
            )
        )

    pruned_tree$node.label <- node_data$node_label

    message('Performing ASR for ', current_phys, '.')
    tim <- system.time({
        fit <- fitMk(
            tree = pruned_tree, x = annotated_tips, model = 'ER',
            pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
        )
        asr <- ancr(object = fit, tips = TRUE)
    })

    res <- asr$ace
    node_rows <- length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
    rownames(res)[node_rows] <- pruned_tree$node.label

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
    new_taxa_for_ncbi_tree <- new_taxa_from_nodes |>
        relocate(NCBI_ID, Rank, Attribute, Score, Evidence)
    new_taxa_for_ncbi_tree_list <- split(
        new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
    )

    ## Perform taxonomic pooling and inheritance (round 2)
    message('Performing second taxPool and Inh pooling for', current_phys)
    ncbi_tree$Do(function(node) {
        cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
        cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
        if (cond1 && cond2) {
            node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
        }
    })
    ncbi_tree$Do(
        function(node) {
            taxPool(
                node = node,
                grp = Attribute_group_var,
                typ = Attribute_type_var
            )
        },
        traversal = 'post-order'
    )
    ncbi_tree$Do(inh2, traversal = 'pre-order')

    ## Extracting files
    result <- ncbi_tree$Get(
        attribute = 'attribute_tbl', simplify = FALSE,
        filterFun = function(node) {
            node$name != 'ArcBac' && !is.null(node$attribute_tbl)
        }
    ) |>
        bind_rows() |>
        discard(~ all(is.na(.x)))
    min_thr <- 1 / length(unique(dat$Attribute))
    add_taxa_1 <- dat |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    add_taxa_2 <- new_taxa_for_ncbi_tree |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    final_result <- bind_rows(list(result, add_taxa_1, add_taxa_2)) |>
        filter(Score > min_thr)

    output[[i]] <- final_result

    message('Cleaning nodes for ', current_phys, '.')
    tim <- system.time({
        ncbi_tree$Do(cleanNode)
    })
    print(tim)


}
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)






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
# no_annotated_tips_chr <- tip_data_annotated |>
#     filter(!tip_label %in% rownames(annotated_tips)) |>
#     pull(tip_label) |>
#     unique()
# no_annotated_tips <- matrix(
#     data = rep(rep(1/ncol(annotated_tips), ncol(annotated_tips)), length(no_annotated_tips_chr)),
#     nrow = length((no_annotated_tips_chr)),
#     byrow = TRUE,
#     dimnames = list(
#         rownames = no_annotated_tips_chr, colnames = colnames(annotated_tips)
#     )
# )

# input_matrix <- rbind(annotated_tips, no_annotated_tips)
# input_matrix <- input_matrix[tree$tip.label,]


## add scode for ASR

tree <- ape::keep.tip(tree, tip = rownames(annotated_tips))
tree <- reorder(tree, 'postorder')
input_matrix <- annotated_tips

## ------------------------------------------------------------------------------------------------------------------------
system.time({
    fit <- fitMk(
        tree = tree, x = input_matrix, model = 'ER', pi = 'fitzjohn',
        lik.func = 'pruning', logscale = TRUE
    )
})
asr <- ancr(object = fit, tips = TRUE)
res <- asr$ace
node_rows <- length(tree$tip.label) + 1:tree$Nnode
rownames(res)[node_rows] <- tree$node.label

## ------------------------------------------------------------------------------------------------------------------------
# new_taxa_from_tips <- res[rownames(no_annotated_tips),] |>
#     as.data.frame() |>
#     tibble::rownames_to_column(var = 'tip_label') |>
#     left_join(tip_data, by = 'tip_label') |>
#     mutate(
#         Rank = taxizedb::taxid2rank(taxid, db = 'ncbi')
#     ) |>
#     filter(Rank %in% c('genus', 'species', 'strain')) |>
#     mutate(
#         NCBI_ID = case_when(
#             Rank == 'genus' ~ paste0('g__', taxid),
#             Rank == 'species' ~ paste0('s__', taxid),
#             Rank == 'strain' ~ paste0('t__', taxid)
#         )
#     ) |>
#     select(-ends_with('_taxid'), -tip_label, -taxid, -accession) |>
#     relocate(NCBI_ID, Taxon_name, Rank) |>
#     pivot_longer(
#         names_to = 'Attribute', values_to = 'Score', cols = 4:last_col()
#     ) |>
#     mutate(
#         Evidence = 'asr',
#         Attribute_source = NA,
#         Confidence_in_curation = NA,
#         taxid = sub('\\w__', '', NCBI_ID),
#         Attribute_type = Attribute_type_var,
#         Attribute_group = Attribute_group_var,
#         Frequency = case_when(
#             Score == 1 ~ 'always',
#             Score > 0.9 ~ 'usually',
#             Score >= 0.5 ~ 'sometimes',
#             Score > 0 & Score < 0.5 ~ 'rarely',
#             Score == 0 ~ 'never'
#         )
#     )





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
    filterFun = function(node) node$name != 'ArcBac' && !is.null(node$attribute_tbl)
    ) |>
    bind_rows() |>
    discard(~ all(is.na(.x)))
min_thr <- 1 / length(unique(phys_data_ready$Attribute))
add_taxa_1 <- phys_data_ready |>
    filter(!NCBI_ID %in% unique(output$NCBI_ID)) |>
    discard(~ all(is.na(.x)))
add_taxa_2 <- new_taxa_for_ncbi_tree |>
    filter(!NCBI_ID %in% unique(output$NCBI_ID)) |>
    discard(~ all(is.na(.x)))
final_output <- bind_rows(list(output, add_taxa_1, add_taxa_2)) |>
    filter(Score > min_thr)
