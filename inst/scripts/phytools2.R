## Script for bugphyzz exports

## Setup ####
library(logr)
library(bugphyzz)
library(taxPPro)
library(data.tree)
library(phytools)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ape)

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

## Import physiology data from bugphyzz ####
phys_names <- c(
    ## multistate-intersection
    'aerophilicity',

    ## multistate-union
    # 'antimicrobial resistance',

    ## binary
    'acetate producing'

    ## numeric/range
    # 'growth temperature'
)

msg <- paste0(
    'Importing ', length(phys_names), ' physiologies for propagation: ',
    paste0(phys_names, collapse = ', '), '.'
)
log_print(msg, blank_after = TRUE)
phys <- physiologies(phys_names)
v <- map_int(phys, nrow)
v <- sort(v)
phys <- phys[names(v)]

## Preparing data for propagation ####
msg <- ('Preparing data for propagation...')
log_print(msg, blank_after = TRUE)
tim <- system.time({
    phys_data_ready <- vector('list', length(phys))
    taxidWarnings <- vector('list', length(phys))
    for (i in seq_along(phys_data_ready)) {
        name <- names(phys)[i]
        msg <- paste0('Preparing ', name, '.')
        log_print(msg)
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
    phys_data_ready <- list_flatten(phys_data_ready)
})
log_print('Total time preparing data for propagation was: ')
log_print(tim, blank_after = TRUE)

taxidWarnings <- discard(taxidWarnings, is.null)
if (!is.null(taxidWarnings)) {
    msg <- 'Some NCBI IDs (taxids) need to be updated:'
    log_print(msg, blank_after = TRUE)
    log_print(taxidWarnings, blank_after = TRUE)
}

## Prepare tree data ####
msg <- paste0('Preparing tree data (NCBI and LTP).')
log_print(msg)
tim <- system.time({
    data('tree_list')
    ncbi_tree <- as.Node(tree_list)

    ltp <- ltp()
    tree <- reorder(ltp$tree, 'postorder')
    tip_data <- ltp$tip_data

    tx <- grep('_taxid$', colnames(tip_data), value = TRUE)
    nodes <- flatten(map(tx, ~ split(tip_data, factor(tip_data[[.x]]))))
    nodes <- map(nodes, ~ .x[['tip_label']])
    node_names <- map_int(nodes, ~ getMRCATaxPPro(tree, .x))
    node_names <- node_names[!is.na(node_names)]
    nodes_df <- data.frame(
        node = unname(node_names),
        node_label = names(node_names)
    ) |>
        group_by(node) |>
        mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
        ungroup() |>
        distinct()
})
log_print(tim, blank_after = TRUE)

## Propagation step ####
start_time <- Sys.time()
msg <- paste0('Performing propagation. It started at ', start_time, '.')
log_print(msg, blank_after = TRUE)

output <- vector('list', length(phys_data_ready))
for (i in seq_along(phys_data_ready)) {
    time1 <- Sys.time()

    ## Define variables
    current_phys <- names(phys_data_ready)[i]
    names(output)[i] <- current_phys
    dat <- phys_data_ready[[i]]
    Attribute_group_var <- unique(dat$Attribute_group)
    Attribute_group_var <- Attribute_group_var[!is.na(Attribute_group_var)]
    Attribute_type_var <- unique(dat$Attribute_type)
    Attribute_type_var <- Attribute_type_var[!is.na(Attribute_type_var)]

    dat_n_tax <- length(unique(dat$NCBI_ID))
    msg <- paste0(
        current_phys, ' has ', dat_n_tax, ' taxa.'
    )
    log_print(msg, blank_after = TRUE)

    ##  Mapping annotations to NCBI tree ####
    msg <- paste0(
        'Mapping source annotations to the NCBI tree for ', current_phys, '.'
    )
    log_print(msg)
    node_list <- split(
        x = dat, f = factor(dat$NCBI_ID)
    )

    tim <- system.time({
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
    })
    log_print(tim, blank_after = TRUE)

    ## Taxonomic pooling (round 1 of propagation) ####
    msg <- paste0(
        'Performing taxonomic pooling (round 1 of propagation) for ',
        current_phys, '.'
    )
    log_print(msg)
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
    })
    log_print(tim, blank_after = TRUE)

    ## Inheritance (round 1 of propagation) ####
    msg <- paste0(
        'Performing inhertiance1 (round 1 of propagation) for ',
        current_phys, '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(inh1, traversal = 'pre-order')
    })
    log_print(tim, blank_after = TRUE)

    new_dat <- ncbi_tree$Get(
        'attribute_tbl', filterFun = function(node) {
            grepl('^[gst]__', node$name)
        }
    ) |>
        discard(~ all(is.na(.x))) |>
        bind_rows() |>
        arrange(NCBI_ID, Attribute) |>
        filter(!NCBI_ID %in% dat$NCBI_ID) |>
        bind_rows(dat) # After this chunk is run, new_data also includes dat

    if (all(!new_dat$taxid %in% tip_data$taxid)) {
        msg <- paste0(
            'Not enough data for ASR. Skipping ASR for ', current_phys,
            '. Stopped after the first round of propagation.'
        )
        log_print(msg, blank_after = TRUE)
        output[[i]] <- new_dat
        next
    }

    ## Annotate pruned tree ####
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
    pruned_node_data <- data.frame(
        node = length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
    )

    tx <- grep('_taxid$', colnames(pruned_tip_data), value = TRUE)
    nodes <- tx |>
        map(~ split(pruned_tip_data, factor(pruned_tip_data[[.x]]))) |>
        flatten() |>
        map(~ .x[['tip_label']])
    node_names <- map_int(nodes, ~ getMRCATaxPPro(pruned_tree, .x))
    node_names <- node_names[!is.na(node_names)]
    nodes_df <- data.frame(
        node = unname(node_names),
        node_label = names(node_names)
    ) |>
        group_by(node) |>
        mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
        ungroup() |>
        distinct() |>
        arrange(node)
    pruned_node_data <- left_join(pruned_node_data, nodes_df, by = 'node') |>
        mutate(
            node_label = ifelse(
                is.na(node_label), paste0('n', as.character(node)), node_label
            )
        )
    pruned_tree$node.label <- pruned_node_data$node_label

    msg <- paste0(
        'Performing ASR for (round 2 of propagation) ', current_phys, '.'
    )
    log_print(msg)
    tim <- system.time({
        fit <- fitMk(
            tree = pruned_tree, x = annotated_tips, model = 'ER',
            pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
        )
        asr <- ancr(object = fit, tips = TRUE)
    })
    log_print(tim, blank_after = TRUE)

    res <- asr$ace
    node_rows <- length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
    rownames(res)[node_rows] <- pruned_tree$node.label

    ## Get annotations for nodes
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

    ## Perform taxonomic pooling and inheritance (propagation round 3)
    ## Mapping new internal nodes to the NCBI taxonomy tree ####
    msg <- paste0(
        'Mapping annotations for third round of propagation for ', current_phys,
        '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(function(node) {
            cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
            cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
            if (cond1 && cond2) {
                node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
            }
        })
    })
    log_print(tim, blank_after = TRUE)

    ## Taxonomic pooling (propagation round 3) ####
    msg <- paste0(
        'Performing taxonomic pooling (round 3 of propagation) for ',
        current_phys, '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(
            function(node_var) {
                taxPool(
                    node = node_var,
                    grp = Attribute_group_var,
                    typ = Attribute_type_var
                )
            },
            traversal = 'post-order'
        )

    })
    log_print(tim, blank_after = TRUE)

    ## Inheritance (propagation round 3) ####
    msg <- paste0(
        'Performing inheritance (round 3 of propagation) for ', current_phys
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(inh2, traversal = 'pre-order')
    })
    log_print(tim, blank_after = TRUE)

    ## Extracting files ####
    result <- ncbi_tree$Get(
        attribute = 'attribute_tbl', simplify = FALSE,
        filterFun = function(node) {
            node$name != 'ArcBac' && !is.null(node$attribute_tbl)
        }
    ) |>
        bind_rows() |>
        discard(~ all(is.na(.x)))
    min_thr <- 1 / length(unique(dat$Attribute))

    msg <- paste0(
        'Minimum threshold for positives in ', current_phys, ' was ',
        min_thr, '.'
    )
    log_print(msg, blank_after = TRUE)

    add_taxa_1 <- dat |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    add_taxa_2 <- new_taxa_for_ncbi_tree |>
        filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
        discard(~ all(is.na(.x)))
    final_result <- bind_rows(list(result, add_taxa_1, add_taxa_2)) |>
        filter(Score > min_thr)

    output[[i]] <- final_result

    msg <- paste0('Cleaning nodes for ', current_phys, '.')
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(cleanNode)
    })
    log_print(tim)

    time2 <- Sys.time()
    time3 <- round(difftime(time2, time1, units = 'min'))
    nrow_fr <- nrow(final_result)
    msg <- paste0(
        'Number of rows for ', current_phys, ' were' , nrow_fr , '.',
        'It took ', time3[[1]], ' mins.'
    )
    log_print(msg, blank_after = TRUE)
    log_print('', blank_after = TRUE)
}
end_time <- Sys.time()
elapsed_time <- round(difftime(end_time, start_time, units = 'min'))

msg <- paste0(
    'Total elapsed time for propagtion was ', elapsed_time[[1]], ' min.'
)
log_print(msg, blank_after = TRUE)

## Exporting annotations as a single tsv file ####
final_obj <- bind_rows(output)
msg <- paste0('Writing final output file.')
log_print(msg, blank_after = TRUE)
final_obj_fname <- paste0('bugphyzz_export_', Sys.Date(), '.tsv')
write.table(
    x = final_obj, file = final_obj_fname, sep = '\t', row.names = FALSE
)

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)
log_close()
