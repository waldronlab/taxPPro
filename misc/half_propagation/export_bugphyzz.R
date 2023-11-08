
suppressMessages({
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
    library(bugphyzzExports)
})

logfile <- "log_file"
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

phys_names <- c(
    ## multistate-intersection
     'aerophilicity',
     'gram stain',
     'biosafety level',
     'COGEM pathogenicity rating',
     'shape',
     'spore shape',
     'arrangement',

    ## multistate-union
     # 'habitat'
#    'isolation site',
#    'disease association',
     # 'antimicrobial resistance',

    # 'growth medium',

    ## binary
     'plant pathogenicity',
     'acetate producing',
     'sphingolipid producing',
     'lactate producing',
     'butyrate producing',
     'hydrogen gas producing',
     'pathogenicity human',
     'motility',
     'biofilm forming',
     'extreme environment',
     'animal pathogen',
     'antimicrobial sensitivity',

    ## numeric/range
     'growth temperature',
     'optimal ph',
     'width',
     'length',
     'genome size',
     'coding genes',
     'mutation rate per site per generation',
     'mutation rate per site per year'
)

msg <- paste0(
    'Importing ', length(phys_names), ' physiologies for propagation: ',
    paste0(phys_names, collapse = ', '), '.'
)
log_print(msg, blank_after = TRUE)
bugphyzz_data <- physiologies(phys_names)
v_order <- sort(map_int(bugphyzz_data, nrow))
bugphyz_data <- bugphyzz_data[names(v_order)]

msg <- paste0(
    'Searching for attributes of type range. They will be converted to type ',
    'multistate-intersection based on thresholds.'
)
log_print(msg, blank_after = TRUE)

phys <- vector('list', length(bugphyzz_data))
for (i in seq_along(phys)) {
    attribute_type <- bugphyzz_data[[i]]$Attribute_type |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    dat_name <- names(bugphyzz_data)[i]
    names(phys)[i] <- dat_name
    if (attribute_type == 'range' && dat_name %in% names(THRESHOLDS())) {
        msg <- paste0(
            dat_name, " is of type range and we have a threshold for it.",
            ' Converting ', dat_name, ' to multistate-intersection.'
        )
        log_print(msg)
        res <- rangeToLogicalThr(bugphyzz_data[[i]], THRESHOLDS()[[dat_name]])
        res$Attribute_type <- 'multistate-intersection'
        phys[[i]] <- res

    } else if (attribute_type == 'range' && !dat_name %in% names(THRESHOLDS())) {
        msg <- paste0(
            dat_name, " is of type range, but we don't have a threshold for it.",
            " Skipping ", dat_name, '.'
        )
        log_print(msg)
        next

    } else {
        phys[[i]] <- bugphyzz_data[[i]]

    }
}

for (i in seq_along(phys)) {
   physName <-  names(phys)[i]
   if (is.null(phys[[i]])) {
       msg <- paste0(physName, ' will be discarded now. Not in thresholds.')
       log_print(msg)
   }
}
log_print("", blank_after = TRUE)
phys <- discard(phys, is.null)

msg <- ('Preparing data for propagation...')
log_print('', blank_after = TRUE)
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
    ## list_flatten ensures that data from the attributes of type
    ## multistate-union are separated into data.frames
})


for (i in seq_along(phys_data_ready)) {
   physName <-  names(phys_data_ready)[i]
   if (is.null(phys_data_ready[[i]])) {
       msg <- paste0(phyName, ' will be discarded now. No taxids.')
       log_print(msg)
   }
}
log_print("", blank_after = TRUE)
phys_data_ready <- discard(phys_data_ready, is.null)

phys_data_ready <- map(phys_data_ready, ~ {
    attribute_type <- .x |>
        pull(Attribute_type) |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    if (attribute_type %in% c('binary', 'multistate-union')) {
        ## This step is used to include FALSE values,
        ## which are necessary for the ASR step below
        ## In the case of multistate-intersection, FALSE values are inferred because they're mutally exclusive (need to elaborate more here).
        return(completeBinaryData(.x))
    }
    return(.x)
})

log_print('', blank_after = TRUE)
log_print('Total time preparing data for propagation was: ')
log_print(tim, blank_after = TRUE)

taxidWarnings <- discard(taxidWarnings, is.null)
if (!is.null(taxidWarnings)) {
    msg <- 'Some NCBI IDs (taxids) might need to be updated:'
    log_print(msg, blank_after = TRUE)
    log_print(taxidWarnings, blank_after = TRUE)
}

msg <- paste0('Preparing tree data (NCBI and LTP).')
log_print(msg)
tim <- system.time({
    data('tree_list')
    ncbi_tree <- as.Node(tree_list)
    ltp <- ltp()
    tree <- reorder(ltp$tree, 'postorder')
    tip_data <- ltp$tip_data
    node_data <- ltp$node_data

})
log_print(tim, blank_after = TRUE)

start_time <- Sys.time()
msg <- paste0('Performing propagation. It started at ', start_time, '.')
log_print(msg, blank_after = TRUE)

output <- vector('list', length(phys_data_ready))
for (i in seq_along(phys_data_ready)) {
    time1 <- Sys.time()
    dat <- phys_data_ready[[i]]

    attribute_group <- dat$Attribute_group |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    attribute_type <- dat$Attribute_type |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    attribute_nms <- dat$Attribute |>
        {\(y) y[!is.na(y)]}() |>
        unique()

    names(output)[i] <- names(phys_data_ready)[i]

    if (attribute_type == 'multistate-union') {
        attrNMS <- unique(sub('--(TRUE|FALSE)$', '', attribute_nms))
        attrGroupMsg <- paste0(attribute_group, '-', attrNMS)
    } else {
        attrGroupMsg <- attribute_group
    }

    dat_n_tax <- length(unique(dat$NCBI_ID))
    msg <- paste0(
        attrGroupMsg, ' has ', format(dat_n_tax, big.mark = ','), ' taxa.'
    )
    log_print(msg, blank_after = TRUE)

    msg <- paste0(
        'Mapping source annotations to the NCBI tree for ', attrGroupMsg, '.'
    )
    log_print(msg)

    node_list <- split(dat, factor(dat$NCBI_ID))
    tim <- system.time({
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
    })
    log_print(tim, blank_after = TRUE)

    msg <- paste0(
        'Performing taxonomic pooling for ',
        attrGroupMsg, '.'
    )
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(
           function(node) {
                taxPool(
                    node = node,
                    grp = attribute_group,
                    typ = attribute_type)
            },
            traversal = 'post-order'
        )
    })
    log_print(tim, blank_after = TRUE)

    msg <- paste0(
        'Performing inheritance (1) for ',
        attrGroupMsg, '.'
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
        bind_rows(dat)

    new_taxids <- new_dat |>
        pull(taxid) |>
        {\(y) y[!is.na(y)]}()

    per <- mean(tip_data$taxid %in% new_taxids) * 100
    if (per < 1) {
        msg <- paste0(
            'Not enough data for ASR. Skipping ASR and inhetiance (2) for ', attrGroupMsg,
            '. Stopped after the first round of propagation.'
        )
        log_print(msg, blank_after = TRUE)

        output[[i]] <- new_dat ## Here, I include data after taxpool and inhertiance 1, ASR and inheritance 2 are skipped

        msg <- paste0('Cleaning nodes for ', attrGroupMsg, '.')
        log_print(msg)
        tim <- system.time({
            ncbi_tree$Do(cleanNode)
        })
        log_print(tim, blank_after = TRUE)

        time2 <- Sys.time()
        time3 <- round(difftime(time2, time1, units = 'min'))
        nrow_fr <- nrow(new_dat)
        msg <- paste0(
            'Number of rows for ', attrGroupMsg, ' were ' ,
            format(nrow_fr, big.mark = ','), '.',
            ' It took ', time3[[1]], ' mins.'
        )
        log_print(msg, blank_after = TRUE)
        log_print('', blank_after = TRUE)

        next
    }

    final_result <- new_dat

    output[[i]] <- final_result

    msg <- paste0('Cleaning nodes for ', attrGroupMsg, '.')
    log_print(msg)
    tim <- system.time({
        ncbi_tree$Do(cleanNode)
    })
    log_print(tim, blank_after = TRUE)

    time2 <- Sys.time()
    time3 <- round(difftime(time2, time1, units = 'min'))
    nrow_fr <- nrow(final_result)
    msg <- paste0(
        'Number of rows for ', attrGroupMsg, ' were ' ,
        format(nrow_fr, big.mark = ','), '.',
        ' It took ', time3[[1]], ' mins.'
    )
    log_print(msg, blank_after = TRUE)
    log_print('', blank_after = TRUE)
}
end_time <- Sys.time()
elapsed_time <- round(difftime(end_time, start_time, units = 'min'))

msg <- paste0(
    'Propagation ended at ', elapsed_time,
    '. Total elapsed time for propagtion for ', length(phys_data_ready),
    ' physiologies was ', elapsed_time[[1]], ' min.'
)
log_print(msg, blank_after = TRUE)

## Exporting annotations as a single tsv file ####
final_obj <- bind_rows(output)
final_obj_size <- lobstr::obj_size(final_obj)

msg <- paste0(
    'Size of final object is ',
    gdata::humanReadable(final_obj_size, standard = 'SI')
)
log_print(msg, blank_after = TRUE)

msg <- paste0('Writing final output file.')
log_print(msg, blank_after = TRUE)
final_obj_fname <- paste0('bugphyzz_export_', Sys.Date(), '.tsv')
write.table(
    x = final_obj, file = final_obj_fname, sep = '\t', row.names = FALSE
)

fsize <- gdata::humanReadable(file.size(final_obj_fname), standard = "SI")
msg <- paste0(
    'The size of the tsv file is ', fsize, '. Output file name is ',
    final_obj_fname, '.'
)
log_print(msg, blank_after = TRUE)

si <- sessioninfo::session_info()
log_print(si, blank_after = TRUE)
log_close()
