library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)
library(dplyr)

phys_names <- c('aerophilicity', 'growth temperature')
phys <- physiologies(phys_names, remove_false = TRUE, full_source = FALSE)
phys <- map(phys, ~ select(.x, -Accession_ID, -Genome_ID))
data_ready <- vector('list', length(phys))
for (i in seq_along(data_ready)) {
    message('Preparing ', names(phys)[i])
    names(data_ready)[i] <- names(phys)[i]
    data_ready[[i]] <- tryCatch(
        error = function(e) e,
        {
            prepareDatForPropagation(phys[[i]])
        }
    )
}
vct_lgl <- map_lgl(data_ready, is_error)
if (any(vct_lgl))  {
    message('Removing data with errors.')
    data_ready <- discard(data_ready, is_error)
}

data('tree_list')
tree <- as.Node(tree_list)

propagated <- vector('list', length(data_ready))
for (i in seq_along(propagated)) {
    message('Propagating ', names(data_ready)[i], ' - ', Sys.time())
    names(propagated)[i] <- names(data_ready)[i]
    propagated[[i]] <- tryCatch(
        error = function(e) e,
        {
            propagate(data_tree = tree, df = data_ready[[i]])
        }
    )
}
ncbi_taxonomy <- get_ncbi_taxonomy()
dfs <- map(propagated, ~ toDataFrame(.x, ncbi_tax = ncbi_taxonomy))

for (i in seq_along(dfs)) {
    if (names(dfs)[i] %in% names(phys)) {
        dfs[[i]]$Attribute_group <- unique(phys[[names(dfs[i])]]$Attribute_group)
        dfs[[i]]$Attribute_type <- unique(phys[[names(dfs[i])]]$Attribute_type)
    }
}

output <- vector('list', length(dfs))
for (i in seq_along(output)) {
    data <- dfs[[i]]
    names(output)[i] <- names(dfs)[i]
    attr_type <- unique(data$Attribute_type)
    if (attr_type == 'logical') {
        common_names <- grep('__', names(data), value = TRUE, invert = TRUE)
        unique_names <- grep('__', names(data), value = TRUE)
        attr_names <- unique(sub('__.*$', '', unique_names))
        data <- map(attr_names, ~ {
            names <- c(grep(.x, unique_names, value = TRUE), common_names)
            k <- data[,names]
            names(k) <- sub('.*__', '', names(k))
            k[['Attribute']] <- .x
            k[['Attribute_value']] <- TRUE
            cols <- c(
                'NCBI_ID', 'Taxon_name', 'Attribute', 'Attribute_value',
                'Evidence', 'Score'
            )
            k <- k[,cols]
            k <- k[which(!is.na(k$Evidence)),]
            k <- k[which(k$Evidence %in% c('asr', 'inh')),]
            k$Frequency <- scores2Freq(k$Score)
            k
        })
        names(data) <- attr_names
        data <- bind_rows(data)
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    } else if (attr_type == 'range') {
        cols2 <- c(
            'NCBI_ID', 'Taxon_name', 'Attribute',
            'Attribute_value_min', 'Attribute_value_max',
            'Evidence', 'Score'
        )
        attr_name <- unique(sub('__.*$', '', grep('__', names(data), value = TRUE)))
        colnames(data) <- sub('.*__', '', colnames(data))
        data[['Attribute']] <- sub('_', ' ', attr_name)
        data <- data[,cols2]
        data <- data[which(!is.na(data$Evidence)),]
        data <- data[which(data$Evidence %in% c('asr', 'inh')),]
        data$Frequency <- scores2Freq(data$Score)
        data <- bind_rows(data_ready[[names(dfs)[i]]], data)
        output[[i]] <- data
    }
}

full_dump <- reduce(output, bind_rows)
fname <- paste0("full_dump_bugphyzz_", Sys.Date(), ".csv.bz2")
unlink(fname)
con <- bzfile(fname, "w")
write.csv(full_dump, file = con, quote = TRUE)
close(con)

