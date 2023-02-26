
## A function for removing prefix from a character vector x
removePrefix <- function(x) {
    regex <- '^[dkpcofgst]__'
    sub(regex, '', x)
}

#' Add strains
#'
#' @param node A node.
#'
#' @return An action on the node.
#' @export
#'
addStrains <- function(node) {
    if (node$isLeaf) {
        current_children <- removePrefix(names(node$children))
        current_taxon <- removePrefix(node$name)
        df <- tryCatch(
            error = function(e) NULL, {
                taxizedb::children(current_taxon)[[1]]
            }
        )
        if (is.null(df)) {
            return(NULL)
        }
        df <- df[df$rank %in% 'strain', ]
        if (!nrow(df)) {
            return(NULL)
        }
        children <- df$id
        new_children <- children[!children %in% current_children]
        new_children <- paste0('t__', new_children)
        if (isFALSE(!length(new_children))) {
            for (i in seq_along(new_children)) {
                message(new_children[i])
                node$AddChild(new_children[i])
            }
        }
    }

}


#' Fill NAs
#'
#' \code{fillNAs} changes NAs wit NA\[number\].
#'
#' @param x An atomic vector.
#'
#' @return An atomic vector
#' @export
#'
fillNAs <- function(x) {
    counter <- 1
    for (i in seq_along(x)) {
        if (is.na(x[i])) {
            x[i] <- paste0('NA', counter)
            counter <- counter + 1
        }
    }
    return(x)
}

#' Conver taxname to taxid
#'
#'\code{taxname2taxid} Converts a taxon name from the ncbi to its corresponding
#'taxid
#'
#' @param tax_tbl An table with taxnames.
#'
#' @return A data.frame.
#' @export
#'
taxname2taxid <- function(tax_tbl) {
    if ('kingdom' %in% colnames(tax_tbl)) {
        pos <- which(colnames(tax_tbl) == 'kingdom')
        colnames(tax_tbl)[pos] <- 'superkingdom'
    }
    ranks <- c(
        'superkingdom', 'phylum', 'class', 'order', 'family', 'genus'
    )
    for (i in ranks) {
        df <- tax_tbl[tax_tbl[['Rank']] == i, ]
        vct <- df$NCBI_ID
        names(vct) <- df[[i]]
        tax_tbl[[i]] <- vapply(tax_tbl[[i]], function(.x) vct[.x], character(1))
        # tax_tbl[[i]] <- purrr::map_chr(tax_tbl[[i]], ~ {vct[.x]})
    }
    return(tax_tbl)
}


#' Get NCBI data (tree or table)
#'
#' \code{getNCBI} gets data from the NCBI in tree or table format. Tree format
#' is in the node format of the data.tree package. For now, it only inlcudes
#' taxa from the bacteria superkingdom/domain/kingdom.
#'
#' @param format A character string. Options: 'table' or 'tree'.
#'
#' @return A data.frame or a tree in data.tree format (Node, R6) .
#' @export
#'
getNCBI <- function(format = 'table') {
    ncbi_taxids <- get_ncbi_taxids(keyword = 'b', with_taxids = TRUE)
    new_ncbi_taxids <- taxname2taxid(tax_tbl = ncbi_taxids)
    cond1 <- new_ncbi_taxids$Rank == 'species'
    # cond2 <- new_ncbi_taxids$superkingdom == '2'
    # bacteria <- new_ncbi_taxids[cond1 & cond2,]
    bacteria <- new_ncbi_taxids[cond1,]
    no_cols <- c('species', 'strain', 'Taxon_name', 'Parent_NCBI_ID', 'Rank')
    bacteria <- bacteria[,!colnames(bacteria) %in% no_cols]
    # for (i in seq_along(bacteria)) {
    #   bacteria[[i]] <- fillNAs(bacteria[[i]])
    # }
    bacteria <- tidyr::drop_na(bacteria)
    pathString <- paste(
        'ArcBac', # root
        '|||d__', bacteria$superkingdom,
        '|||p__', bacteria$phylum,
        '|||c__', bacteria$class,
        '|||o__', bacteria$order,
        '|||f__', bacteria$family,
        '|||g__', bacteria$genus,
        '|||s__', bacteria$NCBI_ID,
        sep = ''
    )
    bacteria$pathString <- pathString
    bacteria$NCBI_ID <- paste0('s__', bacteria$NCBI_ID)
    if (format == 'table') {
        return(bacteria)
    } else if (format == 'tree') {
        tree <- data.tree::as.Node(bacteria, pathDelimiter = '|||')
        return(tree)
    }
}

#' Add attributes to a data.tree object
#'
#' \code{addAttributes} Adds attributes to a data.tree object
#'
#' @param data_tree  A data.tree object.
#' @param df A data.frame. Output of \code{prepareData2}
#'
#' @return A clone of the original data.tree object.
#' @export
#'
addAttributes <- function(data_tree, df) {
    attr_type <- unique(df[['Attribute_type']])
    if (attr_type == 'logical') output <- addAttributesLogical(data_tree, df)
    if (attr_type == 'numeric') output <- addAttributesNumeric(data_tree, df)
    if (attr_type == 'range') output <- addAttributesRange(data_tree, df)
    return(output)
}

#' Add attributes to a data.tree (numeric)
#'
#' \code{addAttributesLogical} adds attributes to a data.tree object.
#'
#' @param data_tree A data.tree object.
#' @param df A data.frame from bugphyzz. Output of \code{prepareData2}.
#'
#' @return A clone of the original data.tree.
#' @export
#'
addAttributesLogical <- function(data_tree, df) {
    datatree <- data.tree::Clone(data_tree)
    attr_grp <- unique(df[['Attribute_group']])
    df[['Attribute']] <- gsub(' ', '_', df[['Attribute']])
    attributes <- unique(df[['Attribute']])
    for (i in seq_along(attributes)) {
        attr <- attributes[i]
        attr_df <- df[df[['Attribute']] == attr, ]
        colnames(attr_df) <- paste0(attr_grp, ':', attr, '__', colnames(attr_df))
        ncbi_col <- grep('__NCBI_ID$', colnames(attr_df), value = TRUE)
        evi_col <- grep('__Evidence$', colnames(attr_df), value = TRUE)
        source_col <- grep('__Attribute_source$', colnames(attr_df), value = TRUE)
        score_col <- grep('__Score$', colnames(attr_df), value = TRUE)
        ncbi_list <- split(attr_df, factor(attr_df[[ncbi_col]]))
        datatree$Do(function(node) node[[score_col]] <- ncbi_list[[node$name]][[score_col]])
        datatree$Do(function(node) node[[evi_col]] <- ncbi_list[[node$name]][[evi_col]])
        datatree$Do(function(node) node[[source_col]] <- ncbi_list[[node$name]][[source_col]])
    }
    return(datatree)
}

#' Add attributes to a data.tree object (numeric)
#'
#' \code{addAttributesNumeric} Add numeric attributes to a data.tree object.
#'
#' @param data_tree  A data.tree object.
#' @param df A data.frame. Output of \code{prepareData2}
#'
#' @return A clone of the data.tree object.
#' @export
#'
addAttributesNumeric <- function(data_tree, df) {
    datatree <- data.tree::Clone(data_tree)
    attr_grp <- unique(df[['Attribute_group']])
    colnames(df) <- paste0(attr_grp, '__', colnames(df))
    colnames(df) <- gsub(' ', '_', colnames(df))
    attr_val_col <- grep('__Attribute_value$', colnames(df), value = TRUE)
    ncbi_col <- grep('__NCBI_ID$', colnames(df), value = TRUE)
    evi_col <- grep('__Evidence$', colnames(df), value = TRUE)
    source_col <- grep('__Attribute_source$', colnames(df), value = TRUE)
    score_col <- grep('__Score$', colnames(df), value = TRUE)
    ncbi_list <- split(df, factor(df[[ncbi_col]]))
    datatree$Do(function(node) node[[attr_val_col]] <- ncbi_list[[node$name]][[attr_val_col]])
    datatree$Do(function(node) node[[score_col]] <- ncbi_list[[node$name]][[score_col]])
    datatree$Do(function(node) node[[evi_col]] <- ncbi_list[[node$name]][[evi_col]])
    datatree$Do(function(node) node[[source_col]] <- ncbi_list[[node$name]][[source_col]])
    return(datatree)
}

#' Add attributes to a data.tree object (range)
#'
#' \code{addAttributesRange} adds range attributes to a data.tree object.
#'
#' @param data_tree A data.tree object.
#' @param df A data.frame. Output of \code{prepareData2}
#'
#' @return A clone of the input data.tree object with extra attributes.
#' @export
#'
addAttributesRange <- function(data_tree, df) {
    datatree <- data.tree::Clone(data_tree)
    attr_grp <- unique(df[['Attribute_group']])
    colnames(df) <- paste0(attr_grp, '__', colnames(df))
    colnames(df) <- gsub(' ', '_', colnames(df))
    attr_val_min_col <- grep('__Attribute_value_min$', colnames(df), value = TRUE)
    attr_val_max_col <- grep('__Attribute_value_max$', colnames(df), value = TRUE)
    ncbi_col <- grep('__NCBI_ID$', colnames(df), value = TRUE)
    evi_col <- grep('__Evidence$', colnames(df), value = TRUE)
    source_col <- grep('__Attribute_source$', colnames(df), value = TRUE)
    score_col <- grep('__Score$', colnames(df), value = TRUE)
    ncbi_list <- split(df, factor(df[[ncbi_col]]))
    datatree$Do(function(node) node[[attr_val_min_col]] <- ncbi_list[[node$name]][[attr_val_min_col]])
    datatree$Do(function(node) node[[attr_val_max_col]] <- ncbi_list[[node$name]][[attr_val_max_col]])
    datatree$Do(function(node) node[[score_col]] <- ncbi_list[[node$name]][[score_col]])
    datatree$Do(function(node) node[[evi_col]] <- ncbi_list[[node$name]][[evi_col]])
    datatree$Do(function(node) node[[source_col]] <- ncbi_list[[node$name]][[source_col]])
    return(datatree)
}

#' ASR / Upstream
#'
#' \code{asrUpstream}
#'
#' @param node
#'
#' @return Node/R6. A data.tree object
#' @export
#'
asrUpstream <- function(node) {
    attrs <- grep('__Score$', node$attributesAll, value = TRUE)
    ## If node is leaf
    cond1 <- node$isLeaf
    ## If the node has any attribute
    cond2 <- any(vapply(attrs, function(x) !is.null(node[[x]]) , logical(1)))
    if (cond1 && cond2) {
        ## This conditional is not ASR. Just recalculates existing scores.
        # message(node$name)
        # message('just homogenize leaf')
        scores <- vector('double', length(attrs))
        for (i in seq_along(scores)) {
            if (is.null(node[[attrs[i]]])) {
                scores[[i]] <- 0
            } else {
                scores[[i]] <- node[[attrs[i]]]
            }
        }
        if (sum(scores) > 0) {
            scores <- scores / sum(scores)
        } else {
            scores[scores == 0] <- NA
        }
        names(scores) <- attrs
        for (i in seq_along(scores)) {
            node[[attrs[i]]] <- scores[[attrs[i]]]
        }
    } else if (!cond1 && cond2) {
        ## This conditional is not ASR. Just recalculate existing scores.
        # message(node$name)
        # message('Just homogenize parent node')
        scores <- vector('double', length(attrs))
        for (i in seq_along(scores)) {
            if (is.null(node[[attrs[i]]]) || is.na(node[[attrs[i]]])) {
                scores[[i]] <- 0
            } else {
                scores[[i]] <- node[[attrs[i]]]
            }
        }
        if (sum(scores) > 0) {
            scores <- scores / sum(scores)
        } else {
            scores[scores == 0] <- NA
        }
        names(scores) <- attrs
        for (i in seq_along(scores)) {
            node[[attrs[i]]] <- scores[[attrs[i]]]
        }
    } else if (!cond1 && !cond2) {
        ## This is the real ASR part
        # message(node$name)
        # message('Aggregate and homogenize parent node.')
        ## Aggregate
        for (i in seq_along(attrs)) {
            children <- names(node$children)
            scores <- vector('double', length(children))
            for (j in seq_along(children)) {
                if (is.null(node[[children[j]]][[attrs[i]]]) || is.na(node[[children[j]]][[attrs[i]]])) {
                    scores[[j]] <- 0
                } else {
                    scores[[j]] <- node[[children[j]]][[attrs[i]]]
                }
            }
            node[[attrs[i]]] <- sum(scores)
        }
        ## Homogenize
        scores <- vector('double', length(attrs))
        for (i in seq_along(scores)) {
            if (is.null(node[[attrs[i]]]) || is.na(node[[attrs[i]]])) {
                scores[[i]] <- 0
            } else {
                scores[[i]] <- node[[attrs[i]]]
            }
        }
        if (sum(scores) > 0) {
            scores <- scores / sum(scores)
        } else {
            scores[scores == 0] <- NA
        }
        names(scores) <- attrs
        for (i in seq_along(scores)) {
            node[[attrs[i]]] <- scores[[attrs[i]]]
            attr_evi <- sub('^(.*)__Score', '\\1__Evidence', attrs[i])
            if (all(is.na(scores))) {
                node[[attr_evi]] <- NA
            } else {
                node[[attr_evi]] <- 'asr'
            }
        }
    }
}















#' Inheritance / Downstream (logical)
#'
#' \code{inhDownstreamLogical} each node inherits an attribute from it's
#' ancestor if values are NA or NULL
#'
#' @param node A node
#'
#' @return Node, R6, data.tree.
#' @export
#'
inhDownstreamLogical <- function(node) {
    attrs <- grep('__Score', node$attributesAll, value = TRUE)
    for (i in seq_along(attrs)) {
        if (is.null(node[[attrs[i]]]) || is.na(node[[attrs[i]]])) {
            node[[attrs[i]]] <- node[['parent']][[attrs[i]]]
            attr_evi <- sub('^(.*)__Score', '\\1__Evidence', attrs[i])
            node[[attr_evi]] <- 'inh'
        }
    }
}

#' Inheritance / Downstream (numeric)
#'
#' \code{inhDownstreamNumeric} each node inherits an attribute from it's ancestor
#' if values are NA or NULL
#'
#' @param node A node
#'
#' @return Node, R6, data.tree.
#' @export
#'
inhDownstreamNumeric <- function(node) {
    attr_score <- grep('__Score', node$attributesAll, value = TRUE)
    attr_evi <- grep('__Evidence', node$attributesAll, value = TRUE)
    attr_val <- grep('__Attribute_value', node$attributesAll, value = TRUE)
    if (is.null(node[[attr_score]]) || is.na(node[[attr_score]])) {
        node[[attr_score]] <- node[['parent']][[attr_score]]
        node[[attr_val]] <- node[['parent']][[attr_val]]
        node[[attr_evi]] <- 'inh'
    }
}

#' Inheritance / Downstream (range)
#'
#' \code{inhDownstreamRange} each node inherits an attribute from it's ancestor
#' if values are NA or NULL
#'
#' @param node A node
#'
#' @return Node, R6, data.tree.
#' @export
#'
inhDownstreamRange <- function(node) {
    attr_score <- grep('__Score', node$attributesAll, value = TRUE)
    attr_evi <- grep('__Evidence', node$attributesAll, value = TRUE)
    attr_val_min <- grep('__Attribute_value_min', node$attributesAll, value = TRUE)
    attr_val_max <- grep('__Attribute_value_max', node$attributesAll, value = TRUE)
    if (is.null(node[[attr_score]]) || is.na(node[[attr_score]])) {
        node[[attr_score]] <- node[['parent']][[attr_score]]
        node[[attr_val_min]] <- node[['parent']][[attr_val_min]]
        node[[attr_val_max]] <- node[['parent']][[attr_val_max]]
        node[[attr_evi]] <- 'inh'
    }
}

