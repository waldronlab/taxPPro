
## A function for removing prefix from a character vector x
removePrefix <- function(x) {
    regex <- '^[kpcofgst]__'
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
        children <- df$id
        new_children <- children[!children %in% current_children]
        new_children <- paste0('t__', new_children)
        if (isFALSE(!length(new_children))) {
            for (i in seq_along(new_children)) {
                message(new_children)
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
    cond2 <- new_ncbi_taxids$superkingdom == '2'
    bacteria <- new_ncbi_taxids[cond1 & cond2,]
    no_cols <- c('species', 'strain', 'Taxon_name', 'Parent_NCBI_ID', 'Rank')
    bacteria <- bacteria[,!colnames(bacteria) %in% no_cols]
    # for (i in seq_along(bacteria)) {
    #   bacteria[[i]] <- fillNAs(bacteria[[i]])
    # }
    bacteria$NCBI_ID <- paste0('s__', bacteria$NCBI_ID)
    bacteria <- tidyr::drop_na(bacteria)
    if (format == 'table') {
        return(bacteria)
    } else if (format == 'tree') {
        pathString <- paste(
            'k__', bacteria$superkingdom,
            '|||p__', bacteria$phylum,
            '|||c__', bacteria$class,
            '|||o__', bacteria$order,
            '|||f__', bacteria$family,
            '|||g__', bacteria$genus,
            '|||s__', bacteria$NCBI_ID,
            sep = ''
        )
        bacteria$pathString <- pathString
        tree <- data.tree::as.Node(bacteria, pathDelimiter = '|||')
        return(tree)
    }
}
