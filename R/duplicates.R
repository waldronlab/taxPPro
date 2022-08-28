
#' Get duplicates
#'
#' \code{get_duplicates} gets duplicated taxa in a bugphuzz dataset.
#' This information could be useful for identifying conflicts of annotations
#' between different sources.
#'
#' @param df A dataframe.
#' @param cols Columns to look for duplicates.
#'
#' @return A dataframe with duplicated rows. If no duplicates are found, this
#' function returns a NULL value with a message.
#'
#' @export
#'
get_duplicates <- function(
        df, cols = c('Taxon_name', 'NCBI_ID'), verbose = FALSE
) {

    Taxon_name <- Attribute <- NULL

    df <- df |>
        dplyr::filter(Attribute != "")

    index1 <- which(duplicated(df[, cols]))
    index2 <- which(duplicated(df[, cols], fromLast = TRUE))
    index <- unique(sort(c(index1, index2)))

    duplicated_df <- dplyr::arrange(df[index,], Taxon_name)

    if (nrow(duplicated_df) == 0) {
        if (verbose) {
            message('No duplicates were found.')
        }
        return(NULL)
    } else {
        return(duplicated_df)
    }
}

#' Remove taxa duplicates
#' \code{remove_taxa_duplicates} remove taxa that are duplicated
#'
#' @param df A dataframe imported from bugphyzz
#'
#' @return A dataframe without duplicated taxa
#' @export
#'
remove_taxa_duplicates <- function(df) {

    dup <- get_duplicates(df)

    if (is.null(dup))
        return(df)

    no_select_taxa <- dup |>
        dplyr::pull(Taxon_name) |>
        unique()
    return(dplyr::filter(df, !Taxon_name %in% no_select_taxa))
}



#' Get conflicts
#'
#' @param df A data frame from bugphzz
#'
#' @return a dataframe or NULL
#' @export
#'
get_conflicts <- function(df) {
    df <- get_duplicates(df)
    if (is.null(df)) {
        message('No duplicates or conflictes here.')
        return(NULL)
    }
    split <- split(df, factor(df$Taxon_name))
    n_sources <- purrr::map_int(split, ~ length(unique(.x$Attribute_source)))
    n_attributes <- purrr::map_int(split, ~ {
        if (is.logical(.x$Attribute_value)) {
            length(unique(.x$Attribute))
        } else {
            length(unique(.x$Attribute_value))
        }
    })

    output <- split[n_sources > 1 & n_attributes > 1]

    if(!length(output)) {
        message('No duplicates here.')
        return(NULL)
    }

    output <- output |>
        purrr::discard(is.null) |>
        purrr::map(bind_rows) |>
        dplyr::bind_rows()

    return(output)
}


#' Resolve conflicts in a bugphyzz dataset
#'
#' @param df  A dataframe imported with bugphyzz
#'
#' @return A dataframe with resolved conflicts
#' @export
#'
resolve_conflicts <- function(df) {

    conflicts <- get_conflicts(df)

    ## if there are no conflicts, return the same dataset
    if (is.null(conflicts))
        return(df)

    ## Resolve conflicts for character/numerical values only
    if (is.logical(df$Attribute_value)) { ## Maybe change here with data type

        conflicts$Confidence_in_curation <- factor(
            x = conflicts$Confidence_in_curation,
            levels = c('Low', 'Medium', 'High'),
            ordered = TRUE
        )

        conflict_names <- unique(conflicts$Taxon_name)
        df_no_conflicts <- df |>
            dplyr::filter(!Taxon_name %in% conflict_names)

        resolved_conflicts <- conflicts |>
            dplyr::group_by(Taxon_name) |>
            dplyr::slice_max(Confidence_in_curation, with_ties = TRUE)

        conf_split <-
            split(resolved_conflicts, factor(resolved_conflicts$Taxon_name))

        n_rows <- map_int(conf_split, nrow)

        if (all(n_rows >= 2)) {
            msg <- paste0(
                'There were conflicts but none could be solved.',
                'Dropping ', length(n_rows), ' taxa.'
            )
            warning(msg, call. = FALSE)
            return(no_conflicts)
        }

        if (any(n_rows >= 2)) {
            remaining_conflicts <- sum(n_rows >= 2)
            msg <- paste0(
                remaining_conflicts,
                " conflicts couldn't be solved. Dropping them."
            )
            warning(msg, call. = FALSE)
        }

        resolved_conflicts <- conf_split[!n_rows >= 2] |>
            dplyr::bind_rows()

        output <-
            dplyr::bind_rows(df_no_conflicts, resolved_conflicts) |>
            dplyr::mutate(
                Confidence_in_curation = as.character(Confidence_in_curation)
            )

        return(output)

    } else {

        warning(
            'Cannot resolve conflicts of numerics and ranges yet.',
            ' No action taken.', call. = FALSE
        )
        return(df)
    }
}
