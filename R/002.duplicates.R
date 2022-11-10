
#' Get duplicates
#'
#' \code{getDuplicates} gets duplicated taxa in a bugphuzz dataset.
#' This information could be useful for identifying conflicts of annotations
#' between different sources. Dupliation is based on both NCBI_ID and
#' Taxon_name.
#'
#' @param df A data frame imported with bugphyzz.
#' @param cols Columns to look for duplicates.
#'
#' @return A data frame with duplicate rows. If no duplicates are found, this
#' function returns a NULL value with a message.
#'
#' @export
#'
getDuplicates <- function(
        df, cols = c('Taxon_name', 'NCBI_ID'), verbose = FALSE
) {

    if (verbose)
        message('Looking for duplicates.')

    df <- df[!is.na(df$Taxon_name) | df$Taxon_name != 'unknown',]
    df <- df[df$Attribute != '' & df$Attribute_value != FALSE,]

    index1 <- which(duplicated(df[, cols]))

    if (!length(index1)) {
        if (verbose)
            message('No duplicates were found. Returning NULL')
        return(NULL)
    }

    index2 <- which(duplicated(df[, cols], fromLast = TRUE))
    index <- sort(unique(c(index1, index2)))

    dplyr::arrange(df[index,], .data$Taxon_name)
}

#' Get taxa with double annotations
#'
#' \code{getDoubleAnnotations} gets taxa annotated twice from the same source.
#'
#' @param df A data frame.
#'
#' @return A data frame.
#' @export
#'
getDoubleAnnotations <- function(df) {

    dups <- getDuplicates(df)

    if (is.null(dups)) {
        message('No duplicates or conflictes here.')
        return(NULL)
    }

    ## Counts are based on Taxon_name, because NCBI_ID might be missing
    taxa_counts <- dups |>
        dplyr::count(.data$Attribute_source, .data$Taxon_name) |>
        dplyr::arrange(.data$Taxon_name, .data$Attribute_source) |>
        dplyr::filter(.data$n > 1)

    if (!nrow(taxa_counts)) {
        message('No double annotations in this data frame.')
        return(NULL)
    }

    taxa <- unique(taxa_counts$Taxon_name)
    double_annotations_all <- dups[dups$Taxon_name %in% taxa,]

    ## Remove conflicts
    conflicts <- getConflicts(df)
    if (!is.null(conflicts)) {
        conflict_taxa <- unique(conflicts$Taxon_name)
    } else {
        conflict_taxa <- NULL
    }

    ## Remove agreements
    agree <- getAgreements(df)
    if (is.null(conflicts)) {
        agree_taxa <- unique(agree$Taxon_name)
    } else {
        agree_taxa <- NULL
    }

    remove_taxa <- c(conflict_taxa, agree_taxa)
    if (is.null(remove_taxa)) {
        message('No double annotations')
        return(NULL)
    }

    double_annotations_all[!double_annotations_all$Taxon_name %in% remove_taxa,]
}

#' Get agreements
#'
#' \code{getAgreements} gets taxa annotated two or more times from different
#' sources.
#'
#' @param df A data frame.
#'
#' @return A data frame.
#' @export
#'
getAgreements <- function(df) {

    df <- getDuplicates(df)

    if (is.null(df)) {
        message('No duplicates or conflictes here.')
        return(NULL)
    }

    if(is.logical(df$Attribute_value)) {
        attr_col <- 'Attribute'
    } else if (is.numeric(df$Attribute_value)) {
        attr_col <- 'Attribute_value'
    }

    split <- split(df, df$Taxon_name)

    select_lgl <- split |>
        purrr::map_lgl(~ {
             attrs <- unique(.x[[attr_col]])
             sources <- unique(.x[['Attribute_source']])
             lgl_value <- length(attrs) == 1 & length(sources) > 1
             lgl_value
        })

    if (!any(select_lgl)) {
        message("No agreements.")
        return(NULL)
    }

    agreements_df <- split[select_lgl] |>
        dplyr::bind_rows()

    conflicts_df <- getConflicts(df)

    if (is.null(conflicts_df))
        return(agreements_df)

    agree_taxa <- unique(agreements_df$Taxon_name)
    conflict_taxa <- unique(conflicts_df$Taxon_name)

    real_agree_taxa <- agree_taxa[!agree_taxa %in% conflict_taxa]
    agreements_df |>
        dplyr::filter(Taxon_name %in% real_agree_taxa)
}

#' Resolve conflicts
#'
#' \code{resolve_agreements} resolves agreements, returning only one
#' (the highest).
#'
#' @param df A data frame imported from bugphuyzz.
#'
#' @return A data frame
#' @export
#'
resolve_agreements <- function(df) {

    agree_df <- getAgreements(df)

    if (is.null(agree_df)) {
        message('No agreements to solve')
        return(df)
    }

    agree_names <- unique(agree_df$Taxon_name)
    df_no_agreements <- df |>
        dplyr::filter(!Taxon_name %in% agree_names)

    resolved_agreements <- agree_df |>
        dplyr::group_by(Taxon_name) |>
        dplyr::slice_max(Confidence_in_curation, with_ties = FALSE)

    dplyr::bind_rows(df_no_agreements, resolved_agreements) |>
        dplyr::mutate(
            Confidence_in_curation = as.character(Confidence_in_curation)
        )

}

#' Get conflicts
#'
#' \code{getConflicts} gets taxa with two or more annotations from different
#' sources.
#'
#' @param df A data frame from bugphzz
#'
#' @return a data frame or NULL
#' @export
#'
getConflicts <- function(df) {

    df <- getDuplicates(df)

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
        message('No conflicts here.')
        return(NULL)
    }

    output <- output |>
        purrr::discard(is.null) |>
        purrr::map(~ dplyr::bind_rows(.x)) |>
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

    conflicts <- getConflicts(df)

    ## if there are no conflicts, return the same dataset
    if (is.null(conflicts))
        return(df)

    ## Resolve conflicts for character/numerical values only
    if (is.logical(df$Attribute_value)) { ## Maybe change here with data type

        conflicts$Confidence_in_curation <- factor(
            x = conflicts$Confidence_in_curation,
            levels = c('low', 'medium', 'high'),
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

        n_rows <- purrr::map_int(conf_split, nrow)

        if (all(n_rows >= 2)) {
            msg <- paste0(
                'There were conflicts but none could be solved.',
                'Dropping ', length(n_rows), ' taxa.'
            )
            warning(msg, call. = FALSE)
            return(df_no_conflicts)
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

#' Remove taxa duplicates
#' \code{removeDuplicates} removes taxa that are duplicated (output of
#' the `getDuplicates` function).
#'
#' @param df A data frame imported with bugphyzz.
#'
#' @return A data frame without duplicated taxa.
#' @export
#'
removeDuplicates <- function(df) {

    dup <- getDuplicates(df)

    if (is.null(dup))
        return(df)

    no_select_taxa <- unique(dup$Taxon_name)
    df[!df$Taxon_name %in% no_select_taxa,]
}
