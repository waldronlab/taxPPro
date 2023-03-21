
# Get functions -----------------------------------------------------------

#' Get duplicates
#'
#' \code{getDuplicates} gets duplicated taxa in a bugphuzz dataset.
#' This information could be useful for identifying conflicts of annotations
#' between different sources. Duplication is based on both NCBI_ID and
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
        df, cols = 'NCBI_ID', verbose = FALSE
) {
    if (verbose)
        message('Looking for duplicates.')
    index1 <- which(duplicated(df[, cols]))
    if (!length(index1)) {
        if (verbose)
            message('No duplicates were found. Returning NULL')
        return(NULL)
    }
    index2 <- which(duplicated(df[, cols], fromLast = TRUE))
    index <- sort(unique(c(index1, index2)))
    dplyr::arrange(df[index,], .data$NCBI_ID)
}

#' Get agreements
#'
#' \code{getAgreements} gets taxa annotated two or more times with the same
#' annotation, but from different sources.
#'
#' An agreement is defined as taxon with the same annotation from different
#' sources.
#'
#' @param df A data frame.
#'
#' @return A data frame.
#' @export
#'
#' @seealso
#' \code{\link{resolveAgreements}}
#'
getAgreements <- function(df) {

    dup <- getDuplicates(df)

    if (is.null(dup)) {
        message('No duplicates or agreements here.')
        return(NULL)
    }

    attr_col <- chooseColVal(dup)

    agreements <- dup |>
        dplyr::count(
            .data$Taxon_name,
            .data[[attr_col]]
        ) |>
        dplyr::filter(.data$n > 1)


    if (!nrow(agreements)) {
        message("No agreements detected.")
        return(NULL)
    }

    tax_names <- agreements$Taxon_name
    attr_vals <- agreements[[attr_col]]

    dup |>
        dplyr::filter(
            .data$Taxon_name %in% tax_names,
            .data[[attr_col]] %in% attr_vals
        )

}

#' Get conflicts
#'
#' \code{getConflicts} gets taxa with two or more annotations from different
#' sources.
#'
#' @param df A data frame. Output of the `getDuplicates` function.
#'
#' @return a data frame or NULL.
#' @export
#'
getConflicts <- function(df) {

    dup <- getDuplicates(df)

    if (is.null(dup)) {
        message('No duplicates or conflicts here.')
        return(NULL)
    }

    split <- split(dup, factor(dup$Taxon_name))
    n_sources <- purrr::map_int(split, ~ length(unique(.x$Attribute_source)))
    n_attributes <- purrr::map_int(split, ~ {
        attr_col <- chooseColVal(.x)
        length(unique(.x[[attr_col]]))
    })

    conflicts <- split[n_sources > 1 & n_attributes > 1]

    if(!length(conflicts)) {
        message('No conflicts here.')
        return(NULL)
    }

    conflicts |>
        purrr::discard(is.null) |>
        purrr::map(~ dplyr::bind_rows(.x)) |>
        dplyr::bind_rows()

}

#' Get taxa with double annotations
#'
#' \code{getDoubleAnnotations} gets taxa annotated two or more times from the
#' same source.
#'
#' @param df A data frame.
#'
#' @return A data frame.
#' @export
#'
#' @seealso
#' \code{\link{getDuplicates}}
#'
getDoubleAnnotations <- function(df) {

    dups <- getDuplicates(df)

    if (is.null(dups)) {
        message('No duplicates or double annotations here.')
        return(NULL)
    }

    taxa_counts <- dups |>
        dplyr::count(.data$Attribute_source, .data$Taxon_name) |>
        dplyr::arrange(.data$Taxon_name, .data$Attribute_source) |>
        dplyr::filter(.data$n > 1)

    if (!nrow(taxa_counts)) {
        message('No double annotations in this data frame.')
        return(NULL)
    }

    taxa <- unique(taxa_counts$Taxon_name)
    dups[dups$Taxon_name %in% taxa,]
}

# Action functions --------------------------------------------------------

#' Resolve agreements
#'
#' \code{resolveAgreements} resolves agreements.
#'
#' Agreements are defined as a taxon with the same annotation from two or more
#' sources. The agreements are resolved with
#' `dplyr::slice_max(x, with_ties = FALSE)`. So, only one source is kept,
#' the one with the highest 'confidence in curation' value.
#'
#' @param df A data frame imported from bugphuyzz.
#'
#' @return A data frame
#' @export
#'
#' @seealso
#' \code{\link{getAgreements}}
#'
resolveAgreements <- function(df) {

    agreements <- getAgreements(df)
    if (is.null(agreements)) {
        message('No agreements to resolve.')
        return(df)
    }
    attr_col <- chooseColVal(agreements)
    agreements$Confidence_in_curation <- factor(
        x = agreements$Confidence_in_curation,
        levels = c('low', 'medium', 'high'),
        ordered = TRUE
    )
    tax_ids <- agreements$NCBI_ID
    attr_vals <- agreements[[attr_col]]
    index <-
        which(df$NCBI_ID %in% taxids & df[[attr_col]] %in% attr_vals)
    new_df <- df[-index,]
    resolved_agreements <- agreements |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
            order_by = dplyr::desc(.data$Confidence_in_curation),
            with_ties = FALSE
        ) |>
        dplyr::mutate(
            Confidence_in_curation = as.character(.data$Confidence_in_curation)
        )
    dplyr::bind_rows(new_df, resolved_agreements)
}

#' Resolve conflicts in a bugphyzz dataset
#'
#' \code{resolveConflicts} resolves conflicts of categorical values.
#'
#' A conflict is defined as a taxon with two different annotations from two
#' different attribute sources.
#'
#' The conflicts are resolved by giving preference to the source with higher
#' 'confidence in curation' value. This is done with `dplyr::slice_max`. Ties
#' are resolved according to the output of `dplyr::slice_max`.
#'
#' @param df  A dataframe imported with bugphyzz
#'
#' @return A dataframe with resolved conflicts
#' @export
#'
resolveConflicts <- function(df) {

    conflicts <- getConflicts(df)

    if (is.null(conflicts))
        return(df)

    attr_col <- chooseColVal(conflicts)

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
        dplyr::slice_max(
            order_by = dplyr::desc(.data$Confidence_in_curation),
            with_ties = TRUE
        ) |>
        dplyr::ungroup()

    conf_split <- split(
        x = resolved_conflicts, f = factor(resolved_conflicts$Taxon_name)
    )

    unresolved_lgl <- conf_split |>
        purrr::map_lgl(~ length(unique(.x$Attribute_source)) > 1)

    total_unresolved_conflicts <- sum(unresolved_lgl)

    if (all(unresolved_lgl)) {
        msg <- paste0(
            'There were ', length(conflict_names),
            ' conflicts but none could be solved.'
            # ' Dropping ', total_unresolved_conflicts, ' taxa.'
        )
        warning(msg, call. = FALSE)
        return(df)
    }

    if (any(unresolved_lgl)) {
        msg <- paste0(
            'There were ', length(conflict_names),
            ' conflicts, but ', total_unresolved_conflicts,
            " conflicts couldn't be solved."
        )
        warning(msg, call. = FALSE)
    }

    resolved_conflicts <- conf_split[!unresolved_lgl] |>
        dplyr::bind_rows()

    unresolved_conflicts <- conf_split[unresolved_lgl] |>
        dplyr::bind_rows()

    df_no_conflicts |>
        dplyr::bind_rows(resolved_conflicts) |>
        dplyr::bind_rows(unresolved_conflicts) |>
        dplyr::mutate(
            Confidence_in_curation = as.character(
                .data$Confidence_in_curation
            )
        )
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
#' @seealso
#' \code{\link{getDuplicates}}
#'
removeDuplicates <- function(df) {
    dup <- getDuplicates(df)
    if (is.null(dup))
        return(df)
    df[!df$Taxon_name %in% unique(dup$Taxon_name),]
}

# helper functions --------------------------------------------------------

#' Choose column for attribute values
#'
#' \code{chooseColVal} chooses the column where attribute values are
#' actually stored.
#'
#' @param df A data frame imported with bugphyzz.
#'
#' @return A character string.
#' @export
#'
chooseColVal <- function(df) {

    av <- df$Attribute_value

    if(is.logical(av)) {
        attr_col <- 'Attribute'
    } else if (is.numeric(av) || is.character(av)) {
        attr_col <- 'Attribute_value'
    }
    attr_col
}
