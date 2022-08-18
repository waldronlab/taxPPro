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
        df, cols = c('Taxon_name', 'NCBI_ID'), verbpse = FALSE
    ) {

    Taxon_name <- Attribute <- NULL

    df <- df |>
        dplyr::filter(Attribute != "")

    index1 <- which(duplicated(df[, cols]))
    index2 <- which(duplicated(df[, cols], fromLast = TRUE))
    index <- sort(c(index1, index2))

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


#' Filter datset for phylogenetic propagation
#'
#' \code{filter_dataset_for_propagation} filters rows and datasets that
#' have enough information to perfrom ASR and Inheritance propagation of
#' taxa annotations
#'
#' @param df A dataset from bugphyzz.
#' @param df_name A character string. The name of the dataset
#'
#' @importFrom magrittr %>%
#' @return A flitered version of the input df. If no rows are conserved, the
#' output is NULL with a warning.
#' @export
filter_dataset_for_propagation <- function(df, df_name = NULL) {

    ## Columns required for propagation
    columns_for_propagation <- c(
        'Taxon_name', 'NCBI_ID',
        'Rank',
        'Attribute', 'Attribute_value', 'Attribute_source',
        'Evidence', 'Frequency',  # TODO add confidence in source
        'Parent_NCBI_ID', 'Parent_name', 'Parent_rank'
    )

    columns_lgl <- columns_for_propagation %in% colnames(df)

    if (!all(columns_lgl))  {

        missing_cols <- columns_for_propagation[!columns_lgl]

        if (!is.null(df_name)) {
            stop(
                'These columns are required for propagation in dataset ',
                df_name, ':',
                paste(columns_for_propagation, collapse = ', '),
                '. The following columns are missing: ',
                paste(missing_cols, collapse = ', '),
                call. = FALSE
            )

        } else {
            stop(
                'These columns are required for propagation: ',
                paste(columns_for_propagation, collapse = ', '),
                '. The following columns are missing: ',
                paste(missing_cols, collapse = ', '),
                call. = FALSE
            )
        }

    }

    ## Filtering
    df <- df %>%
        dplyr::filter(
            !is.na(Taxon_name) | Taxon_name != 'unknown',
            !is.na(Rank),
            Attribute_value != FALSE,
            !is.na(Attribute_source),
            #Evidence ## Complex
            #Frequency ## Comples,
            !is.na(Parent_NCBI_ID), !is.na(Parent_name), !is.na(Parent_rank),
            Rank %in% .valid_ranks(),
            Parent_rank %in% .valid_parent_ranks()
        ) %>%
        dplyr::distinct()

    ## Check output

    if (nrow(df) > 1) {
        return(df)
    } else {
        if (is.null(df_name)) {
            warning(
                paste('Not enough information for propagation.'),
                call. = FALSE
            )
        } else if (!is.null(df_name)) {
            warning(
                paste(
                    'Not enough information for propagation in dataset ',
                    df_name
                )
            )

        }
        return(NULL)
    }

}


#' Convert confidence intervals to numeric scores
#'
#' \code{.ci_to_scores} converts the keywords in the `confidence_interval`
#' column of a bugphyzz dataset into numeric scores, which added as an
#' additional column named as `Score`.
#'
#' @param x  A dataset from bugphyzz, e.g., aerophilicity.
#'
#' @return A datafraame. The same dataframe with the addtional `Score` column.
#'
#' @export
#'
ci_to_scores <- function(x) {
    Frequency <- NULL
    x %>%
        dplyr::mutate(
            ## TODO Change later Confidence_interval by Frequency
            Confidence_interval = stringr::str_to_lower(Frequency),
            Score = dplyr::case_when(
                Confidence_interval == 'always' ~ 1,
                Confidence_interval == 'usually' ~ 0.8,
                Confidence_interval == 'sometimes' ~ 0.5,
                Confidence_interval == 'rarely' ~ 0.2,
                Confidence_interval == 'never' ~ 0,
                Confidence_interval == 'unknown' ~ NA_real_
            )
        ) %>%
        dplyr::filter(!is.na(Score))
}

.valid_ranks <- function() {
    c(
        'strain', 'species', 'genus', 'family', 'order', 'class', 'phylum',
        'superkingdom'
    )
}

.valid_parent_ranks <- function() {
    c(
        'species', 'genus', 'family'
    )
}

#' Get parents
#'
#' \code{get_parents} Get's the parents of a vector of valid NCBI IDs.
#'
#' @param x A vector of valid NCBI IDs
#'
#' @return A table with Parent_NCBI_ID, Parent_name, and Parent_rank
#' @export
#'
get_parents <- function(x) {

    classification <- taxizedb::classification(x)
    parents_df <- purrr::map(classification, ~ {
        colnames(.x) <- c('Parent_name', 'Parent_rank', 'Parent_NCBI_ID')
        valid_ranks <- .valid_ranks()
        .x <- .x[.x$Parent_rank %in% valid_ranks, ]
        n_rows <- nrow(.x)
        .x <- .x[-n_rows, ]
        utils::tail(.x, 1)
    }) %>%
        dplyr::bind_rows()
    parents_df[,c('Parent_NCBI_ID', 'Parent_name', 'Parent_rank')]
}


#' Upstream
#'
#' \code{upstream} propagates annotations upstream in the taxonomy
#' classification according to the NCBI.
#'
#' @param df A dataframe from bugphyzz.
#' @param rank1 Source rank. A character string.
#' @param rank2 Target rank. A character string
#'
#' @importFrom magrittr %>%
#'
#' @return A dataframe
#' @export
#'
upstream <- function(df, rank1, rank2) {

    rank1_df <- dplyr::filter(df, Rank == rank1) ## TODO Recalculate Frequency
    rank2_df <- dplyr::filter(df, Rank == rank2)

    if (nrow(rank2_df) > 0) {

        both_ranks <- intersect(rank1_df$Parent_NCBI_ID, rank2_df$NCBI_ID)

        rank1_df <- rank1_df |>
            dplyr::filter(!Parent_NCBI_ID %in% both_ranks) |>
            dplyr::select(-NCBI_ID, -Taxon_name, -Rank) |>
            dplyr::rename(
                NCBI_ID = Parent_NCBI_ID, Taxon_name = Parent_name,
                Rank = Parent_rank
            ) %>%
            dplyr::relocate(NCBI_ID, Taxon_name)
        new_rank1_df <- dplyr::bind_cols(rank1_df, get_parents(rank1_df$NCBI_ID))
        new_rank1_df$Parent_NCBI_ID <- as.integer(new_rank1_df$Parent_NCBI_ID)
        new_rank1_df$NCBI_ID <- as.character(new_rank1_df$Parent_NCBI_ID)
        new_rank1_df$Evidence <- 'INH-UP'

        return(dplyr::distinct(dplyr::bind_rows(rank2_df, new_rank1_df)))

    } else {

        rank1_df <- rank1_df |>
            # dplyr::filter(!Parent_NCBI_ID %in% both_ranks) |>
            dplyr::select(-NCBI_ID, -Taxon_name, -Rank) |>
            dplyr::rename(
                NCBI_ID = Parent_NCBI_ID, Taxon_name = Parent_name,
                Rank = Parent_rank
            ) %>%
            dplyr::relocate(NCBI_ID, Taxon_name)
        new_rank1_df <- dplyr::bind_cols(rank1_df, get_parents(rank1_df$NCBI_ID))
        new_rank1_df$Parent_NCBI_ID <- as.integer(new_rank1_df$Parent_NCBI_ID)
        new_rank1_df$NCBI_ID <- as.character(new_rank1_df$Parent_NCBI_ID)
        new_rank1_df$Evidence <- 'INH-UP'

        return(dplyr::distinct(new_rank1_df))
    }
}


downstream <- function(df, rank1, rank2) {

    df <- dplyr::filter(df, Rank == rank1)
    taxa <- split(df, df$NCBI_ID)

    output <- vector('list', length(taxa))
    for (i in seq_along(output)) {
        taxon_df <- taxa[[i]]
        taxid <- dplyr::pull(taxon_df, NCBI_ID)
        new_taxa <- taxizedb::children(taxid, db = 'ncbi')[[1]] |>
            dplyr::rename(NCBI_ID = id, Taxon_name = name, Rank = rank) |>
            dplyr::filter(Rank == rank2)
        z <- taxon_df |>
            dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
            dplyr::rename(
                Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name,
                Parent_rank = Rank
            )
        message(taxid)
        output[[i]] <- dplyr::bind_cols(new_taxa, z[rep(1, nrow(z)), ]) |>
            dplyr::mutate(Evidence = 'INH-DOWN')
    }
    return(output)


    # j <- par[[1]]
    # taxid <- dplyr::pull(j, NCBI_ID)
    # new_taxa <- taxizedb::children(taxid, db = 'ncbi')[[1]] |>
    #     dplyr::rename(NCBI_ID = id, Taxon_name = name, Rank = rank) |>
    #     dplyr::filter(Rank == 'species')
    # z <- j |>
    #     dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
    #     dplyr::rename(
    #         Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name, Parent_rank = Rank
    #     )
    #
    # dplyr::bind_cols(new_taxa, z[rep(1, 11), ]) |>
    #     dplyr::mutate(Evidence = 'INH-DOWN')

}



