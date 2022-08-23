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
    classification <- purrr::map(classification, dplyr::distinct)
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

#' Get children taxa
#'
#' @param x A vector of NCBI IDs
#'
#' @return A dataframe with all children taxa
#' @export
#'
get_children <- function(x) {
    taxizedb::children(x, db = 'ncbi') |>
        purrr::map(~ tibble::as_tibble(.x)) |>
        dplyr::bind_rows(.id = 'Parent_NCBI_ID') |>
        dplyr::rename(
            NCBI_ID = id, Taxon_name = name, Rank = rank
        )
}

#' Calculate parent scores
#'
#' \code{calcParentScore} calculates the parent scores in a bugphyzz dataset.
#'
#' @param df A data frame.
#' @param wt If use weights on the dplyr::counts function
#'
#' @return A data frame.
#' @export
#'
calcParentScore <- function(df) {

    .calc_parent_score <- function(df, attr_val_col) {
        df <- df[, c(attr_val_col, 'Score')]
        attr <- val <- NULL
        colnames(df) <- c('attr', 'val')
        df |> dplyr::mutate(
            total = sum(val),
            prop = val / total
        ) |>
            dplyr::count(attr, wt = prop, name = 'Score') |>
            dplyr::mutate(Score = round(Score, 1)) |>
            dplyr::filter(Score >= 0.5)
    }

    Parent_name <- data <- data2 <- NULL


    if (is.logical(df$Attribute_value)) {
        attr_val_col <- 'Attribute'
    } else {
        attr_val_col <- 'Attribute_value'
    }

    df <- df |>
        dplyr::group_by(Parent_name, Parent_NCBI_ID, Parent_rank) |>
        tidyr::nest() |>
        dplyr::ungroup() |>
        dplyr::mutate(
            data2 = purrr::map(data, .calc_parent_score, attr_val_col)
        ) |>
        dplyr::select(-data) |>
        tidyr::unnest(data2) |>
        dplyr::rename(
            Taxon_name = Parent_name, NCBI_ID = Parent_NCBI_ID,
            Rank = Parent_rank
        ) |>
        dplyr::mutate(Evidence = 'asr-tax')

    pos <- which(colnames(df) == 'attr')
    colnames(df)[pos] <- attr_val_col

    parents <- tibble::as_tibble(get_parents(df$NCBI_ID))
    dplyr::bind_cols(df, parents) |>
        dplyr::mutate(
            NCBI_ID = as.character(NCBI_ID),
            Parent_NCBI_ID = as.character(Parent_NCBI_ID)
        )
}

#' Code used to annotate parent nodes according to the NCBI taxonomy
#'
#' @param df A data frame from bugphyzz
#'
#' @return A datta frame with parent nodes annotations
#' @export
#'
upstream <- function(df) {

    Parent_NCBI_ID <- NCBI_ID <- NULL

    df_filtered <- df |>
        filter_dataset_for_propagation() |>
        remove_taxa_duplicates() |>
        ci_to_scores() |>
        dplyr::distinct() |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )

    if (nrow(df) == 0) {
        warning('Nothing to do here.', call. = FALSE)
        return(df_filtered)
    }

    split_by_rank <- split(df_filtered, factor(df_filtered$Rank))

    if (!any(c('strain', 'species') %in% names(split_by_rank))) {
        warning('Nothing to do here.', call. = FALSE)
        return(df_filtered)
    }

    strains <- new_species <- new_genera <- NULL

    ## upstream species
    if ('strain' %in% names(split_by_rank)) {
        message('Getting new species with asr-tax. (Step 1 - upstream).')
        strains <- split_by_rank$strain
        new_species <- calcParentScore(strains) |>
            dplyr::filter(Rank == 'species') |>
            dplyr::distinct()
    }

    ## upstream genera
    if ('species' %in% names(split_by_rank)) {
        message('Getting new genera with asr-tax. (Step 2 = upstream).')
        species <- split_by_rank$species
        if (!is.null(new_species)) {
            new_species <- new_species[!new_species$NCBI_ID %in% species$NCBI_ID,]
            new_species <- dplyr::bind_rows(species, new_species) |>
                purrr::discard(~all(is.na(.x)))
        }
        new_genera <- calcParentScore(species) |>
            dplyr::filter(Rank == 'genus') |>
            dplyr::distinct()
    }

    new_upstream <- list(new_species, new_genera) |>
        purrr::discard(is.null) |>
        dplyr::bind_rows()

    new_upstream <- new_upstream[!new_upstream$NCBI_ID %in% df_filtered$NCBI_ID,]
    dplyr::bind_rows(df_filtered, new_upstream) |>
        dplyr::distinct()
}

#' Get annotations for children node (downstream)
#'
#' @param df A data frame from bugphyzz
#'
#' @return A data frame with children nodes annotated
#' @export
#'
downstream <- function(df) {
    Parent_NCBI_ID <- NCBI_ID <- NULL

    df_filtered <- df |>
        remove_taxa_duplicates() |>
        dplyr::distinct() |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )

    if (nrow(df) == 0) {
        warning('Nothing to do here.', call. = FALSE)
        return(NULL)
    }

    split_by_rank <- split(df_filtered, factor(df_filtered$Rank))

    if (!any(c('species', 'genus') %in% names(split_by_rank))) {
        warning('Nothing to do here.', call. = FALSE)
        return(NULL)
    }

    new_species <- new_strains <- NULL

    ## downstream genus to species
    if ('genus' %in% names(split_by_rank)) {
        message('Getting new species with inh-tax. (Step 3 - downstream)')
        genus <- split_by_rank$genus
        new_species <- get_children(genus$NCBI_ID) |>
            dplyr::filter(Rank == 'species') |>
            dplyr::mutate(Evidence = 'inh-tax') |>
            dplyr::distinct()

        genus <- genus |>
            dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
            dplyr::rename(
                Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name,
                Parent_rank = Rank
            )

        new_species <- dplyr::left_join(
            new_species, genus, by = 'Parent_NCBI_ID',
            suffix = c('', '.y')
        ) |>
            dplyr::select(-tidyselect::ends_with('.y'))

        if ('species' %in% names(split_by_rank)) {
            species <- split_by_rank$species
            new_species <-
                new_species[!new_species$NCBI_ID %in% species$NCBI_ID,]
        }

    }

    ## Downstream species to strain
    if ('species' %in% names(split_by_rank)) {
        message('Getting new strains with inh-tax (Step 4 - downstream).')
        species <- split_by_rank$species
        new_strains <- get_children(species$NCBI_ID) |>
            dplyr::filter(Rank == 'strain') |>
            dplyr::mutate(Evidence = 'inh-tax') |>
            dplyr::distinct()

        species <- species |>
            dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
            dplyr::rename(
                Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name,
                Parent_rank = Rank
            )

        new_strains <- dplyr::left_join(
            new_strains, species, by = 'Parent_NCBI_ID',
            suffix = c('', '.y')
        ) |>
            dplyr::select(-tidyselect::ends_with('.y'))

        if ('strain' %in% names(split_by_rank)) {
            strain <- split_by_rank$strain
            new_strains <-
                new_strains[!new_strains$NCBI_ID %in% strain$NCBI_ID,]
        }
    }

    new_downstream <- list(new_species, new_strains) |>
        purrr::discard(is.null) |>
        dplyr::bind_rows()

    new_downstream <- new_downstream[!new_downstream$NCBI_ID %in% df_filtered$NCBI_ID,]
    dplyr::bind_rows(df_filtered, new_downstream) |>
        dplyr::distinct()
}


#' Propagate annotations
#'
#' @param df A data frame from bugphyzz.
#'
#' @return A data frame with extended annotations.
#' @export
#'
propagate <- function(df) {

    df_filtered <- df |>
        filter_dataset_for_propagation() |>
        remove_taxa_duplicates() |>
        ci_to_scores() |>
        dplyr::distinct() |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )

    no_filtered <- df[!df$NCBI_ID %in% df_filtered$NCBI_ID,] |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )

    propagated <- df_filtered |>
        upstream() |>
        downstream()

    propagated <- propagated[!propagated$NCBI_ID %in% no_filtered$NCBI_ID,]

    dplyr::bind_rows(no_filtered, propagated) |>
        dplyr::distinct()

}


# upstream <- function(df, rank1, rank2) {
#
#     rank1_df <- dplyr::filter(df, Rank == rank1) ## TODO Recalculate Frequency with mean
#     rank2_df <- dplyr::filter(df, Rank == rank2)
#
#     if (nrow(rank2_df) > 0) {
#
#         both_ranks <- intersect(rank1_df$Parent_NCBI_ID, rank2_df$NCBI_ID)
#
#         rank1_df <- rank1_df |>
#             dplyr::filter(!Parent_NCBI_ID %in% both_ranks) |>
#             dplyr::select(-NCBI_ID, -Taxon_name, -Rank) |>
#             dplyr::rename(
#                 NCBI_ID = Parent_NCBI_ID, Taxon_name = Parent_name,
#                 Rank = Parent_rank
#             ) %>%
#             dplyr::relocate(NCBI_ID, Taxon_name)
#         new_rank1_df <- dplyr::bind_cols(rank1_df, get_parents(rank1_df$NCBI_ID))
#         new_rank1_df$Parent_NCBI_ID <- as.integer(new_rank1_df$Parent_NCBI_ID)
#         new_rank1_df$NCBI_ID <- as.character(new_rank1_df$Parent_NCBI_ID)
#         new_rank1_df$Evidence <- 'INH-UP'
#
#         return(dplyr::distinct(dplyr::bind_rows(rank2_df, new_rank1_df)))
#
#     } else {
#
#         rank1_df <- rank1_df |>
#             # dplyr::filter(!Parent_NCBI_ID %in% both_ranks) |>
#             dplyr::select(-NCBI_ID, -Taxon_name, -Rank) |>
#             dplyr::rename(
#                 NCBI_ID = Parent_NCBI_ID, Taxon_name = Parent_name,
#                 Rank = Parent_rank
#             ) %>%
#             dplyr::relocate(NCBI_ID, Taxon_name)
#         new_rank1_df <- dplyr::bind_cols(rank1_df, get_parents(rank1_df$NCBI_ID))
#         new_rank1_df$Parent_NCBI_ID <- as.integer(new_rank1_df$Parent_NCBI_ID)
#         new_rank1_df$NCBI_ID <- as.character(new_rank1_df$Parent_NCBI_ID)
#         new_rank1_df$Evidence <- 'INH-UP'
#
#         return(dplyr::distinct(new_rank1_df))
#     }
# }
#
#
# downstream <- function(df, rank1, rank2) {
#
#     df <- dplyr::filter(df, Rank == rank1)
#     taxa <- split(df, df$NCBI_ID)
#
#     output <- vector('list', length(taxa))
#     for (i in seq_along(output)) {
#         taxon_df <- taxa[[i]]
#         taxid <- dplyr::pull(taxon_df, NCBI_ID)
#         new_taxa <- taxizedb::children(taxid, db = 'ncbi')[[1]] |>
#             dplyr::rename(NCBI_ID = id, Taxon_name = name, Rank = rank) |>
#             dplyr::filter(Rank == rank2)
#         z <- taxon_df |>
#             dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
#             dplyr::rename(
#                 Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name,
#                 Parent_rank = Rank
#             )
#         message(taxid)
#         output[[i]] <- dplyr::bind_cols(new_taxa, z[rep(1, nrow(z)), ]) |>
#             dplyr::mutate(Evidence = 'INH-DOWN')
#     }
#     return(output)
# }


