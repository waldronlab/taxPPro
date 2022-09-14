
#' Filter dataset propagation
#'
#' \code{filter_dataset_for_propagation} filters rows that
#' have enough information to perform ASR and Inheritance propagation of
#' taxa annotations.
#'
#' @param df A dataset from bugphyzz.
#' @param df_name A character string. The name of the dataset.
#'
#' @return A filtered version of the input dataframe.
#' If no rows are kept, the output is NULL with a warning.
#'
#' @export
#'
filter_dataset_for_propagation <- function(df, df_name = NULL) {

    ## Columns required for propagation
    columns_for_propagation <- c(
        'Taxon_name', 'NCBI_ID', 'Rank',
        'Attribute', 'Attribute_value', 'Attribute_source',
        'Evidence', 'Frequency', 'Confidence_in_curation',
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
    df <- df |>
        dplyr::filter(
            !is.na(Taxon_name) | Taxon_name != 'unknown',
            !is.na(Rank),
            Attribute_value != FALSE,
            !is.na(Attribute_source),
            ## Evidence # Complex
            ## Frequency # Complex
            !is.na(Parent_NCBI_ID), !is.na(Parent_name), !is.na(Parent_rank),
            Rank %in% .valid_ranks(),
            Parent_rank %in% .valid_parent_ranks()
            ## Should propagation to ranks higher than genus be included?
        ) |>
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

#' Get parents
#'
#' \code{get_parents} get's the parent taxon of a vector of valid NCBI IDs.
#'
#' @param x A vector of valid NCBI IDs
#'
#' @return A table with Parent_NCBI_ID, Parent_name, and Parent_rank.
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

#' Get children
#'
#' \code{get_children} gets the immediate descendants of a taxon.
#'
#' @param x A vector of NCBI IDs
#'
#' @return A data frame with all children taxa
#'
#' @export
#'
get_children <- function(x) {
    id <- name <- rank <- NULL
    taxizedb::children(x, db = 'ncbi') |>
        purrr::map(~ tibble::as_tibble(.x)) |>
        dplyr::bind_rows(.id = 'Parent_NCBI_ID') |>
        dplyr::rename(
            NCBI_ID = id, Taxon_name = name, Rank = rank
        )
}

#' Calculate parent score
#'
#' \code{calcParentScore} calculates the parent score based on its immediate
#' descendants.
#'
#' @param df A data frame.
#' @param wt whether use weights on the dplyr::counts function
#'
#' @return A data frame.
#'
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

#' Annotate parent nodes (upstream)
#'
#' \code{upstream} annotates parent node based on their immediate descendants
#' (children).
#'
#' @param df A data frame from bugphyzz.
#'
#' @return A data frame with parent nodes annotations (new rows).
#'
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

#' Get annotations for descendants (downstream)
#'
#' \code{downstream} annotates the descendants fo parent nodes.
#'
#' @param df A data frame from bugphyzz.
#'
#' @return A data frame with child nodes annotated (new rows).
#'
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

    new_downstream <-
        new_downstream[!new_downstream$NCBI_ID %in% df_filtered$NCBI_ID,]
    dplyr::bind_rows(df_filtered, new_downstream) |>
        dplyr::distinct()
}

#' Propagate annotations based on custom ASR
#'
#' \code{propagate} propagates annotations upstrean and downstream the
#' taxonomy classification of the taxa.
#'
#' @param df A data frame from bugphyzz.
#'
#' @return A data frame with extended annotations.
#'
#' @export
#'
propagate <- function(df) {

    df_filtered <- df |>
        resolve_conflicts() |>
        resolve_agreements() |>
        filter_dataset_for_propagation() |>
        remove_taxa_duplicates() |>
        ci_to_scores() |>
        dplyr::distinct() |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )

    if (!nrow(df_filtered)) {
        warning('NO propagation for this dataset')
        return(df)
    }

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
    x |>
        dplyr::mutate(
            Frequency = stringr::str_to_lower(Frequency),
            Score = dplyr::case_when(
                Frequency == 'always' ~ 1,
                Frequency == 'usually' ~ 0.8,
                Frequency == 'sometimes' ~ 0.5,
                Frequency == 'rarely' ~ 0.2,
                Frequency == 'never' ~ 0,
                Frequency == 'unknown' ~ NA_real_
            )
        ) |>
        dplyr::filter(!is.na(Score))
}

scores_to_ci <- function(x) {
    Score <- NULL
    x |>
        dplyr::mutate(
            Frequency = dplyr::case_when(
                Score ~ "always"

            )
        )
    ## TODO
}

.valid_ranks <- function() {
    c('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum',
      'superkingdom')
}

.valid_parent_ranks <- function() c('species', 'genus', 'family')
