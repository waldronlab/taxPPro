
#' Propagate annotations
#'
#' \code{propagate} propagates annotations.
#'
#' @param df A data frame from bugphyzz.
#' @param asr_method A character string. The method that should be used for
#' propagation. Option: mv, majority vote; tx, NCBI taxonomy; ph, phylogenetic
#' tree.
#'
#' @return A data frame with extended annotations.
#'
#' @export
#'
propagate <- function(df, asr_method = 'mv', prop) {

    df_filtered <- preSteps(df)

    if (!nrow(df_filtered)) {
        warning('NO propagation for this dataset')
        return(df)
    }

    if (asr_method == 'mv') {
        no_filtered <- df[!df$NCBI_ID %in% df_filtered$NCBI_ID,] |>
            dplyr::mutate(
                Parent_NCBI_ID = as.character(Parent_NCBI_ID),
                NCBI_ID = as.character(NCBI_ID)
            )

        if (prop == 'both') {
            propagated <- df_filtered |>
                upstream() |>
                downstream()
        } else if (prop == 'upstream') {
            propagated <- df_filtered |>
                upstream()
        } else if (prop == 'downstream') {
            propagated <- df_filtered |>
                downstream()
        }

        propagated <- propagated[!propagated$NCBI_ID %in% no_filtered$NCBI_ID,]

        output <- dplyr::bind_rows(no_filtered, propagated) |>
            dplyr::distinct()

        return(output)

    }

}

#' Propagate upstream
#'
#' \code{propagateUpstream}
#'
#' @param df A data frame.
#' @param max.tax.level A character string. Maximum taxonomic rank/level.
#'
#' @return A data frame with asr evidence.
#' @export
#'
propagateUpstream <- function(df, max.tax.level) {

    df_filtered <- preSteps(df)

    if (!nrow(df_filtered)) {
        warning('No propagation upstream.', call. = FALSE)
        return(df)
    }

    split_by_rank <- split(df_filtered, factor(df_filtered$Rank))

    valid_ranks <- .validRanks()
    valid_ranks <- valid_ranks[1:which(valid_ranks == max.tax.level)]

    for (i in seq_along(valid_ranks)) {

        current_rank <- valid_ranks[i]

        if (current_rank %in% names(split_by_rank)) {
            message('Current rank: ', current_rank)
            next_pos <- i + 1
            if (next_pos > length(valid_ranks))
                break()
            next_rank <- valid_ranks[next_pos]
            message('Getting next rank: ', next_rank)
            parent_scores <- getParentScores(split_by_rank[[current_rank]])
            parent_scores <- parent_scores |>
                dplyr::filter(.data$Rank == next_rank)

            if (next_rank %in% names(split_by_rank)) {

                split_by_rank[[next_rank]] <-
                    .replaceParents(split_by_rank[[next_rank]], parent_scores)

            } else {
                split_by_rank[[next_rank]] <- parent_scores
            }

        } else {
            next()
        }
    }

    split_by_rank |>
        dplyr::bind_rows() |>
        dplyr::distinct() |>
        as.data.frame()
}

propagateDownstream <- function(df) {

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
        removeDuplicates() |>
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

    ## downstream family to genus
    if ('family' %in% names(split_by_rank)) {
        message('Getting new genera with inh-tax. (Step 4 - downstream)')
        family <- split_by_rank$family
        new_genera <- getChildren(family$NCBI_ID) |>
            dplyr::filter(Rank == 'genus') |>
            dplyr::mutate(Evidence = 'inh-tax') |>
            dplyr::distinct()

        family <- family |>
            dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
            dplyr::rename(
                Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name,
                Parent_rank = Rank
            )

        new_genera <- dplyr::left_join(
            new_genera, family, by = 'Parent_NCBI_ID',
            suffix = c('', '.y')
        ) |>
            dplyr::select(-tidyselect::ends_with('.y'))

        if ('genus' %in% names(split_by_rank)) {
            genus <- split_by_rank$genus
            new_genera <-
                new_genera[!new_genera$NCBI_ID %in% genus$NCBI_ID,]
        }

    }

    ## downstream genus to species
    if ('genus' %in% names(split_by_rank)) {
        message('Getting new species with inh-tax. (Step 5 - downstream)')
        genus <- split_by_rank$genus
        new_species <- getChildren(genus$NCBI_ID) |>
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
        message('Getting new strains with inh-tax (Step 6 - downstream).')
        species <- split_by_rank$species
        new_strains <- getChildren(species$NCBI_ID) |>
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

    new_downstream <- list(new_genera, new_species, new_strains) |>
        purrr::discard(is.null) |>
        dplyr::bind_rows()

    new_downstream <-
        new_downstream[!new_downstream$NCBI_ID %in% df_filtered$NCBI_ID,]
    dplyr::bind_rows(df_filtered, new_downstream) |>
        dplyr::distinct()
}

# Helper functions --------------------------------------------------------

#' Get parents
#'
#' \code{getParents} get's the parent taxon of a vector of valid NCBI IDs.
#'
#' @param x A vector of valid NCBI taxids
#'
#' @return A data frame with Parent_NCBI_ID, Parent_name, and Parent_rank.
#' @export
#'
getParents <- function(x) {
    classification <- taxizedb::classification(x)
    classification <- classification[!is.na(classification)]
    classification <- purrr::map(classification, dplyr::distinct)
    parents_df <- purrr::map(classification, ~ {
        colnames(.x) <- c('Parent_name', 'Parent_rank', 'Parent_NCBI_ID')
        valid_ranks <- .validRanks()
        .x <- .x[.x$Parent_rank %in% valid_ranks, ]
        n_rows <- nrow(.x)
        .x <- .x[-n_rows, ]
        utils::tail(.x, 1)
    }) |>
        dplyr::bind_rows()
    parents_df$NCBI_ID <- names(classification)
    parents_df[,c('NCBI_ID', 'Parent_NCBI_ID', 'Parent_name', 'Parent_rank')]
}

#' Get children
#'
#' \code{getChildren} gets the immediate descendants of a taxon.
#'
#' @param x A vector of valid NCBI taxids.
#'
#' @return A data frame with all children taxa.
#'
#' @export
#'
getChildren <- function(x) {
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
#' \code{getParentScores} calculates the parent score based on its immediate
#' descendants.
#'
#' @param df A data frame.
#' @param wt whether use weights on the dplyr::counts function.
#'
#' @return A data frame.
#'
#' @export
#'
getParentScores <- function(df) {

    attr_val_col <- chooseColVal(df)

    asr_scores <- df |>
        dplyr::group_by(
            .data$Parent_name,
            .data$Parent_NCBI_ID,
            .data$Parent_rank
        ) |>
        tidyr::nest() |>
        dplyr::ungroup() |>
        dplyr::mutate(
            data2 = purrr::map(# 'data' is the name of the new column (nested)
                .data$data, ~ .calcParentScore(.x, attr_val_col)
            )
        ) |>
        dplyr::select(-.data$data) |>
        tidyr::unnest(.data$data2) |>
        dplyr::rename(
            Taxon_name = .data$Parent_name,
            NCBI_ID = .data$Parent_NCBI_ID,
            Rank = .data$Parent_rank
        ) |>
        dplyr::mutate(Evidence = 'asr')

    pos <- which(colnames(asr_scores) == 'attr')
    colnames(asr_scores)[pos] <- attr_val_col

    parents <- tibble::as_tibble(getParents(asr_scores$NCBI_ID))
    parents$Parent_NCBI_ID <- as.character(parents$Parent_NCBI_ID)
    asr_scores$NCBI_ID <- as.character(asr_scores$NCBI_ID)
    output <- dplyr::right_join(asr_scores, parents, by = 'NCBI_ID')

    if (attr_val_col == 'Attribute') {
        output$Attribute_value <- TRUE
    } else if (attr_val_col == 'Attribute_value') {
        output$Attribute <- unique(df$Attribute)
    }

    output$Frequency <- .scores2Freq(output$Score)

    return(output)

}

getChildrenScores <- NULL

# Helper functions --------------------------------------------------------

## Replace parents
.replaceParents <- function(df, parent_scores) {
    new_taxa_lgl <- !parent_scores$Taxon_name %in% df$Taxon_name
    new_taxa <- parent_scores[new_taxa_lgl,]

    if (!length(new_taxa)) {
        return(df)
    } else {
        return(dplyr::bind_rows(df, new_taxa))
    }
}

## Replace children
.replaceChildren <- NULL

## Function to calculate parent score
.calcParentScore <- function(df, attr_col) {

    pos <- which(colnames(df) == attr_col)
    colnames(df)[pos] <- 'attr'

    pos <- which(colnames(df) == 'Score')
    colnames(df)[pos] <- 'val'

    df_cat <- df |>
        dplyr::group_by(.data[['attr']]) |>
        dplyr::summarise(
            Attribute_source = paste0(.data$Attribute_source, collapse = '|'),
            Confidence_in_curation = paste0(
                .data$Confidence_in_curation, collapse = '|'
            ),
            Original_NCBI_ID = paste0(.data$NCBI_ID, collapse = '|'),
            Original_Taxon_name = paste0(.data$Taxon_name, collapse = '|'),
            Original_Score = paste0(.data$val, collapse = '|'),
            Original_Evidence = paste0(.data$Evidence, collapse = '|')
        ) |>
        dplyr::ungroup()

    new_scores <- df |> dplyr::mutate(
        n = dplyr::n(),
        total = sum(.data$val),
        prop = ifelse(n == 1, .data$val, .data$val / .data$total)
        ) |>
        dplyr::count(.data$attr, wt = .data$prop, name = 'Score') |>
        dplyr::mutate(Score = round(.data$Score, 1))
        # dplyr::filter(.data$Score >= 0.5) # TODO maybe do/don't apply filter
    dplyr::left_join(new_scores, df_cat, by = 'attr')
}

## Function for converting scores to frequency
## Input is a vector
.scores2Freq <- function(Score) {
     dplyr::case_when(
         Score == 1 ~ 'always',
         Score >= 0.7 & Score < 1 ~ 'usually',
         Score >= 0.4 & Score < 0.7 ~ 'sometimes',
         Score >= 0.1 & Score < 0.4 ~ 'rarely',
         Score < 0.1 ~ 'never'
         )
}

## Function with valid ranks for taxa
.validRanks <- function() {
    c('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum',
      'superkingdom')
}

## Function with valid ranks for parents
.validParentRanks <- function() c('species', 'genus', 'family')
