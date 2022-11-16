
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

    df_filtered <- preSteps(df)

    if (nrow(df) == 0) {
        warning('Nothing to do here.', call. = FALSE)
        return(df_filtered)
    }

    split_by_rank <- split(df_filtered, factor(df_filtered$Rank))

    if (!any(c('strain', 'species', 'genus') %in% names(split_by_rank))) {
        warning('Nothing to do here.', call. = FALSE)
        return(df_filtered)
    }

    strains <- new_species <- new_genera <- NULL

    ## upstream species
    if ('strain' %in% names(split_by_rank)) {
        message('Getting new species with asr-tax. (Step 1 - upstream).')
        strains <- split_by_rank$strain
        new_species <- getParentScores(strains) |>
            dplyr::filter(.data$Rank == 'species') |>
            dplyr::distinct()
    }

    ## upstream genera
    if ('species' %in% names(split_by_rank)) {
        message('Getting new genera with asr-tax. (Step 2 - upstream).')
        species <- split_by_rank$species
        if (!is.null(new_species)) {
            new_species <- new_species[!new_species$NCBI_ID %in% species$NCBI_ID,]
            new_species <- dplyr::bind_rows(species, new_species) |>
                purrr::discard(~all(is.na(.x)))
        }
        new_genera <- getParentScores(species) |>
            dplyr::filter(Rank == 'genus') |>
            dplyr::distinct()
    }

    ## Add code for upstream family
    if ('genus' %in% names(split_by_rank)) {
        message('Getting new families with asr-tax. (Step 3 - upstream).')
        genera <- split_by_rank$genus
        if (!is.null(new_genera)) {
            new_genera <- new_genera[!new_genera$NCBI_ID %in% genera$NCBI_ID,]
            new_genera <- dplyr::bind_rows(genera, new_genera) |>
                purrr::discard(~all(is.na(.x)))
        }
        new_families <- getParentScores(genera) |>
            dplyr::filter(Rank == 'family') |>
            dplyr::distinct()
    }

    new_upstream <- list(new_species, new_genera, new_families) |>
        purrr::discard(is.null) |>
        dplyr::bind_rows()

    new_upstream <- new_upstream[!new_upstream$NCBI_ID %in% df_filtered$NCBI_ID,]
    output <- dplyr::bind_rows(df_filtered, new_upstream) |>
        dplyr::distinct()
    return(output)
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
        new_genera <- get_children(family$NCBI_ID) |>
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
        message('Getting new strains with inh-tax (Step 6 - downstream).')
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
#' \code{get_parents} get's the parent taxon of a vector of valid NCBI IDs.
#'
#' @param x A vector of valid NCBI taxids
#'
#' @return A data frame with Parent_NCBI_ID, Parent_name, and Parent_rank.
#' @export
#'
get_parents <- function(x) {
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
    }) %>%
        dplyr::bind_rows()
    parents_df$NCBI_ID <- names(classification)
    parents_df[,c('NCBI_ID', 'Parent_NCBI_ID', 'Parent_name', 'Parent_rank')]
}

#' Get children
#'
#' \code{get_children} gets the immediate descendants of a taxon.
#'
#' @param x A vector of valid NCBI taxids.
#'
#' @return A data frame with all children taxa.
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

    ## TODO not necessarily in this line, but find a way to concatenate
    ## the sources and original score values.

    ## TODO also add concatenated version of confidence in curation.

    ## TODO maybe add Accession_ID and Genome_ID

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
            data2 = purrr::map(
                ## 'data' is the name of the new column with nested data
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
        dplyr::mutate(Evidence = 'asr-tax')

    pos <- which(colnames(asr_scores) == 'attr')
    colnames(asr_scores)[pos] <- attr_val_col

    parents <- tibble::as_tibble(get_parents(asr_scores$NCBI_ID))
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


# Helper functions --------------------------------------------------------

## Function to calculate parent score
.calcParentScore <- function(df, attr_col) {
    df <- df[, c(attr_col, 'Score')]
    colnames(df) <- c('attr', 'val')
    df |> dplyr::mutate(
        total = sum(.data$val),
        prop = .data$val / .data$total
    ) |>
        dplyr::count(.data$attr, wt = .data$prop, name = 'Score') |>
        dplyr::mutate(Score = round(.data$Score, 1)) |>
        ## TODO maybe don't apply filter.
        dplyr::filter(.data$Score >= 0.5)
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
