## These file contains some of the functions used before propagation
## Other functions used here are contained in the 002.duplicates.R file. For
## example, resolveConflicts and resolveDuplicates.

#' Preparare data for propagation part 1
#'
#' \code{prepareData} performs the previous steps before propagating
#' annotations.
#' These steps include: 1) filter data, 2) resolve conflicts, 3) resolve
#' agreements, 4) remove duplicates, 5) convert frequency to scores,
#' 6) remove duplicate lines, 7) ensure Taxid and Parent_tax_id are are
#' character vectors.
#'
#' @param df A data.frame.
#' @param tax.id.type A character string. Valid values are Taxon_name and
#' NCBI_ID. Default is NCBI_ID.
#' @param remove_false If TRUE, Attribute values with FALSE are removed.
#' Default is FALSE.
#' NCBI_ID.
#'
#' @return A data.frame.
#' @export
#'
prepareData <- function(df, tax.id.type = 'NCBI_ID', remove_false = FALSE) {
    df <- df |>
        filterData(
            tax.id.type = tax.id.type, remove_false = remove_false
        ) |>
        freq2Scores() |>
        resolveAgreements() |>
        resolveConflicts()
    if (!nrow(df)) {
        return(NULL)
    } else {
        return(df)
    }
}

#' Filter data for propagation
#'
#' \code{filterData} removes rows that don't have enough information to
#' perform propagation with ASR and Inheritance.
#'
#' @param df A dataset from bugphyzz.
#' @param df_name A character string. The name of the dataset.
#' @param tax.id.type A character string. Valid options: NCBI_ID, Taxon_name.
#' This may be transformed to include other taxonomy IDs.
#' @param remove_false If TRUE, row with FALSE attribute values are removed.
#' Default is FALSE.
#'
#' @return A filtered version of the input data frame.
#' If no rows are kept, the output is NULL with a warning.
#'
#' @export
#'
filterData <- function(df, df_name = NULL, tax.id.type, remove_false = TRUE) {

    columns_for_propagation <- c(
        'Taxon_name', 'NCBI_ID', 'Rank',
        'Attribute', 'Attribute_value', 'Attribute_source',
        'Attribute_value_min', 'Attribute_value_max',
        'Evidence', 'Frequency', 'Confidence_in_curation',
        'Parent_NCBI_ID', 'Parent_name', 'Parent_rank'
    )

    if ('Attribute_value' %in% colnames(df)) {
        columns_for_propagation <- columns_for_propagation[!columns_for_propagation %in% c('Attribute_value_min', 'Attribute_value_max')]
    } else {
        columns_for_propagation <- columns_for_propagation[columns_for_propagation != 'Attribute_value']
    }

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

    if (tax.id.type == 'NCBI_ID') {
        df <- df[!is.na(df$NCBI_ID) | df$NCBI_ID == 'unknown',]
    } else if (tax.id.type == 'Taxon_name') {
        df <- df[!is.na(df$Taxon_name) | df$Taxon_name == 'unknown',]
    } else {
        stop(
            'At the moment, only NCBI_ID or Taxon_name are valid values for',
            'the tax.id.type argument.'
        )
    }

    if (remove_false && 'Attribute_value' %in% colnames(df)) {
       df <- df[which(df$Attribute_value != FALSE),]
    }

    df <- df |>
        dplyr::filter(
            ## Id-related
            !is.na(Rank),
            Rank %in% .validRanks(),

            ## Attribute-related
            !is.na(Attribute) | Attribute != '',
            # Attribute_value != FALSE, ## this was removed above
            !is.na(Attribute_source),

            ## Evidence-related
            !is.na(Evidence),
            !is.na(Frequency),
            !is.na(Confidence_in_curation),

            ## Parent-related
            !is.na(Parent_NCBI_ID),
            !is.na(Parent_name),
            !is.na(Parent_rank),
            Parent_rank %in% .validParentRanks()
        ) |>
        dplyr::distinct()

    if (!nrow(df)) {
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
        return(invisible(NULL))
    }
    return(df)
}

#' Convert frequency values to numeric scores
#'
#' \code{freq2Scores} converts the keywords in the `Frequency`
#' column of a bugphyzz dataset into numeric scores, which are added in a
#' additional column named `Score`.
#'
#' @param df  A data frame imported with `bugphyzz::physiologies`.
#'
#' @return A data frame. The same data frame with the additional `Score` column.
#'
#' @export
#'
freq2Scores <- function(df) {
    attr_type <- unique(df$Attribute_type)
    if (attr_type %in% c('numeric', 'range')) {
        df$Frequency <- ifelse(
            df$Frequency == 'unknown', 'always', df$Frequency
        )
    }
    output <- df |>
        dplyr::mutate(
            Frequency = tolower(.data$Frequency),
            Score = dplyr::case_when(
                Frequency == 'always' ~ 1,
                Frequency == 'usually' ~ 0.8,
                Frequency == 'sometimes' ~ 0.5,
                Frequency == 'rarely' ~ 0.2,
                Frequency == 'never' ~ 0,
                Frequency == 'unknown' ~ NA_real_
                # Frequency == 'unknown' ~ 0
            )
        ) |>
        dplyr::filter(!is.na(.data$Score)) |>
        dplyr::distinct()
    return(output)
}

#' Prepare data for propagation part 2
#'
#' \code{prepareData2} prepares data for propagation. Must be used after
#' \code{prepareData}.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
prepareData2 <- function(df) {
    dict <- c(genus = 'g__', species = 's__', strain = 't__')
    df$NCBI_ID <- paste0(dict[df$Rank], df$NCBI_ID)
    attr_type <- unique(df$Attribute_type)
    if (attr_type == 'logical') {
        select_cols <- c(
            'NCBI_ID', 'Attribute', 'Evidence', 'Attribute_source', 'Score',
            'Attribute_group', 'Attribute_type'
        )
        output <- df |>
            dplyr::select(dplyr::all_of(select_cols)) |>
            dplyr::distinct() |>
            dplyr::filter(grepl('^[gst]__', .data$NCBI_ID)) |>
            dplyr::distinct()
    } else if (attr_type == 'numeric') {
        df <- numericToRange(df)
        df$Attribute_type <- 'range'
        select_cols <- c(
            'NCBI_ID', 'Attribute_value_min', 'Attribute_value_max',
            'Evidence', 'Attribute_source', 'Score', 'Attribute_group',
            'Attribute_type'
        )
        output <- df |>
            dplyr::select(dplyr::all_of(select_cols)) |>
            dplyr::distinct() |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate_at(
                .vars = c('Attribute_value_min', 'Attribute_value_max', 'Score'),
                .funs = ~ mean(.x)
            ) |>
            dplyr::mutate_at(
                .vars = c('Evidence', 'Attribute_source'),
                .funs = ~ paste0(unique(.x), collapse = '|')
            ) |>
            dplyr::filter(grepl('^[gst]__', .data$NCBI_ID)) |>
            dplyr::distinct()
        # select_cols <- c(
        #     'NCBI_ID', 'Attribute_value', 'Evidence', 'Attribute_source',
        #     'Score', 'Attribute_group', 'Attribute_type'
        # )
        # output <- df |>
        #     dplyr::select(dplyr::all_of(select_cols)) |>
        #     dplyr::distinct() |>
        #     dplyr::group_by(.data$NCBI_ID) |>
        #     dplyr::mutate_at(
        #         .vars = c('Attribute_value', 'Score'),
        #         .funs = ~ mean(.x)
        #     ) |>
        #     dplyr::mutate_at(
        #         .vars = c('Evidence', 'Attribute_source'),
        #         .funs = ~ paste0(unique(.x), collapse = '|')
        #     ) |>
        #     dplyr::filter(grepl('^[gst]__', .data$NCBI_ID)) |>
        #     dplyr::distinct()
    } else if (attr_type == 'range') {
        select_cols <- c(
            'NCBI_ID', 'Attribute_value_min', 'Attribute_value_max',
            'Evidence', 'Attribute_source', 'Score', 'Attribute_group',
            'Attribute_type'
        )
        output <- df |>
            dplyr::select(dplyr::all_of(select_cols)) |>
            dplyr::distinct() |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate_at(
                .vars = c('Attribute_value_min', 'Attribute_value_max', 'Score'),
                .funs = ~ mean(.x)
            ) |>
            dplyr::mutate_at(
                .vars = c('Evidence', 'Attribute_source'),
                .funs = ~ paste0(unique(.x), collapse = '|')
            ) |>
            dplyr::filter(grepl('^[gst]__', .data$NCBI_ID)) |>
            dplyr::distinct()
    }
    return(output)
}

#' Convert numeric atttributes to range
#'
#' \code{numericToRange} convers numeric attributes to range attributes.
#' This function is meant to be used in \code{prepareData2}.
#'
#' @param df A data frame.
#'
#' @return A data.frame
#' @export
#'
numericToRange <- function(df) {
    df <- df |>
        dplyr::filter(!is.na(NCBI_ID), NCBI_ID != 'unknown') |>
        dplyr::group_by(NCBI_ID) |>
        dplyr::mutate(
            Attribute_value_min = min(.data$Attribute_value),
            Attribute_value_max = max(.data$Attribute_value)
        ) |>
        dplyr::select(-.data$Attribute_value) |>
        dplyr::distinct()
    return(df)
}
