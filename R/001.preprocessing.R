## These file contains some of the functions used before propagation
## Other functions are contained in the 002.duplicates.R file.

#' \code{preSteps} perform the previous steps before propagating annotations.
#' These steps include: 1) filter data, 2) resolve conflicts, 3) resolve
#' agreements, 4) remove duplicates, 5) convert frequency to scores, 6)
#' remove duplicate lines, 7) ensure Taxid and Parent_tax_id are are character
#' vectors.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
preSteps <- function(df, tax.id.type) {
    NCBI_ID <- Parent_NCBI_ID <- NULL
    df |>
        filterData(tax.id.type = tax.id.type) |>
        resolve_conflicts() |>
        resolve_agreements() |>
        remove_taxa_duplicates() |>
        freq2Scores() |>
        dplyr::distinct() |>
        dplyr::mutate(
            Parent_NCBI_ID = as.character(Parent_NCBI_ID),
            NCBI_ID = as.character(NCBI_ID)
        )
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
#'
#' @return A filtered version of the input data frame.
#' If no rows are kept, the output is NULL with a warning.
#'
#' @export
#'
filterData <- function(df, df_name = NULL, tax.id.type) {

    ## Columns required for propagation
    columns_for_propagation <- c(
        'Taxon_name', 'NCBI_ID', 'Rank', # id-related
        'Attribute', 'Attribute_value', 'Attribute_source', # attribute-related
        'Evidence', 'Frequency', 'Confidence_in_curation', # evidence-related
        'Parent_NCBI_ID', 'Parent_name', 'Parent_rank' # parent-related
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


    ## For now, only NCBI_ID and Taxon_name are valid. Other columns might
    ## be added later

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

    df <- df |>
        dplyr::filter(
            ## Id-related
            !is.na(Rank),
            Rank %in% .validRanks(),

            ## Attribute-related
            !is.na(Attribute) | Attribute != '',
            Attribute_value != FALSE,
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
        return(NULL)
    }

    return(df)
}

#' Convert frequency values to numeric scores
#'
#' \code{freq2Scores} converts the keywords in the `Frequency`
#' column of a bugphyzz dataset into numeric scores, which are added in a
#' additional column named `Score`.
#'
#' @param x  A data frame imported with `bugphyzz::physiologies`.
#'
#' @return A data frame. The same data frame with the additional `Score` column.
#'
#' @export
#'
freq2Scores <- function(x) {
    x |>
        dplyr::mutate(
            Frequency = tolower(.data$Frequency),
            Score = dplyr::case_when(
                Frequency == 'always' ~ 1,
                Frequency == 'usually' ~ 0.8,
                Frequency == 'sometimes' ~ 0.5,
                Frequency == 'rarely' ~ 0.2,
                Frequency == 'never' ~ 0,
                Frequency == 'unknown' ~ NA_real_
            )
        ) |>
        dplyr::filter(!is.na(.data$Score))
}
