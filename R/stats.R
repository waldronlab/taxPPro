
#' Counts per attribute and rank
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#' @export
#'
counts_per_attribute_and_rank <- function(x) {
    NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
    x |>
        dplyr::filter(
            !is.na(NCBI_ID), NCBI_ID != 'unknown',
            !is.na(Taxon_name),
            Attribute != '', Attribute_value != FALSE,
            !is.na(Rank), Rank != '',
            Rank %in% c('genus', 'species', 'strain')
            ) |>
        dplyr::count(Attribute, Rank)
}


#' Counts per  rank
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#' @export
#'
counts_per_rank <- function(x) {
    NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
    x |>
        dplyr::filter(
            !is.na(NCBI_ID), NCBI_ID != 'unknown',
            !is.na(Taxon_name),
            Attribute != '', Attribute_value != FALSE,
            !is.na(Rank), Rank != '',
            Rank %in% c('genus', 'species', 'strain')
        ) |>
        dplyr::count(Rank)
}

#' Total counts
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#' @export
#'
counts_total <- function(x) {
    NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
    x |>
        dplyr::filter(
            !is.na(NCBI_ID), NCBI_ID != 'unknown',
            !is.na(Taxon_name),
            Attribute != '', Attribute_value != FALSE,
            !is.na(Rank), Rank != '',
            Rank %in% c('genus', 'species', 'strain')
        ) |>
        dplyr::count()
}

