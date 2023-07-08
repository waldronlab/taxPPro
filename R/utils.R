
#' Scores to frequency
#'
#' \code{scores2Freq} performs the opposite of \code{freq2scores}.
#'
#' @param Score A numeric vector.
#'
#' @return A character vector
#' @export
#'
scores2Freq <- function(Score) {
    dplyr::case_when(
        Score == 1 ~ 'always',
        Score >= 0.7 & Score < 1 ~ 'usually',
        Score >= 0.4 & Score < 0.7 ~ 'sometimes',
        # Score >= 0.1 & Score < 0.4 ~ 'rarely',
        # Score < 0.1 ~ 'never'
        Score >= 0.1 & Score < 0.4 ~ 'unknown',
        Score < 0.1 ~ 'unknown'
    )
}

## Function with valid ranks for taxa
.validRanks <- function() {
    c('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum',
      'superkingdom')
}

## Function with valid ranks for parents
.validParentRanks <- function() c('species', 'genus', 'family')

#' Taxize classification to data frame
#'
#' \code{classif2Table} converts from classification format (`taxize` packatge)
#' to a data frame.
#'
#' @param x Output from the `taxize::classification` function.
#' @param ranks Ranks to be selected. Default is given by `.validRanks()`.
#'
#' @return A data frame.
#' @export
#'
classif2Table <- function(x, ranks) {

  if (missing(ranks)) {
    valid_ranks <- .validRanks()
  } else {
    valid_ranks <- ranks
  }

  df_filtered <- x |>
    dplyr::select(rank, id) |>
    dplyr::filter(rank %in% valid_ranks)

  new_df <- data.frame(x = df_filtered$id) |>
    t() |>
    as.data.frame(row.names = 1L)
  colnames(new_df) <- df_filtered$rank
  new_df
}

