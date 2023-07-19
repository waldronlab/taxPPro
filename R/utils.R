#' Convert frequency values to numeric scores
#'
#' \code{freq2Scores} converts the keywords in the `Frequency`
#' column of a bugphyzz dataset into numeric scores, which are added in a
#' additional column named `Score`.
#'
#' @param x A  character vector.
#'
#' @return A numeric vector.
#'
#' @export
#'
freq2Scores <- function(x) {
  # attr_type <- unique(df$Attribute_type)
  # if (attr_type %in% c('numeric', 'range')) {
  #   df$Frequency <- ifelse(
  #     df$Frequency == 'unknown', 'always', df$Frequency
  #   )
  # }
  x <- tolower(x)
  dplyr::case_when(
    x == 'always' ~ 1,
    x == 'usually' ~ 0.9,
    x == 'sometimes' ~ 0.5,
    x == 'unknown' ~ 0.1
  )
}

#' Scores to frequency
#'
#' \code{scores2Freq} performs the opposite of \code{freq2scores}.
#'
#' @param x A numeric vector.
#'
#' @return A character vector
#' @export
#'
scores2Freq <- function(x) {
    dplyr::case_when(
        x == 1 ~ 'always',
        x >= 0.9 & x < 1 ~ 'usually',
        ## often is missing here, which is 70%. That's why I
        ## previously lowered the threshold to 0.8.
        ## Now sometimes range is between 0.5 and 0.8999999999999
        ## Which means that we're not really using grammarist ranges (our source)
        x >= 0.5 & x < 0.9 ~ 'sometimes',
        x > 0 & x < 0.5 ~ 'rarely',
        x == 0  ~ 'never'
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
