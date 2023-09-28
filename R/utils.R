#' Frequency to scores
#'
#' \code{freq2Scores} converts the keywords in the "Frequency"
#' column of a bugphyzz dataset into numeric scores, which are interpreted
#' as probabilities.
#'
#' @param x A character vector.
#'
#' @return A numeric vector.
#'
#' @export
#'
freq2Scores <- function(x) {
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
#' \code{scores2Freq} converts numeric scores, stored in the "Score" column in
#' a bugphyzz dataset, into frequency keywords.
#'
#' @param x A numeric vector.
#'
#' @return A character vector.
#'
#' @export
#'
scores2Freq <- function(x) {
    dplyr::case_when(
        x == 1 ~ 'always',
        x >= 0.9 & x < 1 ~ 'usually',
        x >= 0.5 & x < 0.9 ~ 'sometimes',
        x > 0 & x < 0.5 ~ 'rarely',
        x == 0  ~ 'never'
    )
}

#' Confidence in curation to factor
#'
#' \code{conf2Fct} converts a vector from the Confidence_in_curation
#' column in a bugphyzz dataset to a factor.
#'
#' @param x A character vector.
#'
#' @return A factor.
#'
#' @export
#'
conf2Fct <- function(x) {
    factor(x, levels = c('low', 'medium', 'high'), ordered = TRUE)
}

#' Add rank prefix to NCBI_ID
#'
#' \code{addRankPrefix} adds a prefix (`'[kpcofgst]__'`) to an NCBI ID (taxid)
#'
#' @param id A charecter vector of NCBI IDs (taxids).
#' @param rank A character vector of ranks
#'
#' @return A character vector of IDs with prefix.
#' @export
#'
addRankPrefix <- function(id, rank) {
    highest_ranks <- c('kingdom', 'superkingdom', 'domain')
    case_when(
        rank %in% highest_ranks ~ paste0('k__', id),
        rank == 'phylum' ~ paste0('p__', id),
        rank == 'class' ~ paste0('c__', id),
        rank == 'order' ~ paste0('o__', id),
        rank == 'family' ~ paste0('f__', id),
        rank == 'genus' ~ paste0('g__', id),
        rank == 'species' ~ paste0('s__', id),
        rank == 'strain' ~ paste0('t__', id),
        TRUE ~ NA
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


myFun <- function(x, adjF = 0.1) {
    n <- length(x)
    target_scores <- rep(1 / n, n)
    score_diff <- x - target_scores
    res <- x - adjF * score_diff
    return(res)
}

#' Clean node
#'
#' {cleanNode} deletes the `attribute_tbl` slot in a data.tree node.
#'
#' @param node A node in a data.tree object.
#'
#' @return NULL (assigned to the node)
#' @export
#'
cleanNode <- function(node) {
    node$attribute_tbl <- NULL
}
