#' Print data.tree attributes
#'
#' \code{printDataTreeAttributes} prints all of the attributes in a data.tree.
#'
#' @param data_tree A data.tree object with attributes.
#' @param limit The number of nodes to be displayed. Default = 100.
#'
#' @return A data.frame or data on the console. Can use \code{View}.
#' @export
#'
printDataTreeAttributes <- function(data_tree, limit = 100) {
  attrs <- as.list(data_tree$attributesAll)
  lim = list(limit = limit)
  args <- c(list(data_tree), attrs)
  args <- c(args, lim)
  do.call('print', args = args)
}

#' toDataFrame
#'
#' \code{toDataFrame}
#'
#' @param data_tree A data.tree.
#' @param ncbi_tax NCBI taxonomy. Output of \code{get_ncbi_taxonomy}.
#'
#' @return A data.frame.
#' @export
#'
toDataFrame <- function(data_tree, ncbi_tax) {
  args <- as.list(data_tree$attributesAll)
  args <- c(list(x = data_tree, row.names = NULL, optional = FALSE), args)
  df <- do.call('as.data.frame', args)
  df$levelName <- stringr::str_squish(sub('.*-', '', df$levelName))
  df <- df[df$levelName != 'ArcBac',]
  df <- tidyr::separate(
    df, col = 'levelName', into = c('Rank', 'NCBI_ID'), sep = '__'
  )
  dict <- c(
    d = 'domain', p = 'phylum', c = 'class', o = 'order', f = 'family',
    g = 'genus', s = 'species', t = 'strain'
  )
  df$Rank <- dict[df$Rank]
  # ncbi_tax <- get_ncbi_taxonomy()
  if ('Rank' %in% colnames(ncbi_tax)) {
    ncbi_tax$Rank <- NULL
  }
  pos <- which(colnames(ncbi_tax) == 'kingdom')
  names(ncbi_tax)[pos] <- 'domain'
  output <- dplyr::left_join(df, ncbi_tax, by = 'NCBI_ID')
  return(output)
}

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








