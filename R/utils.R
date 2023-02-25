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
