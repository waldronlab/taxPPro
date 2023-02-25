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
