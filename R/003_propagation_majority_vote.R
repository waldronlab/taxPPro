
#' Propagate annotations
#'
#' \code{propagate} propagates annotations.
#'
#' @param df A data frame from bugphyzz.
#' @param max.tax.level A character string. The method that should be used for
#' propagation.
#' @param direction A character string. 'upstream', 'downstream', or 'both'.
#'
#' @return A data frame with propagated annotations.
#'
#' @export
#'
#' @examples
#'
#' library(bugphyzz)
#' aer <- physiologies('aerophilicity')[[1]]
#' aer_prop <- propagate(aer, max.tax.level = 'genus', direction = 'both')
#'
propagate <- function(df, max.tax.level = 'genus', direction = 'both') {

  valid_ranks <- .validRanks()

  if (!max.tax.level %in% valid_ranks) {
    msg <- paste0(
      'Invalid max.tax.level value.',
      ' Value must be one of ', paste0(valid_ranks, collapse = ', ')
    )
    stop(msg, call. = FALSE)
  }

  valid_direction_options <- c('upstream', 'downstream', 'both')

  if (!direction %in% valid_direction_options) {
    msg <- paste0(
      'Invalid direction value.',
      ' Value must be one of ',
      paste0(valid_direction_options, collapse = ', ')
    )
    stop(msg, call. = FALSE)
  }

  df_filtered <- preSteps(df, 'Taxon_name')

  if (direction == 'upstream' || direction == 'downstream') {
    output <- propagateAnnotations(
      df = df_filtered,
      max.tax.level = max.tax.level,
      direction = direction
    )
  } else if (direction == 'both') {
    output <- df_filtered |>
      propagateAnnotations(
        max.tax.level = max.tax.level,
        direction = 'upstream'
      ) |>
      propagateAnnotations(
        max.tax.level = max.tax.level,
        direction = 'downstream'
      )
  }

  return(output)
}

#' Propagate annotations
#'
#' \code{propagateAnnotations}
#'
#' @param df A data frame.
#' @param max.tax.level A character string. Maximum taxonomic rank/level.
#' @param direction A character string. 'upstream' or 'downstream'.
#'
#' @return A data frame with asr/inh evidence.
#' @export
#'
propagateAnnotations <- function(df, max.tax.level, direction) {

  valid_ranks <- .validRanks()

  if (!max.tax.level %in% valid_ranks) {
    msg <- paste0(
      'Invalid max.tax.level value.',
      ' Value must be one of ', paste0(valid_ranks, collapse = ', ')
    )
    stop(msg, call. = FALSE)
  }

  valid_direction_options <- c('upstream', 'downstream', 'both')
  if (!direction %in% valid_direction_options) {
    msg <- paste0(
      'Invalid direction value.',
      ' Value must be one of ',
      paste0(valid_direction_options, collapse = ', ')
    )
    stop(msg, call. = FALSE)
  }

  split_by_rank <- split(df, factor(df$Rank))

  if (direction == 'upstream') {
    valid_ranks <- valid_ranks[1:which(valid_ranks == max.tax.level)]
  } else if (direction == 'downstream') {
    pos <- which(valid_ranks == max.tax.level)
    valid_ranks <- valid_ranks[pos:1]
  }

  for (i in seq_along(valid_ranks)) {

    current_rank <- valid_ranks[i]

    if (current_rank %in% names(split_by_rank)) {
      message('Current rank: ', current_rank)
      next_pos <- i + 1
      if (next_pos > length(valid_ranks))
        break

      next_rank <- valid_ranks[next_pos]
      message('Getting next rank: ', next_rank)

      new_scores <- .getScores(split_by_rank[[current_rank]], direction)

      if (is.null(new_scores))
        break

      new_scores <- new_scores |>
        dplyr::filter(.data$Rank == next_rank)

      if (next_rank %in% names(split_by_rank)) {

        split_by_rank[[next_rank]] <- .replaceTaxa(
          split_by_rank[[next_rank]], new_scores
        )

      } else {
        split_by_rank[[next_rank]] <- new_scores
      }

    } else {
      next
    }
  }

  split_by_rank |>
    dplyr::bind_rows() |>
    dplyr::distinct() |>
    as.data.frame()
}

# Helper functions --------------------------------------------------------

#' Get parents
#'
#' \code{getParents} get's the parent taxon of a vector of valid NCBI IDs.
#'
#' @param x A vector of valid NCBI taxids
#'
#' @return A data frame with Parent_NCBI_ID, Parent_name, and Parent_rank.
#' @export
#'
getParents <- function(x) {
  classification <- taxizedb::classification(x)
  classification <- classification[!is.na(classification)]
  classification <- purrr::map(classification, dplyr::distinct)
  parents_list <- purrr::map(classification, ~ {
    colnames(.x) <- c('Parent_name', 'Parent_rank', 'Parent_NCBI_ID')
    valid_ranks <- .validRanks()
    .x <- .x[.x$Parent_rank %in% valid_ranks, ]
    n_rows <- nrow(.x)
    .x <- .x[-n_rows, ]
    utils::tail(.x, 1)
  })

  lgl_vct <- purrr::map_int(parents_list, ~ nrow(.x) > 0)

  classification <- classification[lgl_vct]

  if (!length(classification))
    return(NULL)

  parents_list <- parents_list[lgl_vct]
  parents_df <- parents_list |>
    dplyr::bind_rows()

  if (!nrow(parents_df))
    return(NULL)
  parents_df$NCBI_ID <- names(classification)
  parents_df[,c('NCBI_ID', 'Parent_NCBI_ID', 'Parent_name', 'Parent_rank')]
}

#' Get children
#'
#' \code{getChildren} gets the immediate descendants of a taxon.
#'
#' @param x A vector of valid NCBI taxids.
#'
#' @return A data frame with all children taxa.
#'
#' @export
#'
getChildren <- function(x) {
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
      data2 = purrr::map(# 'data' is the name of the new column (nested)
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
    dplyr::mutate(Evidence = 'asr')

  pos <- which(colnames(asr_scores) == 'attr')
  colnames(asr_scores)[pos] <- attr_val_col

  ## These are the new parents of the original parents
  parents <- getParents(asr_scores$NCBI_ID)
  if (is.null(parents))
    return(NULL)
  parents <- tibble::as_tibble(parents)
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

#' Get scores for children taxa
#'
#' \code{getChildrenScores} works
#'
#' @param df A data frame imported from bugphyzz and output of the `preSteps`
#' function
#'
#' @return A data frame with scores for children taxa.
#' @export
#'
getChildrenScores <- function(df) {

  ## Current taxa names must become the new parents
  parent_taxa <- df[,!startsWith(colnames(df), 'Parent_')]

  pos <- which(colnames(parent_taxa) == 'NCBI_ID')
  colnames(parent_taxa)[pos] <- 'Parent_NCBI_ID'

  pos <- which(colnames(parent_taxa) == 'Taxon_name')
  colnames(parent_taxa)[pos] <- 'Parent_name'

  pos <- which(colnames(parent_taxa) == 'Rank')
  colnames(parent_taxa)[pos] <- 'Parent_rank'

  pos <- which(colnames(parent_taxa) == 'Evidence')
  colnames(parent_taxa)[pos] <- 'Parent_evidence'

  children_taxa <- getChildren(df$NCBI_ID)

  children_taxa |>
    dplyr::left_join(parent_taxa, by = 'Parent_NCBI_ID') |>
    dplyr::mutate(Evidence = 'inh')

}

# Helper functions --------------------------------------------------------

## Helper function to get scores
.getScores <- function(df, direction) {
  if (direction == 'upstream') {
    return(getParentScores(df))
  } else if (direction == 'downstream') {
    return(getChildrenScores(df))
  }
}

## Replace parents
.replaceTaxa <- function(df, taxa_scores) {
  new_taxa_lgl <- !taxa_scores$Taxon_name %in% df$Taxon_name
  new_taxa <- taxa_scores[new_taxa_lgl,]

  if (!length(new_taxa)) {
    return(df)
  } else {
    return(dplyr::bind_rows(df, new_taxa))
  }
}

## Replace children
.replaceChildren <- NULL

## Function to calculate parent score
.calcParentScore <- function(df, attr_col) {

  pos <- which(colnames(df) == attr_col)
  colnames(df)[pos] <- 'attr'

  pos <- which(colnames(df) == 'Score')
  colnames(df)[pos] <- 'val'

  df_cat <- df |>
    dplyr::group_by(.data[['attr']]) |>
    dplyr::summarise(
      Attribute_source = paste0(.data$Attribute_source, collapse = '|'),
      Confidence_in_curation = paste0(
        .data$Confidence_in_curation, collapse = '|'
      ),
      Original_NCBI_ID = paste0(.data$NCBI_ID, collapse = '|'),
      Original_Taxon_name = paste0(.data$Taxon_name, collapse = '|'),
      Original_Score = paste0(.data$val, collapse = '|'),
      Original_Evidence = paste0(.data$Evidence, collapse = '|')
    ) |>
    dplyr::ungroup()

  new_scores <- df |> dplyr::mutate(
    n = dplyr::n(),
    total = sum(.data$val),
    prop = ifelse(n == 1, .data$val, .data$val / .data$total)
  ) |>
    dplyr::count(.data$attr, wt = .data$prop, name = 'Score') |>
    dplyr::mutate(Score = round(.data$Score, 1))
  # dplyr::filter(.data$Score >= 0.5) # TODO maybe do/don't apply filter
  dplyr::left_join(new_scores, df_cat, by = 'attr')
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
