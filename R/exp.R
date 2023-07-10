#' My function
#'
#' \code{myFun} will be used instead of prepareDataForPropagation.
#'
#' @param df A data.frame.
#'
#' @return A named list.
#' @export
#'
myFun <- function(df) {
    df <- df[which(!is.na(df$Rank)), ]
    df <- df[which(!is.na(df$Evidence)), ]
    df <- df[which(!is.na(df$Frequency)), ]
    df <- df[which(!is.na(df$Confidence_in_curation)), ]
    df$Score <- freq2Scores(df$Frequency)
    df$NCBI_ID[which(is.na(df$NCBI_ID))] <- 'unknown'

    ## Original annotations with TAXID
    original <- df[which(df$NCBI_ID != 'unknown'),]
    if (nrow(original) > 0) {
        original <- original |>
            dplyr::select(
                -.data$Parent_NCBI_ID, -.data$Parent_name, -.data$Parent_rank
            ) |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate(
                Taxon_name = paste(unique(.data$Taxon_name), collapse = ';')
            ) |>
            dplyr::ungroup() |>
            dplyr::distinct() |>
            purrr::discard(~ all(is.na(.x))) |>
            dplyr::mutate(
                NCBI_ID = sub('^(\\w)\\w+(__.*)$', '\\1\\2', paste0(Rank, '__', NCBI_ID))
            ) |>
            dplyr::select(-.data$Rank) |>
            removeAccessionAndGenomeID() |>
            dplyr::distinct() |>
            resolveAgreements() |>
            resolveConflicts() |>
            dplyr::distinct() |>
            as.data.frame()
    }

    ## First step of ASR for entries without taxids (they can't be mapped to the data.tree structure)
    early_asr <- df[which(df$NCBI_ID == 'unknown'),]
    if (nrow(early_asr) > 0) {
        early_asr <- calcParentScores(early_asr) |>
            removeAccessionAndGenomeID() |>
            dplyr::distinct() |>
            purrr::discard(~ all(is.na(.x))) |>
            dplyr::mutate(
                NCBI_ID = sub('^(\\w)\\w+(__.*)$', '\\1\\2', paste0(Rank, '__', NCBI_ID))
            ) |>
            dplyr::select(-.data$Rank) |>
            dplyr::distinct() |>
            as.data.frame()
    }

    output <- list(
        original = original, early_asr = early_asr
    )
    return(output)
}

#' Remove Accession_ID and Genome_ID columns
#'
#' \code{removeAcccessionAndGenomeID} removes Accession_ID and Genome_ID from
#' a bugphyzz dataset. The reason is that right now these columns can be
#' incomplete or inconsistent in some datasets, or just missing in some others.
#' I think a solution would be to implement a relational database in which we
#' have a data object (data.frame?) with all of the taxids.
#'
#' @param df A data.frame imported from bugphyzz.
#'
#' @return A data.frame.
#' @export
#'
removeAccessionAndGenomeID <- function(df) {
    if ('Accession_ID' %in% colnames(df)) {
        df <- dplyr::select(df, -.data[['Accession_ID']])
    }
    if ('Genome_ID' %in% colnames(df)) {
        df <- dplyr::select(df, -.data[['Genome_ID']])
    }
    dplyr::distinct(df)
}


#' Merge original and early ASR
#'
#' \code{mergeOriginalAndEarlyASR} merges original bugphyzz annotations
#' and early ASR entries (those taxids that couldn't be mapped to the
#' data.tree structures, but their parents could).
#'
#' @param l A list. Output of \code{myFun}.
#'
#' @return A data.frame.
#' @export
#'
mergeOriginalAndEarlyASR <- function(l) {
    df <- l |>
        dplyr::bind_rows() |>
        {\(y) split(y, factor(y$NCBI_ID))}() |>
        purrr::map(~ {
            if (!all(.x$Evidence == 'asr')) {
                output <- filter(.x, Evidence != 'asr')
            } else {
                output <- .x
            }
            return(output)
        }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(Evidence = forcats::fct_relevel(Evidence, 'asr')) |>
        dplyr::arrange(Evidence)
}
