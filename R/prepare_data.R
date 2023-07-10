#' Prepare data for propagation
#'
#' \code{prepareDataForPropagation} prepares data for propagation.
#'
#' @param df A data.frame.
#'
#' @return A named list.
#' @export
#'
prepareDataForPropagation <- function(df) {
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
        early_asr <- early_asr |>
            dplyr::group_by(.data$Parent_NCBI_ID) |>
            dplyr::mutate(
                Taxon_name = paste(unique(.data$Taxon_name), collapse = ';'),
                Parent_name = paste(unique(.data$Parent_name), collapse = ';')
            ) |>
            dplyr::ungroup() |>
            calcParentScores() |>
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
        dplyr::arrange(Evidence) |>
        dplyr::mutate(Evidence = as.character(Evidence)) |>
        dplyr::distinct()
}

#' Calculate parent scores
#'
#' \code{calcParentScores} calculate parent scores with ASR.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
calcParentScores <- function(df) {
    attr_type <- unique(df$Attribute_type)
    attr_group <- unique(df$Attribute_group)
    if (attr_type == 'logical') {
        output <- df |>
            dplyr::group_by(
                .data$Parent_name, .data$Parent_NCBI_ID, .data$Parent_rank,
                .data$Attribute_type, .data$Attribute_group
            ) |>
            dplyr::mutate(
                n = dplyr::n(),
                total = sum(.data$Score),
                prop = .data$Score / .data$total
            ) |>
            dplyr::count(.data$Attribute, wt = .data$prop, name = 'Score') |>
            dplyr::mutate(Score = round(.data$Score, 1)) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Evidence = 'asr',
                Attribute_value = TRUE,
                Attribute_source = NA,
                Note = NA,
                Frequency = taxPPro::scores2Freq(Score),
                Attribute_type = attr_type,
                Attribute_group = attr_group
            )
        colnames(output)[which(colnames(output) == 'Parent_name')] <- 'Taxon_name'
        colnames(output)[which(colnames(output) == 'Parent_NCBI_ID')] <- 'NCBI_ID'
        colnames(output)[which(colnames(output) == 'Parent_rank')] <- 'Rank'
    } else if  (attr_type == 'range') {
        output <- df |>
            dplyr::group_by(
                .data$Parent_name, .data$Parent_NCBI_ID, .data$Parent_rank,
                .data$Attribute_type, .data$Attribute_group
            ) |>
            dplyr::mutate(
                min = min(.data$Attribute_value_min),
                max = max(.data$Attribute_value_max)
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Evidence = 'asr',
                Attribute_value = TRUE,
                Attribute_type = attr_type,
                Attribute_group = attr_group
            )
        select_cols <- c(
            'Parent_name', 'Parent_NCBI_ID', 'Parent_rank',
            'Attribute', 'min', 'max', 'Evidence',
            'Frequency', 'Score', 'Attribute_type', 'Attribute_group'
        )
        output <- output[,select_cols]
        colnames(output)[which(colnames(output) == 'Parent_name')] <- 'Taxon_name'
        colnames(output)[which(colnames(output) == 'Parent_NCBI_ID')] <- 'NCBI_ID'
        colnames(output)[which(colnames(output) == 'Parent_rank')] <- 'Rank'
        colnames(output)[which(colnames(output) == 'min')] <- 'Attribute_value_min'
        colnames(output)[which(colnames(output) == 'max')] <- 'Attribute_value_max'
    }
    output <- unique(output)
    return(output)
}

#' Resolve agreements
#'
#' \code{resolveAgreements} resolve all agreeents in a bugphyzz dataset.
#' Agreements happen when the same annotations comes from two or more
#' sources.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
resolveAgreements <- function(df) {
    attr_type <- unique(df$Attribute_type)
    df$Confidence_in_curation <- ifelse(
        is.na(df$Confidence_in_curation), 'unknown', df$Confidence_in_curation
    )
    df$Confidence_in_curation <- factor(
        x = df$Confidence_in_curation,
        levels = c('unknown', 'low', 'medium', 'high'),
        ordered = TRUE
    )
    if (attr_type == 'logical') {
        output <- df |>
            dplyr::group_by(
                .data[['NCBI_ID']],
                .data[['Attribute']],
                .data[['Attribute_value']]
            ) |>
            dplyr::slice_max(
                order_by = .data$Confidence_in_curation,
                with_ties = FALSE
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Confidence_in_curation = as.character(.data$Confidence_in_curation)
            )

    } else if (attr_type == 'range') {
        output <- df |>
            dplyr::group_by(
                .data[['NCBI_ID']],
                .data[['Attribute']],
                .data[['Attribute_value_min']],
                .data[['Attribute_value_max']]
            ) |>
            dplyr::slice_max(
                order_by = .data$Confidence_in_curation,
                with_ties = FALSE
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Confidence_in_curation = as.character(.data$Confidence_in_curation)
            )
    }
    output$Confidence_in_curation <- ifelse(
        test = output$Confidence_in_curation == 'unknown',
        yes = NA,
        no = output$Confidence_in_curation
    )
    return(output)
}

#' Resolve conflicts
#'
#' \code{resolveConflicts} resolves are conflicting annotations in a bugphyzz
#' dataset. Conflicts happen when a single taxon recieves two or more
#' annotations from different sources.
#'
#' @param df A data.frame.
#'
#' @return A data.frame.
#' @export
#'
resolveConflicts <- function(df) {
    attr_type <- unique(df$Attribute_type)
    df$Confidence_in_curation <- ifelse(
        is.na(df$Confidence_in_curation), 'unknown', df$Confidence_in_curation
    )
    df$Confidence_in_curation <- factor(
        x = df$Confidence_in_curation,
        levels = c('unknown', 'low', 'medium', 'high'),
        ordered = TRUE
    )
    if (attr_type == 'logical') {
        output <- df |>
            dplyr::group_by(
                .data$NCBI_ID,
            ) |>
            dplyr::slice_max(
                order_by = .data$Confidence_in_curation,
                with_ties = TRUE
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Confidence_in_curation = as.character(.data$Confidence_in_curation)
            )
    } else if (attr_type == 'range') {
        ## This part is performed in three steps.
        ## First we get the mean of values from the same source
        output <- df |>
            dplyr::group_by(
                .data$NCBI_ID, .data$Attribute_source
            ) |>
            dplyr::mutate(
                Attribute_value_min = mean(Attribute_value_min, na.rm = TRUE),
                Attribute_value_max = mean(Attribute_value_max, na.rm = TRUE)
            ) |>
            dplyr::ungroup() |>
            dplyr::distinct() |>
            ## Second, we get the source with the highest confidence in
            ## curation allowing ties.
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::slice_max(
                order_by = .data$Confidence_in_curation,
                with_ties = TRUE
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(
                Confidence_in_curation = as.character(.data$Confidence_in_curation)
            ) |>
            ## Third, in case there are still duplicates, because of ties,
            ## we get the mean and combine the names of the attribute sources
            ## into a single one separated by ';''
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate(
                Attribute_value_min = mean(Attribute_value_min, na.rm = TRUE),
                Attribute_value_max = mean(Attribute_value_max, na.rm = TRUE),
                Attribute_source = paste(sort(unique(Attribute_source)), collapse = ';')
            ) |>
            dplyr::distinct()
    }
    return(output)
}
