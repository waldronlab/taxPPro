# filterDirectParents <- function(df) {
#     vct <- dplyr::case_when(
#         df$Rank == 'strain' & df$Parent_rank == 'species' ~ TRUE,
#         df$Rank == 'species' & df$Parent_rank == 'genus' ~ TRUE,
#         df$Rank == 'genus' & df$Parent_rank == 'family' ~ TRUE,
#         df$Rank == 'family' & df$Parent_rank == 'order' ~ TRUE,
#         df$Rank == 'order' & df$Parent_rank == 'class' ~ TRUE,
#         df$Rank == 'class' & df$Parent_rank == 'phylum' ~ TRUE,
#         TRUE ~ FALSE
#     )
#     return(df[vct,])
# }

#' Prepare Data for Propagation
#'
#' \code{prepareDataForPropagation} prepares data to be used in the
#' propagation method.
#'
#' @param df A data.frame.
#'
#' @return A data.frame
#' @export
#'
prepareDatForPropagation <- function(df) {
    df$NCBI_ID[which(is.na(df$NCBI_ID))] <- 'unknown'
    df$Parent_NCBI_ID[which(is.na(df$Parent_NCBI_ID))] <- 'unknown'
    df <- df[df$Parent_NCBI_ID != 'unknown',]
    df <- df[which(!is.na(df$Rank)), ]
    df <- df[which(!is.na(df$Evidence)), ]
    df <- df[which(!is.na(df$Frequency)), ]
    df <- df[which(!is.na(df$Confidence_in_curation)), ]
    df <- df |>
        dplyr::group_by(.data$Parent_NCBI_ID) |>
        dplyr::mutate(
            Parent_name = paste(unique(.data$Parent_name), collapse = ';')
        ) |>
        dplyr::ungroup() |>
        dplyr::distinct()
    df <- unique(df)
    df <- freq2Scores(df)
    df_yesid <- df[which(df$NCBI_ID != 'unknown'),]
    if (nrow(df_yesid) > 0) {
        df_yesid <- df_yesid |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate(
                Taxon_name = paste(unique(.data$Taxon_name), collapse = ';')
            ) |>
            dplyr::ungroup() |>
            dplyr::distinct()
        # df_yesid <- df_yesid[grep(';', df_yesid$Taxon_name),]
    }
    df_noid <- df[which(df$NCBI_ID == 'unknown'),]
    if (nrow(df_noid) > 0) {
        df_noid_asr <- calcParentScores(df_noid)
        df_new <- dplyr::bind_rows(df_yesid, df_noid_asr)
    } else {
        df_new <- df_yesid
    }
    dict <- c(genus = 'g__', species = 's__', strain = 't__')
    df_new$NCBI_ID <- paste0(dict[df_new$Rank], df_new$NCBI_ID)
    df_new <- df_new[which(!startsWith(colnames(df_new), 'Parent'))]
    df_new <- unique(df_new)
    output <- df_new |>
        resolveAgreements() |>
        resolveConflicts() |>
        dplyr::distinct() |>
        as.data.frame()
    return(output)

    # attr_type <- unique(x$Attribute_type)
    # cols <- c(
    #     'NCBI_ID', 'Taxon_name', 'Rank',
    #     'Attribute', 'Attribute_source',
    #     'Evidence', 'Frequency',
    #     'Attribute_type', 'Attribute_group',
    #     'Confidence_in_curation', 'Score'
    # )
    # if (attr_type == 'logical') {
    #     cols <- c(cols, 'Attribute_value')
    # } else if (attr_type == 'range') {
    #     cols <- c(cols, c('Attribute_value_min', 'Attribute_value_max'))
    # }
    #
    # x_new <- unique(x_new[,cols])
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

#' Convert frequency values to numeric scores
#'
#' \code{freq2Scores} converts the keywords in the `Frequency`
#' column of a bugphyzz dataset into numeric scores, which are added in a
#' additional column named `Score`.
#'
#' @param df  A data frame imported with `bugphyzz::physiologies`.
#'
#' @return A data frame. The same data frame with the additional `Score` column.
#'
#' @export
#'
freq2Scores <- function(df) {
    attr_type <- unique(df$Attribute_type)
    if (attr_type %in% c('numeric', 'range')) {
        df$Frequency <- ifelse(
            df$Frequency == 'unknown', 'always', df$Frequency
        )
    }
    output <- df |>
        dplyr::mutate(
            Frequency = tolower(.data$Frequency),
            Score = dplyr::case_when(
                Frequency == 'always' ~ 1,
                Frequency == 'usually' ~ 0.8,
                Frequency == 'sometimes' ~ 0.5,
                # Frequency == 'rarely' ~ 0.2,
                # Frequency == 'never' ~ 0,
                # Frequency == 'unknown' ~ NA_real_
                Frequency == 'unknown' ~ 0.1
            )
        ) |>
        dplyr::filter(!is.na(.data$Score)) |>
        dplyr::distinct()
    return(output)
}
