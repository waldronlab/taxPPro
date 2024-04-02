#' Filter data
#'
#' \code{filterData} is the first step in the propagation process.
#' It ensures that only data that can be used for propagation is retained.
#'
#' @param tbl A data.frame.
#'
#' @return A data.frame.
#' @export
#'
filterData <- function(tbl) {
    discrete_types <- c("binary", "multistate-intersection", "mltistate-union")
    attr_type <- unique(tbl$Attribute_type)
    if (attr_type %in% discrete_types){
        output <- filterDataDiscrete(tbl)
    } else if (attr_type == "range") {
        output <- filterDataNumeric(tbl)
    } else {
        output <- NULL
    }
    return(output)
}

#' Filter data with numeric attributes
#'
#' \code{filterDataNumeric} filters data that could be used for propagation
#' of numeric attributes. These are labeled as "range" when being imported
#' with the \code{physiologies function}. It also creates a Score column
#' based on the Frequency column.
#'
#' @param tbl A data.frame imported with \code{bugphyzz::physiologies}
#'
#' @return A data.frame.
#' @export
#'
filterDataNumeric <- function(tbl) {
    select_cols <- c(
        'NCBI_ID', 'Taxon_name', 'Parent_NCBI_ID',
        'Attribute_value_min', 'Attribute_value_max',
        'Attribute_source', 'Confidence_in_curation',
        'Frequency', 'Score', 'Evidence',
        'Attribute_type', 'Attribute_group'
    )
    filtered_tbl <- tbl |>
        ## Attribute_value_min and Attribute_value_max are the same when
        ## only one temperature is annotated. So neither columns should be
        ## NA.
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::filter(
            ## We need NCBI ID data. If not available, we need the NCBI ID of
            ## the immediate parent. If this is not available then the entry
            ## is dropped.
            !((is.na(.data$NCBI_ID) | .data$NCBI_ID == 'unknown') & is.na(.data$Parent_NCBI_ID))
        ) |>
        ## We need a source and a frequency
        dplyr::filter(!is.na(.data$Attribute_source), !is.na(.data$Frequency)) |>
        dplyr::mutate(Score = freq2Scores(.data$Frequency)) |>
        dplyr::select(tidyselect::all_of(select_cols)) |>
        dplyr::distinct()

    n_rows_dropped <- nrow(tbl) - nrow(filtered_tbl)

    if (!nrow(filtered_tbl)) {
        message(
            "Not enough data. Dropping :",
            format(n_rows_dropped, big.mark = ","), " rows."
        )
        return(NULL)
    }
    message(format(n_rows_dropped, big.mark = ','), ' rows were dropped.')
    return(filtered_tbl)
}

#' Filter data with discrete attributes
#'
#' \code{filterDataDiscrete} filters data that could be used for propagation
#' of discrete attributes. It also creates a Scores column based on
#' the Frequency column.
#'
#' @param tbl A data.frame imported with \code{bugphyzz::physiologies}
#'
#' @return A data.frame.
#' @export
#'
filterDataDiscrete <- function(tbl) {
    phys_name <- unique(tbl$Attribute_group)
    select_cols <- c(
        'NCBI_ID', 'Taxon_name', 'Parent_NCBI_ID',
        'Attribute','Attribute_source', 'Confidence_in_curation',
        'Frequency', 'Score', 'Evidence',
        'Attribute_type', 'Attribute_group'
    )
    attributes_fname <- system.file(
        'extdata', 'attributes.tsv', package = 'bugphyzz'
    )

    ## This ensures that only valid attributes are imported
    attributes <- utils::read.table(attributes_fname, sep = '\t', header = TRUE)
    rgx <- paste0('\\b', phys_name, '\\b')
    valid_attributes <- attributes |>
        dplyr::filter(grepl(rgx, .data$attribute_group)) |>
        dplyr::pull(.data$attribute) |>
        unique()

    n_rows <- nrow(tbl)

    tbl <- tbl |>
        dplyr::filter(!is.na(.data$Attribute)) |>
        dplyr::filter(!is.na(.data$Attribute_value))

    allAttrs <- unique(tbl$Attribute)
    invalid_attrs <- allAttrs[which(!allAttrs %in% valid_attributes)]

    if (length(invalid_attrs) > 0) {
        invalid_attrs <- paste0(invalid_attrs, collapse = "---")
        wng_msg <- paste0(
            "Invalid attributes in ", phys_name, ": ", invalid_attrs
        )
        warning(wng_msg, call. = FALSE)
    }

    attr_type <- unique(tbl$Attribute_type)
    if (attr_type == 'multistate-intersection') {
        tbl <- tbl |>
            dplyr::filter(.data$Attribute %in% valid_attributes) |>
            dplyr::filter(.data$Attribute_value == TRUE)
    } else if (attr_type %in% c('binary', 'multistate-union')) {
        tbl <- tbl |>
            dplyr::filter(.data$Attribute %in% valid_attributes) |>
            dplyr::mutate(
                Attribute = paste0(.data$Attribute, '--', .data$Attribute_value)
            )
    }

    phys_data <- tbl |>
        tibble::as_tibble() |>
        dplyr::filter(
            !((is.na(.data$NCBI_ID) | .data$NCBI_ID == 'unknown') & is.na(.data$Parent_NCBI_ID))
        ) |>
        dplyr::filter(!is.na(.data$Attribute_source), !is.na(.data$Frequency)) |>
        dplyr::mutate(Score = freq2Scores(.data$Frequency)) |>
        dplyr::select(tidyselect::all_of(select_cols)) |>
        dplyr::distinct()
    n_dropped_rows <- n_rows - nrow(phys_data)
    message(format(n_dropped_rows, big.mark = ','), ' rows were dropped.')
    return(phys_data)
}

#' Get data ready for propagation
#'
#' @param tbl A data.frame.
#'
#' @return A data.frame.
#' @export
#'
getDataReady <- function(tbl) {
    if (is.null(tbl)) {
        return(NULL)
    }

    if (!nrow(tbl)) {
        return(NULL)
    }

    attr_type <- unique(tbl$Attribute_type)
    if (attr_type == 'binary') {
        set_with_ids <- getSetWithIDs(tbl) |>
                purrr::discard(~ all(is.na(.x)))
        set_without_ids <- getSetWithoutIDs(tbl, set_with_ids) |>
                purrr::discard(~ all(is.na(.x)))
        if (is.null(set_with_ids) && is.null(set_without_ids))
            return(NULL)
        dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
        output <- completeBinaryData(dataset)
    } else if (attr_type == 'multistate-intersection') {
        set_with_ids <- getSetWithIDs(tbl) |>
                purrr::discard(~ all(is.na(.x)))
        set_without_ids <- getSetWithoutIDs(tbl, set_with_ids = set_with_ids) |>
                purrr::discard(~ all(is.na(.x)))
        if (is.null(set_with_ids) && is.null(set_without_ids))
            return(NULL)
        dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
        output <- dataset |>
            tidyr::complete(NCBI_ID, Attribute, fill = list(Score = 0)) |>
            dplyr::arrange(NCBI_ID, Attribute)
    } else if (attr_type == 'range') { # all numeric are converted to range when imported with the physiologies function
        set_with_ids <- getSetWithIDs(tbl) |>
            purrr::discard(~ all(is.na(.x)))
        set_without_ids <- getSetWithoutIDs(tbl, set_with_ids = set_with_ids) |>
            purrr::discard(~ all(is.na(.x)))
        if (is.null(set_with_ids) && is.null(set_without_ids))
            return(NULL)
        dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
        output <- dataset |>
            dplyr::arrange(NCBI_ID)
    } else if (attr_type == 'multistate-union') {
        tbl$Attribute_group_2 <- sub('--(TRUE|FALSE)', '', tbl$Attribute)
        l <- split(tbl, factor(tbl$Attribute_group_2))
        output <- vector('list', length(l))
        for (i in seq_along(output)) {
            set_with_ids <- getSetWithIDs(l[[i]]) |>
                purrr::discard(~ all(is.na(.x)))
            set_without_ids <- getSetWithoutIDs(l[[i]], set_with_ids) |>
                purrr::discard(~ all(is.na(.x)))
            dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
            if (!nrow(dataset))
                next
            names(output)[i] <- names(l)[i]
            output[[i]] <- completeBinaryData(dataset)
        }
    }
    return(output)
}

#' Complete binary data
#'
#' \code{completeBinaryData} completes Attributes and scores for binary
#' attributes
#'
#' @param tbl A data.frame.
#'
#' @return A data.frame.
#'
#' @export
#'
completeBinaryData <- function(tbl) {
        current_attr <- unique(tbl$Attribute)
        if (length(current_attr) == 1 && grepl('--TRUE$', current_attr)) {
            extra_level <- sub('--TRUE$', '--FALSE', current_attr)
            tbl$Attribute <- factor(tbl$Attribute, levels = c(extra_level, current_attr))
        }
        tbl |>
            tidyr::complete(.data$NCBI_ID, .data$Attribute, fill = list(Score = 0)) |>
            dplyr::arrange(.data$NCBI_ID, .data$Attribute)
}

#' Get set with IDs
#'
#' \code{getSetWithIDs}
#'
#' @param tbl A data.frame.
#'
#' @return A data.frame.
#' @export
#'
getSetWithIDs <- function(tbl) {

    ## Check for taxa with NCBI IDs
    lgl_vct <- is.na(tbl$NCBI_ID) | tbl$NCBI_ID == 'unknown'
    tbl <- tbl |>
        dplyr::filter(!lgl_vct)
    if (!nrow(tbl))
        return(NULL)

    ## Check for taxa with valid ranks
    valid_ranks <- c('genus', 'species', 'strain')
     tbl <- tbl |>
        dplyr::mutate(
            Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')
        ) |>
        dplyr::filter(.data$Rank %in% valid_ranks)
    if (!nrow(tbl))
        return(NULL)

     ## Check taxa names, resolve conflicts, normalize scores
     if (unique(tbl$Attribute_type) == "range") {
        output <- tbl |>
            dplyr::mutate(
                Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')
            ) |>
            dplyr::distinct() |>
            dplyr::mutate(
                Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)
            ) |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::slice_max(
                .data$Confidence_in_curation, n = 1, with_ties = TRUE
            ) |>
            dplyr::mutate(
                Attribute_value = mean(cbind(.data$Attribute_value_min, .data$Attribute_value_max))
            ) |>
            # dplyr::mutate(Attribute_value = mean(Attribute_value)) |>
            dplyr::slice_max(
                .data$Confidence_in_curation, n = 1, with_ties = FALSE
            ) |>
            dplyr::ungroup() |>
            dplyr::select(-Attribute_value_min, -Attribute_value_max) |>
            dplyr::distinct() |>
            dplyr::select(-.data$Parent_NCBI_ID) |>
            dplyr::mutate(taxid = .data$NCBI_ID) |>
            dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
            dplyr::filter(!is.na(.data$NCBI_ID)) |>
            dplyr::distinct()
            # dplyr::relocate(tidyselect::all_of(.orderedColumns()))
     } else {
        output <- tbl |>
            dplyr::mutate(
                Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')
            ) |>
            dplyr::distinct() |>
            dplyr::mutate(
                Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)
            ) |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::slice_max(
                .data$Confidence_in_curation, n = 1, with_ties = TRUE
            ) |>
            dplyr::ungroup() |>
            dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
            dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate(
                Total_score = sum(.data$Score),
                Score = .data$Score / .data$Total_score
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(Frequency = scores2Freq(.data$Score)) |>
            dplyr::select(-.data$Parent_NCBI_ID, -.data$Total_score) |>
            dplyr::mutate(taxid = .data$NCBI_ID) |>
            dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
            dplyr::filter(!is.na(.data$NCBI_ID)) |>
            dplyr::distinct() |>
            dplyr::relocate(tidyselect::all_of(.orderedColumns()))
     }

    return(output)
}

#' Get set without IDs
#'
#' \code{getSetWithoutIDs}
#'
#' @param tbl A data.frame.
#' @param set_with_ids  A data.frame
#'
#' @return A data.frame.
#' @export
#'
getSetWithoutIDs <- function(tbl, set_with_ids = NULL) {
    attribute_type_var <- unique(tbl$Attribute_type)
    attribute_group_var <- unique(tbl$Attribute_group)

    if (is.null(set_with_ids))
        set_with_ids <- data.frame(NCBI_ID = 'not a real NCBI_ID')

    valid_ranks <- c('genus', 'species', 'strain')
    lgl_vct <- is.na(tbl$NCBI_ID) | tbl$NCBI_ID == 'unknown'
    if (!any(lgl_vct))
        return(NULL)

    if (unique(tbl$Attribute_type == "range")) {
        output <- tbl |>
            dplyr::filter(lgl_vct) |>
            dplyr::select(
                -.data$NCBI_ID, -.data$Taxon_name, -.data$Frequency
            ) |>
            dplyr::relocate(NCBI_ID = .data$Parent_NCBI_ID) |>
            dplyr::distinct() |>
            dplyr::mutate(Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')) |>
            dplyr::filter(Rank %in% valid_ranks) |>
            dplyr::mutate(Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')) |>
            dplyr::mutate(Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)) |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::slice_max(.data$Confidence_in_curation, n = 1, with_ties = TRUE) |>
            dplyr::mutate(
                Attribute_value = mean(cbind(.data$Attribute_value_min, .data$Attribute_value_max))
            ) |>
            dplyr::slice_max(
                .data$Confidence_in_curation, n = 1, with_ties = FALSE
            ) |>
            dplyr::ungroup() |>
            dplyr::select(-Attribute_value_min, -Attribute_value_max) |>
            dplyr::distinct() |>
            dplyr::mutate(Frequency = scores2Freq(.data$Score)) |>
            dplyr::mutate(
                Evidence = 'tax',
                Attribute_group = attribute_group_var,
                Attribute_type = attribute_type_var,
                Attribute_source = NA,
                Confidence_in_curation = NA
            ) |>
            dplyr::mutate(taxid = .data$NCBI_ID) |>
            dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
            dplyr::filter(!is.na(.data$NCBI_ID)) |>
            dplyr::filter(!.data$NCBI_ID %in% unique(set_with_ids$NCBI_ID)) |>
            dplyr::distinct() |>
            dplyr::arrange(.data$NCBI_ID)
            # dplyr::relocate(tidyselect::all_of(.orderedColumns()))
    } else {
        output <- tbl |>
            dplyr::filter(lgl_vct) |>
            dplyr::select(
                -.data$NCBI_ID, -.data$Frequency
            ) |>
            dplyr::relocate(NCBI_ID = .data$Parent_NCBI_ID) |>
            dplyr::distinct() |>
            dplyr::mutate(Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')) |>
            dplyr::filter(Rank %in% valid_ranks) |>
            dplyr::mutate(Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)) |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::slice_max(.data$Confidence_in_curation, n = 1, with_ties = TRUE) |>
            dplyr::ungroup() |>
            dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
            dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = TRUE) |>
            dplyr::ungroup() |>
            dplyr::group_by(.data$NCBI_ID) |>
            dplyr::mutate(
                Total_score = sum(.data$Score),
                Score = .data$Score / .data$Total_score
            ) |>
            dplyr::ungroup() |>
            dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
            dplyr::mutate(Score = sum(.data$Score, na.rm = TRUE)) |>
            dplyr::ungroup() |>
            dplyr::select(-.data$Taxon_name) |>
            dplyr::distinct() |>
            dplyr::mutate(Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')) |>
            dplyr::mutate(Frequency = scores2Freq(.data$Score)) |>
            dplyr::mutate(
                Evidence = 'tax',
                Attribute_group = attribute_group_var,
                Attribute_type = attribute_type_var,
                Attribute_source = NA,
                Confidence_in_curation = NA
            ) |>
            dplyr::ungroup() |>
            dplyr::select(-.data$Total_score) |>
            dplyr::mutate(taxid = .data$NCBI_ID) |>
            dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
            dplyr::filter(!is.na(.data$NCBI_ID)) |>
            dplyr::filter(!.data$NCBI_ID %in% unique(set_with_ids$NCBI_ID)) |>
            dplyr::distinct() |>
            dplyr::arrange(.data$NCBI_ID, .data$Attribute) |>
            dplyr::relocate(tidyselect::all_of(.orderedColumns()))
    }
    return(output)
}

.orderedColumns <- function() {
    c(
        'NCBI_ID', 'Taxon_name', 'Rank',
        'Attribute', 'Attribute_source', 'Confidence_in_curation',
        'Evidence', 'Frequency', 'Score',
        'Attribute_group', 'Attribute_type',
        'taxid'
    )
}
