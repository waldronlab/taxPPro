#' Filter data
#'
#' \code{filterData}
#'
#' @param tbl A data.frame.
#'
#' @return A data.frame
#' @export
#'
filterData <- function(tbl) {
    types <- bugphyzz:::.DISCRETE_ATTRIBUTE_TYPES()
    attr_type <- unique(tbl$Attribute_type)
    if (attr_type %in% types){
        output <- filterDataDiscrete(tbl)
    } else {
        output <- NULL
    }
    return(output)
}

#' Filter attributes with discrete type
#'
#' \code{filterDataDiscrete} filters discrete data.
#'
#' @param tbl A data.frame.
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
    # valid_ranks <- c('genus', 'species', 'strain')
    attributes_fname <- system.file(
        'extdata', 'attributes.tsv', package = 'bugphyzz'
    )
    attributes <- utils::read.table(attributes_fname, sep = '\t', header = TRUE)
    rgx <- paste0('\\b', phys_name, '\\b')
    valid_attributes <- attributes |>
        # dplyr::filter(.data$attribute_group == phys_name) |>
        dplyr::filter(grepl(rgx, .data$attribute_group)) |>
        dplyr::pull(.data$attribute) |>
        unique()

    n_rows <- nrow(tbl)

    tbl <- tbl |>
        dplyr::filter(!is.na(.data$Attribute)) |>
        dplyr::filter(!is.na(.data$Attribute_value))

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
    attr_type <- unique(tbl$Attribute_type)
    if (attr_type == 'binary') {
        set_with_ids <- getSetWithIDs(tbl) |>
                purrr::discard(~ all(is.na(.x)))
        set_without_ids <- getSetWithoutIDs(tbl, set_with_ids) |>
                purrr::discard(~ all(is.na(.x)))
        dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
        output <- completeBinaryData(dataset)
    } else if (attr_type == 'multistate-intersection') {
        set_with_ids <- getSetWithIDs(tbl) |>
                purrr::discard(~ all(is.na(.x)))
        set_without_ids <- getSetWithoutIDs(tbl, set_with_ids = set_with_ids) |>
                purrr::discard(~ all(is.na(.x)))
        dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
        output <- dataset |>
            tidyr::complete(NCBI_ID, Attribute, fill = list(Score = 0)) |>
            dplyr::arrange(NCBI_ID, Attribute)
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
            # output[[i]] <- completeBinaryData(dataset)
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
        # dplyr::filter(!is.na(.data$NCBI_ID)) |>
        dplyr::distinct() |>
        dplyr::arrange(.data$NCBI_ID, .data$Attribute) |>
        dplyr::relocate(tidyselect::all_of(.orderedColumns()))

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

    output <- tbl |>
        dplyr::filter(lgl_vct) |>
        dplyr::select(
            -.data$NCBI_ID, -.data$Taxon_name, -.data$Frequency
        ) |>
        dplyr::relocate(NCBI_ID = .data$Parent_NCBI_ID) |>
        distinct() |>
        dplyr::mutate(Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')) |>
        dplyr::filter(Rank %in% valid_ranks) |>
        dplyr::mutate(Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')) |>
        dplyr::mutate(Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(.data$Confidence_in_curation, n = 1, with_ties = TRUE) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
        dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::mutate(
            Total_score = sum(.data$Score),
            Score = .data$Score / .data$Total_score
        ) |>
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

#' Perform taxonomic pooling
#'
#' \code{taxPool} performs taxonomic pooling at the strain, species, and
#' genus ranks. Only use in the data.tree NCBI tree.
#'
#' @param node A node in a data.tree object.
#' @param grp  A character string. Attribute group.
#' @param typ A character string. Attribute_type.
#'
#' @return A table for the attribute of a node.
#' @export
#'
taxPool <- function(node, grp, typ) {
    if (!node$isLeaf) {
        children_names <- names(node$children)
        attribute_tbls <- children_names |>
            purrr::map(~ node[[.x]]$attribute_tbl) |>
            purrr::discard(is.null)
        not_all_children_tbls_are_null <- length(attribute_tbls) > 0
        node_attribute_tbl_is_null <- is.null(node$table)
        # node_is_gst <- grepl('^[gst]__', node$name)
        node_is_gst <- grepl('^[st]__', node$name)
        conds <- node_attribute_tbl_is_null &
            not_all_children_tbls_are_null &
            node_is_gst
        if (conds) {
            res_tbl <- attribute_tbls |>
                purrr::discard(is.null) |>
                dplyr::bind_rows() |>
                dplyr::select(
                    .data$NCBI_ID, .data$Attribute, .data$Score
                ) |>
                dplyr::mutate(
                    NCBI_ID = node$name,
                    taxid = node$taxid,
                    Taxon_name = node$Taxon_name,
                    Rank = node$Rank,
                    Evidence = 'tax',
                    Attribute_group = grp,
                    Attribute_type = typ
                ) |>
                dplyr::group_by(.data$NCBI_ID) |>
                dplyr::mutate(
                    Total_score = sum(.data$Score),
                    Score = .data$Score / .data$Total_score
                ) |>
                dplyr::ungroup() |>
                dplyr::select(-.data$Total_score) |>
                dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
                dplyr::mutate(Score = sum(.data$Score)) |>
                dplyr::ungroup() |>
                dplyr::distinct() |>
                dplyr::mutate(
                    Frequency = dplyr::case_when(
                        .data$Score == 1 ~ 'always',
                        .data$Score > 0.9 ~ 'usually',
                        .data$Score >= 0.5 ~ 'sometimes',
                        .data$Score > 0 & .data$Score < 0.5 ~ 'rarely',
                        .data$Score == 0 ~ 'never'
                    )
                ) |>
                dplyr::mutate(
                    Attribute_source = NA,
                    Confidence_in_curation = NA
                ) |>
                dplyr::distinct()
            node$attribute_tbl <- res_tbl
        }
    }
}

#' Inheritance first round
#'
#' \code{inh1} First round of inheritance. Only use with Do in a data.tree
#' object.
#'
#' @param node  A node in a data.tree.
#' @param adjF Adjustment factor for penalty. Default is 1.0.
#' @param evidence_label A character string, either inh1 or inh2.
#'
#' @return A table for a node attribute.
#' @export
#'
inh1 <- function(node,  adjF = 0.1, evidence_label = 'inh') {
    if (node$isRoot)
        return(NULL)
    if (is.null(node$parent$attribute_tbl))
        return(NULL)
    if (is.null(node$attribute_tbl) && grepl('^[st]__', node$name)) {
        df <- node$parent$attribute_tbl
        n <- nrow(df)
        df <- df |>
            dplyr::mutate(
                target_scores = rep(1 / n, n),
                score_diff = .data$Score - .data$target_scores,
                Score = .data$Score - adjF * .data$score_diff,
                NCBI_ID = node$name,
                Evidence = evidence_label,
                Taxon_name = node$Taxon_name,
                Rank = node$Rank,
                taxid = node$taxid,
            ) |>
            dplyr::select(-.data$target_scores, -.data$score_diff) |>
            dplyr::relocate(.data$NCBI_ID) |>
            dplyr::mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                taxid = node$taxid,
                Taxon_name = node$Taxon_name,
                Rank = node$Rank,
                Frequency = dplyr::case_when(
                    .data$Score == 1 ~ 'always',
                    .data$Score > 0.9 ~ 'usually',
                    .data$Score >= 0.5 ~ 'sometimes',
                    .data$Score > 0 & .data$Score < 0.5 ~ 'rarely',
                    .data$Score == 0 ~ 'never'
                )
            )
        node$attribute_tbl <- df
    }
}

#' Inheritance second round
#'
#' \code{inh2} Second round of inheritance. Only use with Do in a data.tree
#' object.
#'
#' @param node  A node in a data.tree.
#' @param adjF Adjustment factor for penalty. Default is 1.0.
#'
#' @return A table for a node attribute.
#' @export
#'
inh2 <- function(node, adjF = 0.1) {
    cond1 <- !node$isRoot
    cond2 <- is.null(node$attribute_tbl)
    cond3 <- !is.null(node$parent$attribute_tbl)
    if (cond1 && cond2 && cond3) {
        tbl <- node$parent$attribute_tbl
        n <- nrow(tbl)
        res <- tbl |>
            dplyr::mutate(
                target_scores = rep(1 / n, n),
                score_diff = .data$Score - .data$target_scores,
                Score = .data$Score - adjF * .data$score_diff,
                NCBI_ID = node$name,
                Evidence = 'inh2'
            ) |>
            dplyr::select(-.data$target_scores, -.data$score_diff) |>
            dplyr::relocate(.data$NCBI_ID) |>
            dplyr::mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                taxid = node$taxid,
                Taxon_name = node$Taxon_name,
                Rank = node$Rank,
                Frequency = dplyr::case_when(
                    .data$Score == 1 ~ 'always',
                    .data$Score > 0.9 ~ 'usually',
                    .data$Score >= 0.5 ~ 'sometimes',
                    .data$Score > 0 & .data$Score < 0.5 ~ 'rarely',
                    .data$Score == 0 ~ 'never'
                )
            )
        node$attribute_tbl <- res
    }
}

#' Get Most Recent Common Ancestor
#'
#' \code{getMRCATaxPPro} gets the most recent common ancestor using phytools.
#'
#' @param tree A phylo tree object.
#' @param tips Tips in the phylo tree object.
#'
#' @return A data.frame.
#' @export
#'
getMRCATaxPPro <- function(tree, tips) {
    res <- phytools::findMRCA(tree = tree, tips = tips)
    if (is.null(res))
        res <- NA
    res
}

