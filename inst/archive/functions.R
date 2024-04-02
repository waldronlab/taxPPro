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

#' Clean node
#'
#' \code{cleanNode} deletes the `attribute_tbl` slot in a data.tree node.
#'
#' @param node A node in a data.tree object.
#'
#' @return NULL (assigned to the node)
#' @export
#'
cleanNode <- function(node) {
    node$attribute_tbl <- NULL
}
