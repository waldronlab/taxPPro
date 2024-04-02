
#' Get living tree project (LTP) tree and data
#'
#' \code{ltp} Imports a modified version fo the phylogenetic tree of the
#' Living Tree Project.
#'
#' @param remove_gn_nodes Remove rows from node_data with genus taxids.
#' Default is TRUE, which removes them. These are already in tip_data.
#' @param node_names Add names to unnamed nodes (n + node number).
#' This will only be in the tree (not in node_data).
#'
#' @return A list with a the LTP tree (phylo), tips and node data (data.frames),
#' and the names of the genera added to the original tree (character vector).
#' .
#' @export
#'
ltp <- function(remove_gn_nodes = TRUE, node_names = TRUE) {
    tree_fname <- system.file(
        'extdata', 'LTP_all_08_2023.newick', package = 'taxPPro'
    )
    tip_data_fname <- system.file(
        'extdata', 'LTP_all_08_2023.tip_data', package = 'taxPPro'
    )
    node_data_fname <- system.file(
        'extdata', 'LTP_all_08_2023.node_data', package = 'taxPPro'
    )
    tree <- ape::read.tree(tree_fname)

    if (node_names) {
        tree$node.label <- ifelse(
            tree$node.label == "NA",
            paste0("n", ape::Ntip(tree) + 1:(ape::Nnode(tree))),
            tree$node.label
        )
    }

    tip_data <- utils::read.table(
        file = tip_data_fname, header = TRUE, sep = '\t', row.names = NULL
    ) |>
        purrr::modify(as.character) |>
        as.data.frame()
    rownames(tip_data) <- tip_data$tip_label

    node_data <- utils::read.table(
        file = node_data_fname, header = TRUE, sep = '\t', row.names = NULL
    ) |>
        purrr::modify(as.character) |>
        as.data.frame()

    gn_tips <- grep('g__', tip_data$tip_label, value = TRUE)

    if (remove_gn_nodes) {
        node_data <- node_data |>
            dplyr::filter(Rank != 'genus') |>
            purrr::discard(~all(is.na(.x)))
    }

    list(
        tree = tree, tip_data = tip_data, node_data = node_data,
        gn_tips = gn_tips
    )
}

#' Get NSTI
#'
#' \code{getNsti} gets the nsti described in picrust2. Code is based
#' on the source code of picrust2 using the castor package. NSTI
#' values are described in picrust and picrust2.
#'
#' @param tree A phylo object
#' @param annotated_tip_labels A character vector with the names of the tips
#'
#' @return A data.frame with two columns: tip_label and nsti.
#' @export
#'
getNsti <- function(tree, annotated_tip_labels) {
    unknown_tips_index <- which(!tree$tip.label %in% annotated_tip_labels)
    known_tips_index <- which(tree$tip.label %in% annotated_tip_labels)
    res <- castor::find_nearest_tips(
        tree = tree, target_tips = known_tips_index
    )
    data.frame(
        tip_label = tree$tip.label[unknown_tips_index],
        nsti = res$nearest_distance_per_tip[unknown_tips_index]
    )
}

