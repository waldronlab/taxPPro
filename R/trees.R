#' Get metphlan tree or data v31
#'
#' \code{mpa} gets the metaphlan tree or data v31.
#'
#' @param x A character string: tree, tips, or nodes. Default is tree.
#'
#' @return A phylo or data.frame (see the `x` param)
#' @export
#'
mpa <- function(x = 'tree') {
    if (x == 'tree') {
        fname <- system.file(
            'extdata', 'mpav31.newick', package = 'taxPPro'
        )
        output <- ape::read.tree(fname)
    } else if (x == 'tips') {
        fname <- system.file(
            'extdata', 'mpav31_tips.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        )
    } else if (x == 'nodes') {
        fname <- system.file(
            'extdata', 'mpav31_nodes.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        )
    }
    return(output)
}

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
        tree$node.label <- ifelse(tree$node.label == "NA", paste0("n", ape::Ntip(tree) + 1:(ape::Nnode(tree))), tree$node.label)
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
