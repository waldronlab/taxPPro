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

#' Get living tree project (LTP) tree or data
#'
#' \code{ltp} gets the LTP tree or data.
#'
#' @param x A character string. tree, tips, or nodes. Default is tree.
#' @param remove_zero_tips Remove or not tips with zero branch lengths.
#' These tips could cause trouble.
#'
#' @return A phylo or data.frame (see the `x` param)
#' @export
#'
ltp <- function(x = 'tree', remove_zero_tips = TRUE) {
    tree_fname <- system.file(
        'extdata', 'livingTree.newick', package = 'taxPPro'
    )
    tree <- ape::read.tree(tree_fname)
    if (remove_zero_tips) {
        ntips1 <- length(tree$tip.label)
        zero_edges <- which(tree$edge.length == 0)
        nodes_with_zero_edges <- tree$edge[zero_edges, 2]
        remove_tips <- nodes_with_zero_edges[nodes_with_zero_edges <= length(tree$tip.label)]
        remove_tips_names <- tree$tip.label[remove_tips]
        tree <- drop.tip(tree, remove_tips_names)
        ntips2 <- length(tree$tip.label)
        message('Dropping ', ntips1 - ntips2, ' tips with zero length branches.')
        message(ntips2, ' tips remaining in the tree.')
    }
    if (x == 'tree') {
        output <- tree
    } else if (x == 'tips') {
        fname <- system.file(
            'extdata', 'livingTree_tips.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        ) |>
            purrr::modify(as.character)
        rownames(output) <- output$tip_label
        output <- output[tree$tip.label,]
    } else if (x == 'nodes') {
        fname <- system.file(
            'extdata', 'livingTree_nodes.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        ) |>
            purrr::modify(as.character)
    }
    return(output)
}
