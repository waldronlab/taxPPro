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
#' @param remove_zero_tips Remove or not tips with zero branch lengths.
#' These tips could cause trouble.
#'
#' @return A list with a phylo (the ltp tree) and a data.frame (tip_data)
#' objects.
#' @export
#'
ltp <- function(remove_zero_tips = TRUE) {
    tree_fname <- system.file(
        'extdata', 'livingTree.newick', package = 'taxPPro'
    )
    tip_data_fname <- system.file(
        'extdata', 'livingTree_tips.tsv', package = 'taxPPro'
    )
    node_data_fname <- system.file(
        'extdata', 'livingTree_nodes.tsv', package = 'taxPPro'
    )
    tree <- ape::read.tree(tree_fname)
    message('Initial number of tips in LTP tree: ', length(tree$tip.label))
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

    if (remove_zero_tips) {
        ntips1 <- length(tree$tip.label)
        zero_edges <- which(tree$edge.length == 0)
        nodes_with_zero_edges <- tree$edge[zero_edges, 2]
        remove_tips <- nodes_with_zero_edges[nodes_with_zero_edges <= length(tree$tip.label)]
        remove_tips_names <- tree$tip.label[remove_tips]
        tree <- ape::drop.tip(phy = tree, tip = remove_tips_names)
        ntips2 <- length(tree$tip.label)
        tip_data <- tip_data[tree$tip.label,]
        message('Dropping ', ntips1 - ntips2, ' tips with zero length branches.')
    }

    ntips3 <- nrow(tip_data)
    tip_data <- tip_data |>
        dplyr::group_by(taxid) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        as.data.frame()
    rownames(tip_data) <- tip_data$tip_label
    ntips4 <- nrow(tip_data)
    tree <- ape::keep.tip(phy = tree, tip = tip_data$tip_label)
    message('Dropping ', ntips3 - ntips4, ' tips because of duplicated taxids.')
    tip_data <- tip_data[tree$tip.label,]
    message('Tips remaining: ', length(tree$tip.label))
    list(
        tree = tree, tip_data = tip_data, node_data = node_data
    )
    # if (x == 'tree') {
    #     output <- tree
    # } else if (x == 'tips') {
    #     fname <- system.file(
    #         'extdata', 'livingTree_tips.tsv', package = 'taxPPro'
    #     )
    #     output <- utils::read.table(
    #         file = fname, header = TRUE, sep = '\t',
    #         row.names = NULL
    #     ) |>
    #         purrr::modify(as.character)
    #     rownames(output) <- output$tip_label
    #     output <- output[tree$tip.label,]
    # } else if (x == 'nodes') {
    #     fname <- system.file(
    #         'extdata', 'livingTree_nodes.tsv', package = 'taxPPro'
    #     )
    #     output <- utils::read.table(
    #         file = fname, header = TRUE, sep = '\t',
    #         row.names = NULL
    #     ) |>
    #         purrr::modify(as.character)
    # }
    # return(output)
}




ltp2 <- function(remove_zero_tips = TRUE) {
    tree_fname <- system.file(
        'extdata', 'livingTree.newick', package = 'taxPPro'
    )
    tip_data_fname <- system.file(
        'extdata', 'livingTree_tips.tsv', package = 'taxPPro'
    )
    node_data_fname <- system.file(
        'extdata', 'livingTree_nodes.tsv', package = 'taxPPro'
    )

    tree <- ape::read.tree(tree_fname)
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

    message('Initial number of tips in LTP tree: ', length(tree$tip.label))


    if (remove_zero_tips) {

        ntips1 <- length(tree$tip.label)
        remove_tips <- 0
        while(length(remove_tips) > 0) {
            zero_edges <- which(tree$edge.length == 0)
            nodes_with_zero_edges <- tree$edge[zero_edges, 2]
            remove_tips <- nodes_with_zero_edges[nodes_with_zero_edges <= length(tree$tip.label)]
            remove_tips_names <- tree$tip.label[remove_tips]
            tree <- ape::drop.tip(phy = tree, tip = remove_tips_names, trim.internal = FALSE)
            message('an interation complete')

        }
        message('out of the loop')
        ntips2 <- length(tree$tip.label)
        tip_data <- tip_data[tree$tip.label,]
        message('Dropping ', ntips1 - ntips2, ' tips with zero length branches.')
    }

    ntips3 <- nrow(tip_data)
    tip_data <- tip_data |>
        dplyr::group_by(taxid) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        as.data.frame()
    rownames(tip_data) <- tip_data$tip_label
    ntips4 <- nrow(tip_data)
    tree <- ape::keep.tip(phy = tree, tip = tip_data$tip_label)
    message('Dropping ', ntips3 - ntips4, ' tips because of duplicated taxids.')
    tip_data <- tip_data[tree$tip.label,]
    message('Tips remaining: ', length(tree$tip.label))
    list(
        tree = tree, tip_data = tip_data
    )
}











