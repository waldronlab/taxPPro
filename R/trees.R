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
#'
#' @return A phylo or data.frame (see the `x` param)
#' @export
#'
ltp <- function(x = 'tree') {
    if (x == 'tree') {
        fname <- system.file(
            'extdata', 'livingTree.newick', package = 'taxPPro'
        )
        output <- ape::read.tree(fname)
    } else if (x == 'tips') {
        fname <- system.file(
            'extdata', 'livingTree_tips.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        )
    } else if (x == 'nodes') {
        fname <- system.file(
            'extdata', 'livingTree_nodes.tsv', package = 'taxPPro'
        )
        output <- utils::read.table(
            file = fname, header = TRUE, sep = '\t',
            row.names = NULL
        )
    }
    return(output)
}
