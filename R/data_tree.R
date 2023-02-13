
## A function for removing prefix from a character vector x
removePrefix <- function(x) {
    regex <- '^[kpcofgst]__'
    sub(regex, '', x)
}

#' Add child taxa
#'
#' @param node A node.
#'
#' @return An action on the node.
#' @export
#'
addChildren <- function(node) {

    if (node$isLeaf) {
        current_children <- removePrefix(names(node$children))
        current_taxon <- removePrefix(node$name)
        df <- tryCatch(
            error = function(e) NULL, {
                taxizedb::children(current_taxon)[[1]]
            }
        )
        if (is.null(df)) {
            return(NULL)
        }
        df <- df[df$rank %in% 'strain', ]
        children <- df$name
        new_children <- children[!children %in% current_children]
        if (isFALSE(!length(new_children))) {
            for (i in seq_along(new_children)) {
                message(new_children)
                node$AddChild(new_children[i])
            }
        }
    }

}
