
## A function for removing prefix from a character vector x
removePrefix <- function(x) {
    regex <- '^[kpcofgst]__'
    sub(regex, '', x)
}

#' Add strains
#'
#' @param node A node.
#'
#' @return An action on the node.
#' @export
#'
addStrains <- function(node) {

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
        children <- df$id
        new_children <- children[!children %in% current_children]
        new_children <- paste0('t__', new_children)
        if (isFALSE(!length(new_children))) {
            for (i in seq_along(new_children)) {
                message(new_children)
                node$AddChild(new_children[i])
            }
        }
    }

}
