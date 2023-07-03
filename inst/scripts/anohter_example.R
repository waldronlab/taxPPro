library(data.tree)
library(bugphyzz)
library(dplyr)
library(taxPPro)

data("tree_list")
my_tree <- as.Node(tree_list)
tree <- my_tree$clone()

aer <- physiologies('aerophilicity', remove_false = TRUE, full_source = FALSE)[[1]] |>
    # mutate(
    #     letter = sub('^(.).*', '\\1', Rank),
    #     NCBI_ID = paste0(letter, '__', NCBI_ID)
    # ) |>
    prepareDatForPropagation()
l <- split(aer, factor(aer$NCBI_ID))
l[['g__561']] <- NULL


tree$Do(function(node) {
    node[['table']] <- l[[node$name]]
})

print(tree$d__2$p__1224$c__1236$o__91347$f__543$g__561, 'table')
tree$d__2$p__1224$c__1236$o__91347$f__543$g__561[['s__562']]$table
tree$d__2$p__1224$c__1236$o__91347$f__543$g__561$table


calcParent <- function(tbl) {
    tbl |>
        dplyr::group_by(NCBI_ID) |>
        dplyr::mutate(Score =  Score / sum(Score)) |>
        dplyr::mutate(Score = ifelse(is.na(Score), 0, Score)) |>
        dplyr::group_by(Attribute) |>
        dplyr::reframe(
            Score = sum(Score),
            Evidence = paste0(Evidence, collapse = '|'),
        ) |>
        dplyr::mutate(Score = Score / sum(Score)) |>
        dplyr::mutate(Score = ifelse(is.na(Score), 0, Score)) |>
        dplyr::mutate(Evidence = sub('^(\\|*)', '', Evidence)) |>
        dplyr::mutate(Evidence = sub('\\|\\|+', '|', Evidence))
}


myFun <- function(node) {
    if (is.null(node[['table']])) {
        message(node$name)
        children <- names(node$children)
        output <- vector('list', length(children))
        for (i in seq_along(output)) {
            output[[i]] <- node[[children[i]]]$table
        }
        cond <- purrr::map_lgl(output, is.null)
        if (!all(cond)) {
            df <- purrr::discard(output, is.null) |>
                dplyr::bind_rows() |>
                dplyr::select(
                    NCBI_ID, Attribute, Score, Evidence #, Attribute_source
                ) |>
                calcParent()
            node[['table']] <- df
        }
    }
}

tree$Do(myFun, traversal = 'post-order')
