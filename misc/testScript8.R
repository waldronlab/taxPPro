library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)
phys <- physiologies(remove_false = TRUE, full_source = FALSE)
data1 <- map(phys, prepareData)
data1 <- discard(data1, is.null)
data2 <- map(data1, prepareData2)
## isolation site must be filtered or it would take too long
select_attrs <- data2$`isolation site`|>
    count(Attribute) |>
    arrange(-n) |>
    filter(n >= 10) |>
    pull(Attribute)
data2$`isolation site` <- data2$`isolation site` |>
    filter(Attribute %in% select_attrs)
data('tree_list')
tree <- data.tree::as.Node(tree_list)
trees <- vector('list', length(data2))
for (i in seq_along(trees)) {
    message('Propagating ', names(data2)[i], ' - ', Sys.time())
    names(trees)[i] <- names(data2)[i]
    trees[[i]] <- propagate(tree, data2[[i]])
}
ncbi_taxonomy <- get_ncbi_taxonomy()
dfs <- map(trees, ~ toDataFrame(.x, ncbi_taxonomy))
dfs2 <- map(trees, ~ {
    args <- as.list(.x$attributesAll)
    args <- c(list(x = .x, row.names = NULL, optional = FALSE), args)
    df <- do.call('as.data.frame', args)
    return(df)
})
