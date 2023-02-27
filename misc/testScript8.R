library(bugphyzz)
library(purrr)
phys <- physiologies(keyword = 'all', remove_false = TRUE, full_source = FALSE)
data1 <- map(phys, prepareData)
data1 <- discard(data1, is.null)
data2 <- map(data1, prepareData2)
phys_names <- c('aerophilicity', 'width', 'acetate producing')
data2 <- data2[phys_names]
data('tree_list')
tree <- data.tree::as.Node(tree_list)
trees <- vector('list', length(data2))
for (i in seq_along(trees)) {
    message('Propagating ', names(data2)[i])
    trees[[i]] <- propagate(tree, data2[[i]])
}
dfs <- map(trees, toDataFrame)
