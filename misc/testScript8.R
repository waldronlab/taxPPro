library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)
fname <- system.file('extdata/links.tsv', package = 'bugphyzz')
data <- read.table(fname, row.names = NULL, header = TRUE, sep = '\t')
data <- data[,c('physiology', 'sig_type')]
data <- data[data$sig_type == 'numeric', ]
phys_names <- data[['physiology']]
phys <- physiologies(phys_names, remove_false = TRUE, full_source = FALSE)
data1 <- map(phys, prepareData)
data1 <- discard(data1, is.null)
data2 <- map(data1, prepareData2)
data('tree_list')
tree <- data.tree::as.Node(tree_list)
trees <- vector('list', length(data2))
for (i in seq_along(trees)) {
    message('Propagating ', names(data2)[i])
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
