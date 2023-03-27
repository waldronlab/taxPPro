
library(taxPPro)
library(bugphyzz)
library(purrr)
library(dplyr)
library(data.tree)

phys <- physiologies(remove_false = TRUE, full_source = TRUE)
x <- map(phys, ~ {
    tryCatch(
        error = function(e) e, getDataReadyForPropagation(.x)
    )
})

any(map_lgl(x, rlang::is_error))
data('tree_list')
tree <- as.Node(tree_list)

propagated <- vector('list', length(data_ready))
for (i in seq_along(propagated)) {
    message('Propagating ', names(data_ready)[i], ' - ', Sys.time())
    names(propagated)[i] <- names(data_ready)[i]
    propagated[[i]] <- tryCatch(
        error = function(e) e,
        {
            propagate(data_tree = tree, df = data_ready[[i]])
        }
    )
}



















