library(bugphyzz)
library(taxPPro)
library(purrr)
library(rlang)

phys_names <- c('aerophilicity', 'growth temperature')
phys <- physiologies(phys_names, remove_false = TRUE, full_source = FALSE)
data_ready <- vector('list', length(phys))
for (i in seq_along(data_ready)) {
    message('Preparing ', names(phys)[i])
    names(data_ready)[i] <- names(phys)[i]
    data_ready[[i]] <- tryCatch(
        error = function(e) e,
        {
            prepareDatForPropagation(phys[[i]])
        }
    )
}
vct_lgl <- map_lgl(data_ready, is_error)
if (any(vct_lgl))  {
    message('Removing data with errors.')
    data_ready <- discard(data_ready, is_error)
}






