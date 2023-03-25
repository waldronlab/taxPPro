
library(bugphyzz)
library(purrr)
library(dplyr)

phys <- physiologies(remove_false = TRUE, full_source = TRUE)
x <- map(phys, ~ {
    tryCatch(
        error = function(e) e, getDataReadyForPropagation(.x)
    )
})

any(map_lgl(x, rlang::is_error))
