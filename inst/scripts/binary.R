library(bugphyzz)
library(taxPPro)
library(tidyr)
library(dplyr)
library(purrr)

p <- physiologies()

filtered <- vector('list', length(p))
for (i in seq_along(p)) {
    names(filtered)[i] <- names(p)[i]
    message(names(p)[i])
    res <- filterData(p[[i]])
    if (is.null(res))
        next
    filtered[[i]] <- res
}

ap <- p$`acetate producing`
ap_ready <- getDataReady(filterData(ap))

