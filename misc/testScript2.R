library(bugphyzz)
library(taxPPro)

# aer <- physiologies('aerophilicity')[[1]]
# aer_filtered <- aer |>
#     filterData(tax.id.type = 'Taxon_name') |>
#     freq2Scores() |>
#     resolveAgreements() |>
#     resolveConflicts()
# x <- getConflicts(aer_filtered)
# length(unique(x$Taxon_name))

phys <- physiologies()

output <- vector('list', length(phys))
for (i in seq_along(output)) {
    message('>>>>>>', names(phys)[i])
    output[[i]] <- preSteps(phys[[i]])
    names(output)[i] <- names(phys)[i]
}

is <- phys$`isolation site`
is |>
    filterData(tax.id.type = 'Taxon_name') |>
    freq2Scores() |>
    resolveAgreements()

conflicts <- lapply(phys, getConflicts) |>
    purrr::discard(is.null)
