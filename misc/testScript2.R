library(bugphyzz)
library(dplyr)
library(purrr)

# phys <- physiologies()
# duplicates <- phys |>
#     lapply(function(x) filterData(x, tax.id.type = 'NCBI_ID')) |>
#     purrr::discard(is.null) |>
#     lapply(freq2Scores) |>
#     lapply(getDuplicates) |>
#     purrr::discard(is.null)

aer <- physiologies('aerophilicity')[[1]]
dup_filtered <- aer |>
    filterData(tax.id.type = 'Taxon_name') |>
    freq2Scores() |>
    resolveAgreements() |>
    resolveConflicts()

## From duplicates we should get
b <- getAgreements(dup_filtered)
c <- getConflicts(dup_filtered)

















