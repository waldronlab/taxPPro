library(bugphyzz)
library(dplyr)

# phys <- physiologies()
# duplicates <- phys |>
#     lapply(function(x) filterData(x, tax.id.type = 'NCBI_ID')) |>
#     purrr::discard(is.null) |>
#     lapply(freq2Scores) |>
#     lapply(getDuplicates) |>
#     purrr::discard(is.null)

aer <- physiologies('aerophilicity')[[1]]
dup_filtered <- filterData(aer, tax.id.type = 'Taxon_name') |>
    freq2Scores() |>
    getDuplicates()

## From duplicates we should get
a <- getDoubleAnnotations(dup_filtered)
b <- getAgreements(dup_filtered)
c <- getConflicts(dup_filtered)

resolved_conflicts <- resolveConflicts(dup_filtered)
