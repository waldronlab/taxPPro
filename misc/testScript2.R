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
b <- getConflicts(dup_filtered)
c <- getAgreements(dup_filtered)

a2 <- dup |>
    dplyr::count(.data$Attribute_source, .data$Taxon_name) |>
    dplyr::filter(.data$n > 1)

danno <- dup |>
    filter(Taxon_name %in% a2$Taxon_name) |>
    arrange(NCBI_ID, Attribute)

danno2 <- getDoubleAnnotations(dup)


## Resolve agreements, then resolve conflicts. Double annotations don't need
## to be resolved. They can remain as they are in the annotations.


