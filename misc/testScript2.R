
suppressMessages({
    library(bugphyzz)
    library(dplyr)
})

ar <- as_tibble(physiologies('antimicrobial resistance')[[1]])
dup_rows <- which(duplicated(ar[,c('NCBI_ID', 'Taxon_name', 'Attribute')]))
length(dup_rows)

patric_ids <- c(
    '1733.103', '1733.223', '1733.105'
)

ar |>
    filter(PATRIC_ID %in% patric_ids) |>
    mutate(PATRIC_ID = as.character(PATRIC_ID)) |>
    select(NCBI_ID, Taxon_name, PATRIC_ID) |>
    distinct()




x <- ar |>
    propagate(max.tax.level = 'genus', direction = 'upstream')


filtered <- preSteps(ar, 'Taxon_name')



View(ar[dup_rows, ])

ar |>
    filter(
        Taxon_name == 'Mycobacterium tuberculosis G04046',
        Attribute == 'resistance to streptomycin'
    ) |>
    View()

