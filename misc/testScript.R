## A simple script file to test code
library(bugphyzz)
library(dplyr)
library(purrr)

phys <- map(physiologies(), as_tibble)

get_dups <- function(df) {
    double_annotations <- get_double_annotations(df)
    agreements <- get_agreements(df)
    conflicts <- get_conflicts(df)
    all_dups <- list(
        double_annotation = double_annotations,
        agreement = agreements,
        conflict = conflicts
    ) |>
        dplyr::bind_rows(.id = 'Duplicate_type')
    if (!nrow(all_dups)) {
        return(NULL)
    } else {
        return(all_dups)
    }
}

dups_1 <- map(phys, get_duplicates) |>
    discard(is.null)

dups <- map(phys, get_dups) |>
    discard(is.null)


sum(!map_int(dups_1, nrow) == map_int(dups, nrow))





