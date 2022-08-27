## A simple script file to test code
library(bugphyzz)
library(dplyr)
library(purrr)

phys = map(physiologies(), as_tibble)

aer = phys$aerophilicity

aer_conflicts <- get_conflicts(aer)

aer_conflicts$Confidence_in_curation <- factor(
    x = aer_conflicts$Confidence_in_curation,
    levels = c('Low', 'Medium', 'High'),
    # labels = c(1, 2, 3),
    ordered = TRUE
)


aer_conflicts |>
    dplyr::group_by(Taxon_name) |>
    dplyr::slice_max(Confidence_in_curation)



aer |>
    resolve_conflicts()

h <- phys$habitat
x <- resolve_conflicts(h)
x


