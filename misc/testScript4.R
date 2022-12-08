
library(bugphyzz)
library(purrr)
library(taxPPro)
library(dplyr)
library(purrr)

phys <- physiologies()
phys <- map(phys, as_tibble)
remove_dat <- c(
    'habitat', 'antimicrobial resistance', 'isolation site', 'width',
    'length', 'disease association',
    'animal pathogen', 'biofilm forming', 'growth temperature'
)

phys <- phys[!names(phys) %in% remove_dat]

phys_prop <- map(phys, ~ propagate(.x, 'genus', 'both'))

