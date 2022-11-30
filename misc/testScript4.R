
library(bugphyzz)
library(purrr)
library(taxPPro)
library(dplyr)

phys <- physiologies()
phys <- map(phys, as_tibble)
remove_dat <- c(
    'habitat', 'antimicrobial resistance', 'isolation site', 'width',
    'length', 'disease association',
    'animal pathogen', 'biofilm forming', 'growth temperature'
)

# phys <- phys[!names(phys) %in% remove_dat]

phys_char <- phys |>
    keep(~ is.character(.x$Attribute_value))

phys_char |>
    map(head) |>
    bind_rows(.id = 'dataset') |>
    View()

