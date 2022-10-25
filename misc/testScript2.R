library(bugphyzz)
library(taxPPro)
library(dplyr)

aer = physiologies('aerophilicity')[[1]]

aer_plus <- aer |>
    propagate()

count(aer_plus, Evidence)

aer_upstream <- aer |>
    upstream()


x = calcParentScore(strains)
x |>
    filter(.data$Rank == 'species')

aer_downstream <- aer_upstream |>
    downstream()
