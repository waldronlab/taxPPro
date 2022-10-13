library(bugphyzz)
aer = physiologies('aerophilicity')[[1]]

aer_plus <- aer |>
    propagate()

