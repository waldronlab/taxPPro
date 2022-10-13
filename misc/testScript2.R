library(bugphyzz)
library(taxPPro)
aer = physiologies('aerophilicity')[[1]]

aer_plus <- aer |>
    propagate()

