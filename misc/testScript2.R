library(bugphyzz)
# library(taxPPro)
library(dplyr)

phys <- physiologies()
phys_filt <- lapply(phys, preSteps) |>
    purrr::discard(is.null) |>
    purrr::discard(~!nrow(.x))

lapply(phys_filt, function(x) table(x$Score))


aer <- physiologies('aerophilicity')[[1]]
ph <- physiologies('optimal ph')[[1]]

aer_upstream <- propagateUpstream(aer, max.tax.level = 'genus')
ph_upstream <- propagateUpstream(ph, max.tax.level = 'genus')
