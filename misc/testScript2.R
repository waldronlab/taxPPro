library(bugphyzz)
library(taxPPro)
library(dplyr)

aer <- physiologies('aerophilicity')[[1]]
ph <- physiologies('optimal ph')[[1]]

aer_filtered <- preSteps(aer)
ph_filtered <- preSteps(ph)

aer_pscores <- getParentScores(aer_filtered)
ph_pscores <- getParentScores(ph_filtered)

aer_prop <- propagate(aer, prop = 'both')
aer_prop_down <- propagate(aer, prop = 'downstream')

table(aer_prop$Evidence)


## An example with Escherichia

eg_taxon <- 'Escherichia'

aer_filtered |>
    filter(Taxon_name == eg_taxon)

aer_pscores |>
    filter(Taxon_name == eg_taxon) |>
    as.data.frame()

ph_filtered |>
    filter(Taxon_name == eg_taxon)

ph_pscores |>
    filter(Taxon_name == eg_taxon) |>
    as.data.frame()


