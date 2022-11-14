library(bugphyzz)
library(taxPPro)

eg_taxon <- 'Lactobacillus sakei'
aer <- physiologies('aerophilicity')[[1]]
aer[aer$Taxon_name == eg_taxon, ]

aer_filtered <- aer |>
    filterData(tax.id.type = 'Taxon_name') |>
    freq2Scores() |>
    resolveAgreements() |>
    resolveConflicts()

a <- getDoubleAnnotations(aer_filtered)
b <- getAgreements(aer_filtered)
c <- getConflicts(aer_filtered)

aer_filtered[aer_filtered$Taxon_name == eg_taxon,]

