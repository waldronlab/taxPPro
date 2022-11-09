library(bugphyzz)

aer <- physiologies('aerophilicity')[[1]]

aer_filtered_taxids <- filterData(aer, tax.id.type = 'NCBI_ID')
aer_filtered_taxnames <- filterData(aer, tax.id.type = 'Taxon_name')
