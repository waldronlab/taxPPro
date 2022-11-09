library(bugphyzz)

aer <- physiologies('aerophilicity')[[1]]

aer_filtered_taxids <- filterData(aer, tax.id.type = 'NCBI_ID')
aer_filtered_taxnames <- filterData(aer, tax.id.type = 'Taxon_name')

aer_scores_taxids <- freq2Scores(aer_filtered_taxids)
aer_scores_taxnames <- freq2Scores(aer_filtered_taxnames)


x = getDuplicates(aer_filtered_taxids)
