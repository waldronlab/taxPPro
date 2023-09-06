## This script is for the living tree project
##
## This is the url of the file (active as of Aug 30, 2023)
## I also backed up the file in this package (to make things faster).
## url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'

library(ape)
library(taxonomizr)

tree <- read.tree('https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree')
tip_labels <- t$tip.label
node_labesl <- t$node.label
sql <- '~/accessionTaxa.sql' # this takes a while - date Aug 29, 2023
acc <- sub("^'([^,]+).*", "\\1", tip_labels)
taxids <- accessionToTaxa(accessions = acc, sqlFile = sql, version = 'base')
names(taxids) <- acc

tree_data <- data.frame(
    tip_label = tip_labels,
    accesion = acc,
    taxid = unname(taxids)

)

fix(tree_data)







names(taxids)[which(is.na(taxids))]

x <- taxids[which(!is.na(taxids))]
ranks <- taxizedb::taxid2rank(x, db = 'ncbi')

t$tip.label <- taxids


bp <- importBugphyzz()
bp2 <- bp |>
    filter(
        !Evidence %in% c('asr', 'inh'),
        !Rank %in% c('genus', 'strain')
    )
bp_ids <- unique(bp2$NCBI_ID)
sum(bp_ids %in% nodes)
sum(bp_ids %in% taxids)

mean(taxids %in% nodes)
mean(nodes %in% taxids)
