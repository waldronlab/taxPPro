## This script is for the living tree project
##
## This is the url of the file (active as of Aug 30, 2023)
## I also backed up the file in this package (to make things faster).
## url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_all_06_2022.ntree'

library(ape)
library(ggtree)
library(taxonomizr)
library(taxPPro)
library(data.tree)
library(bugphyzz)
library(dplyr)

data('tree_list')
data_tree <- as.Node(tree_list)

nodes <- data_tree$Get(function(node) {
    if (grepl('s__', node$name)) return(node$name)
})
nodes <- unname(nodes)
nodes <- sub('^\\w__', '', nodes)
nodes <- nodes[-1]
nodes <- nodes[!is.na(nodes)]

t_file <- system.file(
    'extdata', 'LTP_all_06_2022.ntree', package = 'taxPPro', mustWork = TRUE
)
t <- read.tree(t_file)
p <- ggtree(t, layout = 'circular', size = 0.025)

tip_labels <- t$tip.label
node_labesl <- t$node.label

sql <- '~/accessionTaxa.sql' # this takes a while - date Aug 29, 2023
acc <- sub("^'([^,]+).*", "\\1", tip_labels)
taxids <- accessionToTaxa(accessions = acc, sqlFile = sql, version = 'base')
names(taxids) <- acc
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
