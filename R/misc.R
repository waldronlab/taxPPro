get_silva_table <- function(x) {
    db_url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_12_2021.csv'
    db <- vroom::vroom(db_url, delim = '\t', col_names = NULL)
}

get_silve_tree <- function(x) {

}

tree_metaphlan <- function() {
    tree_url <- 'https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk'
    tree <- treeio::read.tree(tree_url)
    tree$tip.label <- gsub('_', ' ', sub('^.+s__', '', tree$tip.label))
    tree

}


