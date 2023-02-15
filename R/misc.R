get_silva_table <- function(x) {
    db_url <- 'https://imedea.uib-csic.es/mmg/ltp/wp-content/uploads/ltp/LTP_12_2021.csv'
    db <- vroom::vroom(db_url, delim = '\t', col_names = NULL)
}

#' Get metaplhan tree
#'
#' \code{getMpaTree} gets metplhan tree. Tips are taxids. v4.0.0.
#'
#' @return A phylo object
#' @export
#'
getMpaTree <- function() {
    file_url <- 'https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk'
    treeio::read.tree(file_url)
    # tree_url <- 'https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk'
    # tree$tip.label <- gsub('_', ' ', sub('^.+s__', '', tree$tip.label))
}

#' Get metaphlan table
#'
#' \code{getMpaTable} gets metaphlan table v4.0.0
#'
#' @return A data.frame
#' @export
#'
getMpaTable <- function() {
    ranks <- c(
        'domain', 'phylym', 'class', 'order', 'family', 'genus', 'species'
    )
    file_url <- 'https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB2GTDB.tsv'
    df <- utils::read.table(file_url, header = FALSE, sep = '\t')
    colnames(df) <- c('id', 'taxonomy')
    df <- tidyr::separate(
        data = df, col = 'taxonomy', into = ranks, sep = ';'
    )
    purrr::modify_at(.x = df, .at = ranks, .f = ~ sub('[dpcofgs]__', '', .x))
}
