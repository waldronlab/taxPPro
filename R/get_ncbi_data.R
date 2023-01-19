
#' Get NCBI taxids
#'
#' \code{get_ncbi_taxids} retrieves the taxids from taxPPro. These taxids
#' correspond to the sets of the taxonomy statistics in the NCBI.
#'
#' @param keyword Default 'all'. Options b = base, i = informal,
#' l = unclassified, t = uncultured.
#' @param with_taxids Default TRUE.
#'
#' @return A table of taxids with taxonomic information.
#'
#' @export
#'
get_ncbi_taxids <- function(keyword = 'all', with_taxids = TRUE) {
    ## proc all is 'base' + 'informal' + 'unclassified' + 'uncultured'
    ## exclude all is 'base'
    ## exclude informal is 'base' + 'unclassified' + uncultured'
    ## exclude unclassified informal is 'base' + 'uncultured'
    ## exclude unclassified is 'base' + 'uncultured' +  'informal'
    ## exclude unclassified uncultured is 'base' + 'informal'
    ## exclude uncultured informal is 'base' + 'unclassified'
    ## exclude uncultured 'base' + 'unclassified' + 'informal'
    if (keyword == 'all')
        keyword <- 'bilt'
    files <- c(
        bilt = 'proc_all_ids.txt'
    )
    file_name <- paste0('extdata/', files[keyword])
    file_path <- system.file(file_name, package = 'taxPPro')
    taxids_df <- readr::read_table(
        file_path, col_names = 'taxid',
        col_types = readr::cols(taxid = readr::col_character())
    )
    taxids <- taxids_df$taxid

    if (with_taxids) {
        ncbi_taxonomy <- get_ncbi_taxonomy()
        output <- ncbi_taxonomy[ncbi_taxonomy$NCBI_ID %in% taxids, ]
        return(output)
    } else {
        return(taxids)
    }
}
