
#' Get NCBI taxids
#'
#' \code{get_ncbi_taxids} retrieves the taxids from taxPPro. These taxids
#' correspond to the sets of the taxonomy statistics in the NCBI.
#'
#' @param keyword Default 'all'. Options b = base, i = informal,
#' l = unclassified, t = uncultured.
#' @param full_taxonomy Default TRUE.
#'
#' @return A table of taxids with taxonomic information.
#'
#' @export
#'
get_ncbi_taxids <- function(keyword = 'all', full_taxonomy = TRUE) {
    ## bilt - proc all is 'base' + 'informal' + 'unclassified' + 'uncultured'
    ## b - exclude all is 'base'
    ## blt - exclude informal is 'base' + 'unclassified' + uncultured'
    ## bu - exclude unclassified informal is 'base' + 'uncultured'
    ## bui - exclude unclassified is 'base' + 'uncultured' +  'informal'
    ## exclude unclassified uncultured is 'base' + 'informal'
    ## exclude uncultured informal is 'base' + 'unclassified'
    ## exclude uncultured 'base' + 'unclassified' + 'informal'

    if (keyword == 'all') {
        keyword <- 'bilt'
    } else {
        keyword <- stringr::str_split(keyword, pattern = '') |>
            unlist() |>
            sort() |>
            paste0(collapse = '')
    }

    files <- .files()


    if (!keyword %in% names(files))
        stop('No such keyword.', call. = FALSE)

    file_name <- paste0('extdata/', files[keyword])
    file_path <- system.file(file_name, package = 'taxPPro')
    taxids <- readr::read_table(
        file_path, col_names = 'taxid',
        col_types = readr::cols(taxid = readr::col_character())
    )

    if (full_taxonomy) {
        ncbi_taxonomy <- get_ncbi_taxonomy()
        taxids_df <- ncbi_taxonomy[ncbi_taxonomy$NCBI_ID %in% taxids[[1]], ,]
        return(as.data.frame(taxids_df))
    } else {
        return(as.data.frame(taxids))
    }
}

.files <- function() {
    c(
        bilt = 'proc_all_ids.txt',
        b = 'proc_exclude_all.txt',
        blt = 'proc_exclude_informal.txt',
        bt = 'proc_exclude_unclassified_informal.txt',
        bit = 'proc_exclude_unclassified.txt',
        bi = 'proc_exclude_unclassified_uncultured.txt',
        bl = 'proc_exclude_uncultured_informal.txt',
        bil = 'proc_exclude_uncultured.txt'
    )

}
