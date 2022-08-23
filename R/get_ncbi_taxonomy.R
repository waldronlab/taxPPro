
#' get_ncbi_taxonomy
#'
#' \code{get_ncbi_taxonomy} downloads and imports the full taxonomic
#' classification from the NCBI database
#'
#' @param force_download If FALSE (default) the function uses the ncbi
#' taxonomy in cache (if present). If TRUE the ncbi taxonomy is dowloaded.
#' If a taxonomy was present in the cache, it will be removed and replaced.
#'
#' @importFrom magrittr %>%
#'
#' @return A data frame with the complete NCBI taxonomy
#' @export
#'
get_ncbi_taxonomy <- function(force_download = FALSE) {

    superkingdom <- kingdom <- phylum <- class <- order <- family <-
        genus <- species <- NCBI_ID <- tax_name <- NULL

    ## Untar files
    temp_dir <- tempdir()
    nodes_file <- paste0(temp_dir, "/nodes.dmp")
    rankedlineage_file <- paste0(temp_dir, "/rankedlineage.dmp")
    utils::untar(
        tarfile = .ncbi_taxonomy_dump(force = force_download),
        files = c("rankedlineage.dmp", "nodes.dmp"),
        exdir = temp_dir
    )

    ## Read rankedlineage file
    rankedlineage_col_names <- c(
        "NCBI_ID", "tax_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"
    )

    rankedlineage <- vroom::vroom(
        rankedlineage_file, col_names = FALSE, show_col_types = FALSE
    ) %>%
        purrr::discard(
            ~all(is.na(.x) | all(stringr::str_detect(.x, "\\|")))
        ) %>%
        magrittr::set_colnames(rankedlineage_col_names) %>%
        dplyr::relocate(superkingdom, kingdom, phylum, class, order, family,
                        genus, species, NCBI_ID, tax_name) %>%
        dplyr::mutate(NCBI_ID = as.character(NCBI_ID))

    ## Read nodes file
    nodes <- vroom::vroom(
        nodes_file, delim = "|", show_col_types = FALSE, col_names = FALSE,
        col_types = readr::cols_only(
            X1 = readr::col_character(), X3 = readr::col_character()
        )
    ) %>%
        purrr::map_df(~ stringr::str_remove_all(.x, "\t")) %>%
        magrittr::set_colnames(c("NCBI_ID", "rank")) %>%
        dplyr::mutate(NCBI_ID = as.character(NCBI_ID))

    ## Combine into taxonomy_table
    dplyr::left_join(rankedlineage, nodes, by = "NCBI_ID") %>%
        dplyr::filter(superkingdom %in% c('Archaea', 'Bacteria')) %>%
        purrr::discard(~ all(is.na(.x))) %>%
        dplyr::select(-kingdom) %>%
        dplyr::rename(kingdom = superkingdom)

}


## This function downloads the NCBI taxonomy dump file and stores it in the
## package's cache with BiocFileCache
.ncbi_taxonomy_dump <- function(verbose = FALSE, force = FALSE) {

    bfc <- .get_cache()

    query_search <- BiocFileCache::bfcquery(
        x = bfc, query = "new_taxdump", field = "rname", exact = TRUE
    )

    rid <- query_search$rid

    if (isFALSE(!length(rid)) && force)
        BiocFileCache::bfcremove(bfc, rid)

    if (!length(rid) || force) {

        if (verbose)
            message('Downloading NCBI taxdump. This might take a while.')

        dir <- tempdir()
        data_file <- paste0(dir, '/new_taxdump.tar.gz')
        checksum_file <- paste0(data_file, '.md5')

        data_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz')
        checksum_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5')

        utils::download.file(data_url, data_file)
        utils::download.file(checksum_url, checksum_file)

        actual_checksums <- as.character(tools::md5sum(data_file))
        expected_checksums <- as.character(utils::read.table(checksum_file)[1,1])

        if (isFALSE(expected_checksums == actual_checksums)) {
            stop(
                'The checksum of the NCBI taxonomy dump file and',
                ' the expected md5 aren\' equal. Try again.'
            )
        }

        resources <- BiocFileCache::bfcadd(
            x = bfc, rname = 'new_taxdump', fpath = data_file
        )
        rid <- names(resources)
    }

    BiocFileCache::bfcrpath(bfc, rids = rid)
}

## Function to crate a cache
.get_cache <- function() {
    cache <- tools::R_user_dir('taxPPro', which = 'cache')
    BiocFileCache::BiocFileCache(cache)
}
