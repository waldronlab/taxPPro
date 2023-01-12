
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

    if (force_download) {
        tax_dump <- .ncbi_taxonomy_dump(force = force_download)
    } else {
        tax_dump <- .ncbi_taxonomy_dump()
    }

    message('Extracting files...')
    ## Untar files
    temp_dir <- tempdir()
    nodes_file <- paste0(temp_dir, "/nodes.dmp")
    rankedlineage_file <- paste0(temp_dir, "/rankedlineage.dmp")
    utils::untar(
        tarfile = tax_dump,
        files = c("rankedlineage.dmp", "nodes.dmp"),
        exdir = temp_dir
    )

    delim <- '\t|\t'

    ## Read rankedlineage file
    rankedlineage_col_names <- c(
        "NCBI_ID", "Taxon_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"
    )

    message('Importing ranked lineage...')

    rankedlineage <- vroom::vroom(
        rankedlineage_file, col_names = rankedlineage_col_names,
        show_col_types = FALSE, delim = delim
    ) %>%
        dplyr::mutate(
            superkingdom = stringr::str_remove(superkingdom, '(\\||\t\\|)')
        ) |>
        dplyr::relocate(superkingdom, kingdom, phylum, class, order, family,
                        genus, species, NCBI_ID, Taxon_name) %>%
        dplyr::mutate(NCBI_ID = as.character(NCBI_ID))

    ## Read nodes file
    nodes_col_names <- c('NCBI_ID', 'Parent_NCBI_ID', 'Rank',
                         'embl code', 'ddvision id',
                         'inherited div flag  (1 or 0)',
                         'genetic code id',
                         'inherited GC  flag  (1 or 0)',
                         'mitochondrial genetic code id',
                         'inherited MGC flag  (1 or 0)',
                         'GenBank hidden flag (1 or 0)',
                         'hidden subtree root flag (1 or 0)',
                         'comments',
                         'plastid genetic code id',
                         'inherited PGC flag  (1 or 0)',
                         'specified_species',
                         'hydrogenosome genetic code id',
                         'inherited HGC flag  (1 or 0)')

    message('Importing nodes...')

    nodes <- vroom::vroom(
        nodes_file, delim = delim,
        show_col_types = FALSE, col_names = nodes_col_names,
        col_types = readr::cols_only(
            NCBI_ID = readr::col_character(),
            Parent_NCBI_ID = readr::col_character(),
            Rank = readr::col_character()
            )
        )
    nodes[[ncol(nodes)]] <- sub('\\|', '', nodes[[ncol(nodes)]])

    message('Combining ranked lineages and nodes...')

    ## Combine into a single taxonomy table
    dplyr::left_join(rankedlineage, nodes, by = "NCBI_ID") %>%
        dplyr::filter(
            Taxon_name %in% c('Archaea', 'Bacteria') |
            superkingdom %in% c('Archaea', 'Bacteria')
        ) %>%
        dplyr::select(-kingdom) |>
        dplyr::rename(kingdom = superkingdom) |>
        purrr::discard(~ all(is.na(.x))) |>
        dplyr::mutate(
            strain = ifelse(.data$Rank == 'strain', .data$Taxon_name, NA),
            species = ifelse(.data$Rank == 'species', .data$Taxon_name, .data$species),
            genus = ifelse(.data$Rank == 'genus', .data$Taxon_name, .data$genus),
            family = ifelse(.data$Rank == 'family', .data$Taxon_name, .data$family),
            order = ifelse(.data$Rank == 'order', .data$Taxon_name, .data$order),
            class = ifelse(.data$Rank == 'class', .data$Taxon_name, .data$class),
            phylum = ifelse(.data$Rank == 'phylum', .data$Taxon_name, .data$phylum),
            kingdom = ifelse(.data$Rank == 'superkingdom', .data$Taxon_name, .data$kingdom),
        ) |>
        dplyr::relocate(.data$strain, .after = 'species')
}

## This function downloads the NCBI taxonomy dump file and stores it in the
## package's cache with BiocFileCache
.ncbi_taxonomy_dump <- function(verbose = TRUE, force = FALSE) {

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
