
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
#'
# get_ncbi_taxids <- function(keyword = 'all', full_taxonomy = TRUE) {
#     ## bilt - proc all is 'base' + 'informal' + 'unclassified' + 'uncultured'
#     ## b - exclude all is 'base'
#     ## blt - exclude informal is 'base' + 'unclassified' + uncultured'
#     ## bu - exclude unclassified informal is 'base' + 'uncultured'
#     ## bui - exclude unclassified is 'base' + 'uncultured' +  'informal'
#     ## exclude unclassified uncultured is 'base' + 'informal'
#     ## exclude uncultured informal is 'base' + 'unclassified'
#     ## exclude uncultured 'base' + 'unclassified' + 'informal'
#
#     if (keyword == 'all') {
#         keyword <- 'bilt'
#     } else {
#         keyword <- stringr::str_split(keyword, pattern = '') |>
#             unlist() |>
#             sort() |>
#             paste0(collapse = '')
#     }
#
#     files <- .files()
#
#
#     if (!keyword %in% names(files))
#         stop('No such keyword.', call. = FALSE)
#
#     file_name <- paste0('extdata/', files[keyword])
#     file_path <- system.file(file_name, package = 'taxPPro')
#     taxids <- readr::read_table(
#         file_path, col_names = 'taxid',
#         col_types = readr::cols(taxid = readr::col_character())
#     )
#
#     if (full_taxonomy) {
#         ncbi_taxonomy <- get_ncbi_taxonomy()
#         taxids_df <- ncbi_taxonomy[ncbi_taxonomy$NCBI_ID %in% taxids[[1]], ,]
#         return(as.data.frame(taxids_df))
#     } else {
#         return(as.data.frame(taxids))
#     }
# }

# .files <- function() {
#     c(
#         bilt = 'proc_all_ids.txt',
#         b = 'proc_exclude_all.txt',
#         blt = 'proc_exclude_informal.txt',
#         bt = 'proc_exclude_unclassified_informal.txt',
#         bit = 'proc_exclude_unclassified.txt',
#         bi = 'proc_exclude_unclassified_uncultured.txt',
#         bl = 'proc_exclude_uncultured_informal.txt',
#         bil = 'proc_exclude_uncultured.txt'
#     )
#
# }

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
#'
# get_ncbi_taxonomy <- function(force_download = FALSE) {
#
#     superkingdom <- kingdom <- phylum <- class <- order <- family <-
#         genus <- species <- NCBI_ID <- tax_name <- NULL
#
#     if (force_download) {
#         tax_dump <- .ncbi_taxonomy_dump(force = force_download)
#     } else {
#         tax_dump <- .ncbi_taxonomy_dump()
#     }
#
#     message('Extracting files...')
#     ## Untar files
#     temp_dir <- tempdir()
#     nodes_file <- paste0(temp_dir, "/nodes.dmp")
#     rankedlineage_file <- paste0(temp_dir, "/rankedlineage.dmp")
#     utils::untar(
#         tarfile = tax_dump,
#         files = c("rankedlineage.dmp", "nodes.dmp"),
#         exdir = temp_dir
#     )
#
#     delim <- '\t|\t'
#
#     ## Read rankedlineage file
#     rankedlineage_col_names <- c(
#         "NCBI_ID", "Taxon_name", "species", "genus", "family", "order",
#         "class", "phylum", "kingdom", "superkingdom"
#     )
#
#     message('Importing ranked lineage...')
#
#     rankedlineage <- vroom::vroom(
#         rankedlineage_file, col_names = rankedlineage_col_names,
#         show_col_types = FALSE, delim = delim
#     ) %>%
#         dplyr::mutate(
#             superkingdom = stringr::str_remove(superkingdom, '(\\||\t\\|)')
#         ) |>
#         dplyr::relocate(superkingdom, kingdom, phylum, class, order, family,
#                         genus, species, NCBI_ID, Taxon_name) %>%
#         dplyr::mutate(NCBI_ID = as.character(NCBI_ID))
#
#     ## Read nodes file
#     nodes_col_names <- c('NCBI_ID', 'Parent_NCBI_ID', 'Rank',
#                          'embl code', 'ddvision id',
#                          'inherited div flag  (1 or 0)',
#                          'genetic code id',
#                          'inherited GC  flag  (1 or 0)',
#                          'mitochondrial genetic code id',
#                          'inherited MGC flag  (1 or 0)',
#                          'GenBank hidden flag (1 or 0)',
#                          'hidden subtree root flag (1 or 0)',
#                          'comments',
#                          'plastid genetic code id',
#                          'inherited PGC flag  (1 or 0)',
#                          'specified_species',
#                          'hydrogenosome genetic code id',
#                          'inherited HGC flag  (1 or 0)')
#
#     message('Importing nodes...')
#
#     nodes <- vroom::vroom(
#         nodes_file, delim = delim,
#         show_col_types = FALSE, col_names = nodes_col_names,
#         col_types = readr::cols_only(
#             NCBI_ID = readr::col_character(),
#             Parent_NCBI_ID = readr::col_character(),
#             Rank = readr::col_character()
#             )
#         )
#     nodes[[ncol(nodes)]] <- sub('\\|', '', nodes[[ncol(nodes)]])
#
#     message('Combining ranked lineages and nodes...')
#
#     ## Combine into a single taxonomy table
#     dplyr::left_join(rankedlineage, nodes, by = "NCBI_ID") %>%
#         dplyr::filter(
#             Taxon_name %in% c('Archaea', 'Bacteria') |
#             superkingdom %in% c('Archaea', 'Bacteria')
#         ) %>%
#         dplyr::select(-kingdom) |>
#         dplyr::rename(kingdom = superkingdom) |>
#         purrr::discard(~ all(is.na(.x))) |>
#         dplyr::mutate(
#             strain = ifelse(.data$Rank == 'strain', .data$Taxon_name, NA),
#             species = ifelse(.data$Rank == 'species', .data$Taxon_name, .data$species),
#             genus = ifelse(.data$Rank == 'genus', .data$Taxon_name, .data$genus),
#             family = ifelse(.data$Rank == 'family', .data$Taxon_name, .data$family),
#             order = ifelse(.data$Rank == 'order', .data$Taxon_name, .data$order),
#             class = ifelse(.data$Rank == 'class', .data$Taxon_name, .data$class),
#             phylum = ifelse(.data$Rank == 'phylum', .data$Taxon_name, .data$phylum),
#             kingdom = ifelse(.data$Rank == 'superkingdom', .data$Taxon_name, .data$kingdom),
#         ) |>
#         dplyr::relocate(.data$strain, .after = 'species')
# }

## This function downloads the NCBI taxonomy dump file and stores it in the
## package's cache with BiocFileCache
# .ncbi_taxonomy_dump <- function(verbose = TRUE, force = FALSE) {
#
#     bfc <- .get_cache()
#
#     query_search <- BiocFileCache::bfcquery(
#         x = bfc, query = "new_taxdump", field = "rname", exact = TRUE
#     )
#
#     rid <- query_search$rid
#
#     if (isFALSE(!length(rid)) && force)
#         BiocFileCache::bfcremove(bfc, rid)
#
#     if (!length(rid) || force) {
#
#         if (verbose)
#             message('Downloading NCBI taxdump. This might take a while.')
#
#         dir <- tempdir()
#         data_file <- paste0(dir, '/new_taxdump.tar.gz')
#         checksum_file <- paste0(data_file, '.md5')
#
#         data_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz')
#         checksum_url <- paste0('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5')
#
#         utils::download.file(data_url, data_file)
#         utils::download.file(checksum_url, checksum_file)
#
#         actual_checksums <- as.character(tools::md5sum(data_file))
#         expected_checksums <- as.character(utils::read.table(checksum_file)[1,1])
#
#         if (isFALSE(expected_checksums == actual_checksums)) {
#             stop(
#                 'The checksum of the NCBI taxonomy dump file and',
#                 ' the expected md5 aren\' equal. Try again.'
#             )
#         }
#
#         resources <- BiocFileCache::bfcadd(
#             x = bfc, rname = 'new_taxdump', fpath = data_file
#         )
#         rid <- names(resources)
#     }
#
#     BiocFileCache::bfcrpath(bfc, rids = rid)
# }

## Function to crate a cache
.get_cache <- function() {
    cache <- tools::R_user_dir('taxPPro', which = 'cache')
    BiocFileCache::BiocFileCache(cache)
}

#' Estimate ASR
#'
#' \code{extimateASR} performs ASR.
#'
#' @param x A a named character vector. The values are the observed character
#' states and the names are the tip labels. It must be of the same length as
#' the number of tip labels in the phylogenetic tree.
#' @param phy A phylo object. With tip labels and branch lenghts.
#'
#' @return A dataframe.
#'
# estimateASR <- function(x, phy) {
#
#     ## Ancestral character estimation
#     ace_result <- ape::ace(
#         x = x, phy = phy, type = "discrete",
#         method = "ML", model = "ER", CI = TRUE,
#         marginal = TRUE, kappa = 0
#     )
#
#     ## Extract likelihood as dataframe
#     ace_states <- as.data.frame(round(ace_result$lik.anc, 3))
#
#     ace_states
# }


#' Annotate internal nodes
#'
#' \code{annotateInternalNodes} annotates the internal nodes of the output of
#' the `estimateASR` function.
#'
#' @param ace Output of the estimateASR function.
#' @param df A data frame.
#' @param phy A phylo object.
#'
#' @return A data frame
#'
# annotateInternalNodes <- function(ace, df, phy) {
#     limit <- length(phy$tip.label) + 1
#     all_nodes <- dplyr::arrange(df, node)
#     internal_nodes <- all_nodes[limit:nrow(all_nodes),]
#     new_attribute_data <- cbind(internal_nodes, ace)
#     output <- new_attribute_data[!is.na(new_attribute_data$taxid),]
#     output$Evidence <- "asr"
#     output$Frequency <- dplyr::case_when(
#         output$`1` == 1 ~ "always",
#         output$`1` >= .7 & output$`1` < 1 ~ "usually",
#         output$`1` >= .4 & output$`1` < .7 ~ "sometimes",
#         output$`1` >= .1 & output$`1` < .4 ~ "rarely",
#         output$`1` < .1 ~ "never"
#     )
#     output$Taxon_name <- gsub("_", " ", output$sci_name)
#     output$NCBI_ID <- output$taxid
#     output$Rank <- output$rank
#     output$Attribute_source <- TRUE
#     select_cols <- c(
#         "NCBI_ID", "Taxon_name", "Evidence",
#         "Frequency", "Rank"
#     )
#     output[, select_cols]
# }

#' Create entries for ASR
#'
#' @param x Character states.
#' @param phy A phylo object.
#' @param df Tree metadata (a data frame).
#'
#' @return A data frame
#'
# createEntriesASR <- function(x, phy, df) {
#
#     valid_ranks <- c(
#         "superkingdom", "phylum", "class", "order", "family", "genus",
#         "species", "strain"
#     )
#
#     output <- vector("list", length(x))
#     names(output) <- names(x)
#     for (i in seq_along(output)) {
#         output[[i]] <- x[[i]] |>
#             {\(y) estimateASR(x = y, phy = phy)}() |>
#             {\(y) annotateInternalNodes(ace = y, df = df, phy = phy)}()
#     }
#
#     output %>%
#         dplyr::bind_rows(.id = "Attribute") %>%
#         dplyr::filter(
#             Frequency != "never",
#             Rank %in% valid_ranks
#         ) %>%
#         dplyr::relocate(
#             NCBI_ID, Taxon_name, Attribute, Evidence,
#             Frequency, Rank
#         )
# }

# appendASR <- function(df1, df2) {
#     df1_taxids <- df1$NCBI_ID
#     df2 <- df2[!df2$NCBI_ID %in% df1_taxids, ]
#     dplyr::bind_rows(df1, df2)
# }

#' Counts per attribute and rank
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#'
# counts_per_attribute_and_rank <- function(x) {
#     NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
#     x |>
#         dplyr::filter(
#             !is.na(NCBI_ID), NCBI_ID != 'unknown',
#             !is.na(Taxon_name),
#             Attribute != '', Attribute_value != FALSE,
#             !is.na(Rank), Rank != '',
#             Rank %in% c('genus', 'species', 'strain')
#             ) |>
#         dplyr::count(Attribute, Rank)
# }


#' Counts per  rank
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#'
# counts_per_rank <- function(x) {
#     NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
#     x |>
#         dplyr::filter(
#             !is.na(NCBI_ID), NCBI_ID != 'unknown',
#             !is.na(Taxon_name),
#             Attribute != '', Attribute_value != FALSE,
#             !is.na(Rank), Rank != '',
#             Rank %in% c('genus', 'species', 'strain')
#         ) |>
#         dplyr::count(Rank)
# }

#' Total counts
#'
#' @param x  A data frame.
#'
#' @return A data frame.
#'
# counts_total <- function(x) {
#     NCBI_ID <- Taxon_name <- Attribute <- Attribute_value <- Rank <- NULL
#     x |>
#         dplyr::filter(
#             !is.na(NCBI_ID), NCBI_ID != 'unknown',
#             !is.na(Taxon_name),
#             Attribute != '', Attribute_value != FALSE,
#             !is.na(Rank), Rank != '',
#             Rank %in% c('genus', 'species', 'strain')
#         ) |>
#         dplyr::count()
# }

#' Print data.tree attributes
#'
#' \code{printDataTreeAttributes} prints all of the attributes in a data.tree.
#'
#' @param data_tree A data.tree object with attributes.
#' @param limit The number of nodes to be displayed. Default = 100.
#'
#' @return A data.frame or data on the console. Can use \code{View}.
#'
# printDataTreeAttributes <- function(data_tree, limit = 100) {
#   attrs <- as.list(data_tree$attributesAll)
#   lim = list(limit = limit)
#   args <- c(list(data_tree), attrs)
#   args <- c(args, lim)
#   do.call('print', args = args)
# }

#' toDataFrame
#'
#' \code{toDataFrame}
#'
#' @param data_tree A data.tree.
#' @param ncbi_tax NCBI taxonomy. Output of \code{get_ncbi_taxonomy}.
#'
#' @return A data.frame.
#'
# toDataFrame <- function(data_tree, ncbi_tax) {
#   args <- as.list(data_tree$attributesAll)
#   args <- c(list(x = data_tree, row.names = NULL, optional = FALSE), args)
#   df <- do.call('as.data.frame', args)
#   df$levelName <- stringr::str_squish(sub('.*-', '', df$levelName))
#   df <- df[df$levelName != 'ArcBac',]
#   df <- tidyr::separate(
#     df, col = 'levelName', into = c('Rank', 'NCBI_ID'), sep = '__'
#   )
#   dict <- c(
#     d = 'domain', p = 'phylum', c = 'class', o = 'order', f = 'family',
#     g = 'genus', s = 'species', t = 'strain'
#   )
#   df$Rank <- dict[df$Rank]
#   # ncbi_tax <- get_ncbi_taxonomy()
#   if ('Rank' %in% colnames(ncbi_tax)) {
#     ncbi_tax$Rank <- NULL
#   }
#   pos <- which(colnames(ncbi_tax) == 'kingdom')
#   names(ncbi_tax)[pos] <- 'domain'
#   output <- dplyr::left_join(df, ncbi_tax, by = 'NCBI_ID')
#   return(output)
# }



#' Prepare Data for Propagation
#'
#' \code{prepareDataForPropagation} prepares data to be used in the
#' propagation method.
#'
#' @param df A data.frame.
#' @param resolve Logical value. If TRUE (default) resolve agreements,
#' conflicts, and double annotations.
#'
#' @return A data.frame
#'
# prepareDatForPropagation <- function(df, resolve = TRUE) {
#   df$NCBI_ID[which(is.na(df$NCBI_ID))] <- 'unknown'
#   df$Parent_NCBI_ID[which(is.na(df$Parent_NCBI_ID))] <- 'unknown'
#   df <- df[df$Parent_NCBI_ID != 'unknown',]
#   df <- df[which(!is.na(df$Rank)), ]
#   df <- df[which(!is.na(df$Evidence)), ]
#   df <- df[which(!is.na(df$Frequency)), ]
#   df <- df[which(!is.na(df$Confidence_in_curation)), ]
#   df <- df |>
#     dplyr::group_by(.data$Parent_NCBI_ID) |>
#     dplyr::mutate(
#       Parent_name = paste(unique(.data$Parent_name), collapse = ';')
#     ) |>
#     dplyr::ungroup() |>
#     dplyr::distinct()
#   df <- unique(df)
#   df$Score <- freq2Scores(df$Frequency)
#   df_yesid <- df[which(df$NCBI_ID != 'unknown'),]
#   if (nrow(df_yesid) > 0) {
#     df_yesid <- df_yesid |>
#       dplyr::group_by(.data$NCBI_ID) |>
#       dplyr::mutate(
#         Taxon_name = paste(unique(.data$Taxon_name), collapse = ';')
#       ) |>
#       dplyr::ungroup() |>
#       dplyr::distinct()
#     # df_yesid <- df_yesid[grep(';', df_yesid$Taxon_name),]
#   }
#   df_noid <- df[which(df$NCBI_ID == 'unknown'),]
#   if (nrow(df_noid) > 0) {
#     df_noid_asr <- calcParentScores(df_noid)
#     df_new <- dplyr::bind_rows(df_yesid, df_noid_asr)
#   } else {
#     df_new <- df_yesid
#   }
#   dict <- c(genus = 'g__', species = 's__', strain = 't__')
#   df_new$NCBI_ID <- paste0(dict[df_new$Rank], df_new$NCBI_ID)
#   df_new <- df_new[which(!startsWith(colnames(df_new), 'Parent'))]
#
#   attr_type <- unique(df$Attribute_type)
#   cols <- c(
#     'NCBI_ID', 'Taxon_name', 'Rank',
#     'Attribute', 'Attribute_source',
#     'Evidence', 'Frequency',
#     'Attribute_type', 'Attribute_group',
#     'Confidence_in_curation', 'Score'
#   )
#   if (attr_type == 'logical') {
#     cols <- c(cols, 'Attribute_value')
#   } else if (attr_type == 'range') {
#     cols <- c(cols, c('Attribute_value_min', 'Attribute_value_max'))
#   }
#   df_new <- df_new[,cols]
#   df_new <- unique(df_new)
#
#   if (resolve) {
#     output <- df_new |>
#       resolveAgreements() |>
#       resolveConflicts() |>
#       dplyr::distinct() |>
#       as.data.frame()
#     return(output)
#   } else {
#     return(df_new)
#   }
# }
#
#
