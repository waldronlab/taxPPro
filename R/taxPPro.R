filterData <- function(tbl) {
    attr_type <- unique(tbl$Attribute_type)
    if(attr_type == 'logical'){
        output <- filterDataMulti(tbl)
    }
    return(output)
}

filterDataMulti <- function(tbl) {
    select_cols <- c(
        'NCBI_ID', 'Taxon_name', 'Parent_NCBI_ID',
        'Attribute','Attribute_source', 'Confidence_in_curation',
        'Frequency', 'Score', 'Evidence',
        'Attribute_type', 'Attribute_group'
    )
    # valid_ranks <- c('genus', 'species', 'strain')
    attributes_fname <- system.file(
        'extdata', 'attributes.tsv', package = 'bugphyzz'
    )
    attributes <- utils::read.table(attributes_fname, sep = '\t', header = TRUE)
    valid_attributes <- attributes |>
        dplyr::filter(.data$attribute_group == phys_name) |>
        dplyr::pull(.data$attribute) |>
        unique()
    phys_data <- tbl |>
        tibble::as_tibble() |>
        dplyr::filter(.data$Attribute_value == TRUE) |>
        dplyr::filter(.data$Attribute %in% valid_attributes) |>
        dplyr::filter(
            !((is.na(.data$NCBI_ID) | .data$NCBI_ID == 'unknown') & is.na(.data$Parent_NCBI_ID))
        ) |>
        dplyr::filter(!is.na(.data$Attribute_source), !is.na(.data$Frequency)) |>
        dplyr::mutate(Score = freq2Scores(.data$Frequency)) |>
        dplyr::select(tidyselect::all_of(select_cols)) |>
        dplyr::distinct()
    n_dropped_rows <- nrow(tbl) - nrow(phys_data)
    message(format(n_dropped_rows, big.mark = ','), ' rows were dropped.')
    return(phys_data)
}

getDataReady <- function(tbl) {
    set_With_ids <- getSetWithIDs(tbl)
    set_without_ids <- getSetWithoutIDs(tbl)
    dplyr::bind_rows(set_with_ids, set_without_ids) |>
        tidyr::complete(NCBI_ID, Attribute, fill = list(Score = 0)) |>
        dplyr::arrange(NCBI_ID, Attribute)
}


getSetWithIDs <- function(tbl) {
    valid_ranks <- c('genus', 'species', 'strain')
    lgl_vct <- is.na(tbl$NCBI_ID) | tbl$NCBI_ID == 'unknown'
    tbl |>
        dplyr::filter(!lgl_vct) |>
        dplyr::mutate(
            Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')
        ) |>
        dplyr::filter(.data$Rank %in% valid_ranks) |>
        dplyr::mutate(
            Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')
        ) |>
        dplyr::distinct() |>
        dplyr::mutate(
            Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)
        ) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
            .data$Confidence_in_curation, n = 1, with_ties = TRUE
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
        dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::mutate(
            Total_score = sum(.data$Score),
            Score = .data$Score / .data$Total_score
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(Frequency = scores2Freq(.data$Score)) |>
        dplyr::select(-.data$Parent_NCBI_ID, -.data$Total_score) |>
        dplyr::mutate(taxid = .data$NCBI_ID) |>
        dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
        dplyr::filter(!is.na(.data$NCBI_ID)) |>
        dplyr::distinct() |>
        dplyr::arrange(.data$NCBI_ID, .data$Attribute) |>
        dplyr::relocate(all_of(.orderedColumns()))
}

getSetWithoutIDs <- function(tbl, set_with_ids = NULL) {
    attribute_type_var <- unique(tbl$Attribute_type)
    attribute_group_var <- unique(tbl$Attribute_group)

    if (is.null(set_with_ids))
        set_with_ids <- data.frame(NCBI_ID = 'not a real NCBI_ID')

    valid_ranks <- c('genus', 'species', 'strain')
    lgl_vct <- is.na(tbl$NCBI_ID) | tbl$NCBI_ID == 'unknown'
    tbl |>
        dplyr::filter(lgl_vct) |>
        dplyr::select(
            -.data$NCBI_ID, -.data$Taxon_name, -.data$Frequency
        ) |>
        dplyr::relocate(NCBI_ID = .data$Parent_NCBI_ID) |>
        dplyr::mutate(Rank = taxizedb::taxid2rank(.data$NCBI_ID, db = 'ncbi')) |>
        dplyr::filter(Rank %in% valid_ranks) |>
        dplyr::mutate(Taxon_name = taxizedb::taxid2name(.data$NCBI_ID, db = 'ncbi')) |>
        dplyr::distinct() |>
        dplyr::mutate(Confidence_in_curation =  conf2Fct(.data$Confidence_in_curation)) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(.data$Confidence_in_curation, n = 1, with_ties = TRUE) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID, .data$Attribute) |>
        dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::mutate(
            Total_score = sum(.data$Score),
            Score = .data$Score / .data$Total_score
        ) |>
        dplyr::mutate(Frequency = scores2Freq(.data$Score)) |>
        dplyr::mutate(
            Evidence = 'tax',
            Attribute_group = Attribute_group_var,
            Attribute_type = Attribute_type_var,
            Attribute_source = NA,
            Confidence_in_curation = NA
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.data$Total_score) |>
        dplyr::mutate(taxid = .data$NCBI_ID) |>
        dplyr::mutate(NCBI_ID = addRankPrefix(.data$NCBI_ID, .data$Rank)) |>
        dplyr::filter(!is.na(.data$NCBI_ID)) |>
        dplyr::filter(!.data$NCBI_ID %in% unique(set_with_ids$NCBI_ID)) |>
        dplyr::distinct() |>
        dplyr::arrange(.data$NCBI_ID, .data$Attribute) |>
        dplyr::relocate(tidyselect::all_of(.orderedColumns()))
}

.orderedColumns <- function() {
    c(
        'NCBI_ID', 'Taxon_name', 'Rank',
        'Attribute', 'Attribute_source', 'Confidence_in_curation',
        'Evidence', 'Frequency', 'Score',
        'Attribute_group', 'Attribute_type',
        'taxid'
    )
}





