
#' Import NCBI IDs from the statistics page
#'
#' \code{imrpot_ncbi_ids} imports the NCBI ids as grouped in the statistics
#' page. This data was imported by manually selecting exclusion options in
#' the NCBI taxonomy statistics page: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&taxon=2157%3AArchaea&taxon=2%3ABacteria&period=&from=&to=
#'
#' The data was downloaded on
#' September 2, 2022. Current numbers might be different.
#'
#' @return A list of NCBI taxids.
#' @export
#'
import_ncbi_ids <- function() {
    ## 1. Include all ids.
    all_ids_fname <- system.file(
        'extdata/proc_all_ids.txt', package = 'taxPPro'
    )
    all_ids <-
        read.table(all_ids_fname, header = FALSE)[[1]] |>
        as.character()

    ## 2. Exclude all
    exclude_all_fname <- system.file(
        'extdata/proc_exclude_all.txt', package = 'taxPPro'
    )
    exclude_all <-
        read.table(exclude_all_fname, header = FALSE)[[1]] |>
        as.character()

    ## 3. Exclude informal names
    exclude_informal_fname <- system.file(
        'extdata/proc_exclude_informal.txt', package = 'taxPPro'
    )
    exclude_informal <-
        read.table(exclude_informal_fname, header = FALSE)[[1]] |>
        as.character()

    ## 4. Exclude unclassified
    exclude_unclassified_fname <- system.file(
        'extdata/proc_exclude_unclassified.txt', package = 'taxPPro'
    )
    exclude_unclassified <-
        read.table(exclude_unclassified_fname, header = FALSE)[[1]] |>
        as.character()

    ## 5. Exclude uncultured
    exclude_uncultured_fname <- system.file(
        'extdata/proc_exclude_uncultured.txt', package = 'taxPPro'
    )
    exclude_uncultured <-
        read.table(exclude_uncultured_fname, header = FALSE)[[1]] |>
        as.character()

    ## 6. Exclude unclassified and uncultured
    exclude_unclassified_uncultured_fname <- system.file(
        'extdata/proc_exclude_unclassified_uncultured.txt', package = 'taxPPro'
    )
    exclude_unclassified_uncultured <-
        read.table(exclude_unclassified_uncultured_fname, header = FALSE)[[1]] |>
        as.character()

    ## 7. Exclude unclassified and informal names
    exclude_unclassified_informal_fname <- system.file(
        'extdata/proc_exclude_unclassified_informal.txt', package = 'taxPPro'
    )
    exclude_unclassified_informal <-
        read.table(exclude_unclassified_informal_fname, header = FALSE)[[1]] |>
        as.character()

    ## 8. Exclude uncultured and informal names
    exclude_uncultured_informal_fname <- system.file(
        'extdata/proc_exclude_uncultured_informal.txt', package = 'taxPPro'
    )
    exclude_uncultured_informal <-
        read.table(exclude_uncultured_informal_fname, header = FALSE)[[1]] |>
        as.character()

    list(
        all_ids = all_ids,
        exclude_all = exclude_all,
        exclude_unclassified = exclude_unclassified,
        exclude_uncultured = exclude_uncultured,
        exclude_informal = exclude_informal,
        exclude_unclassified_uncultured = exclude_unclassified_uncultured,
        exclude_unclassified_informal = exclude_unclassified_informal,
        exclude_uncultured_informal = exclude_uncultured_informal
    )
}

