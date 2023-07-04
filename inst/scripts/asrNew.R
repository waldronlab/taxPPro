library(data.tree)
library(bugphyzz)
library(dplyr)
library(taxPPro)
library(tidyr)

data("tree_list")
tree <- as.Node(tree_list)

aer <- physiologies('aerophilicity', remove_false = TRUE, full_source = FALSE)[[1]] |>
    # mutate(
    #     letter = sub('^(.).*', '\\1', Rank),
    #     NCBI_ID = paste0(letter, '__', NCBI_ID)
    # ) |>
    prepareDatForPropagation()
l <- split(aer, factor(aer$NCBI_ID))
l[['g__561']] <- NULL

start_time <- Sys.time()

tree$Do(function(node) {
    if (!is.null(l[[node$name]])) {
        node[['table']] <- l[[node$name]] |>
            dplyr::select(NCBI_ID, Attribute, Score, Evidence) |>
            dplyr::distinct() |>
            tibble::as_tibble()
    }
})

calcParent <- function(tbl) {
    tbl |>
        dplyr::group_by(NCBI_ID) |>
        dplyr::mutate(Score =  Score / sum(Score)) |>
        dplyr::mutate(Score = ifelse(is.na(Score), 0, Score)) |>
        dplyr::group_by(Attribute) |>
        dplyr::reframe(
            Score = sum(Score)
            # Evidence = paste0(Evidence, collapse = '|'),
        ) |>
        dplyr::mutate(Score = Score / sum(Score)) |>
        dplyr::mutate(Score = ifelse(is.na(Score), 0, Score))
        # dplyr::mutate(Evidence = sub('^(\\|*)', '', Evidence)) |>
        # dpyr::mutate(Evidence = sub('\\|\\|+', '|', Evidence))
}

myFun <- function(node) {
    if (is.null(node[['table']])) {
        # message(node$name)
        children <- names(node$children)
        output <- vector('list', length(children))
        for (i in seq_along(output)) {
            output[[i]] <- node[[children[i]]]$table
        }
        cond <- purrr::map_lgl(output, is.null)
        if (!all(cond)) {
            df <- purrr::discard(output, is.null) |>
                dplyr::bind_rows() |>
                dplyr::select(
                    NCBI_ID, Attribute, Score, # Evidence #, Attribute_source
                ) |>
                calcParent() |>
                dplyr::mutate(
                    NCBI_ID = node$name,
                    Evidence = 'asr'
                ) |>
                dplyr::relocate(NCBI_ID)
            node[['table']] <- df
        }
    }
}

myFun2 <- function(node) {
    if (is.null(node[['table']])) {
        df <- node[['parent']][['table']]
        df <- df |>
            dplyr::mutate(
                NCBI_ID = node$name,
                Evidence = 'inh'
            )
        node[['table']] <- df
    }
}

tree$Do(myFun, traversal = 'post-order')
tree$Do(myFun2, traversal = 'pre-order')
tables <- tree$Get(function(node) node[['table']], simplify = FALSE) |>
    purrr::discard(~ all(is.na(.x)))
main_table <- dplyr::bind_rows(tables) |>
    dplyr::filter(Evidence %in% c('asr', 'inh'))

x <- tidyr::complete(
    main_table, NCBI_ID, Attribute, fill = list(Score = 0, Evidence = NA)
    ) |>
    relocate(NCBI_ID)

end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time





# after -------------------------------------------------------------------

mat <- x |>
    select(-Evidence) |>
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score'
    ) |>
    filter(!grepl('^[st]__', NCBI_ID)) |>
    tibble::column_to_rownames(var = 'NCBI_ID') |>
    as.matrix()

data('tree_sp')
