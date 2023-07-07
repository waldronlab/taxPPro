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

aer <- aer |>
    group_by(NCBI_ID) |>
    mutate(Score = Score / sum(Score)) |>
    ungroup()


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

tree$d__2$p__1224$c__1236$o__91347$f__543$g__561$s__562$table

tree$Do(function(node) {
    node[['table']] <- NULL
}

)

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
main_table <- dplyr::bind_rows(tables)
    # dplyr::filter(Evidence %in% c('asr', 'inh'))

x <- tidyr::complete(
    main_table, NCBI_ID, Attribute, fill = list(Score = 0, Evidence = NA)
    ) |>
    relocate(NCBI_ID)

end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time

## Plotting a tree or some branches ####

library(tidytree)
library(ggtree)
library(TreeTools)
library(purrr)

data('tree_sp')
sub_trees <- ape::subtrees(tree_sp)
fam_tree <- keep(sub_trees, ~ .x$node.label[1] == 'f__543')[[1]]
# fam_tree <- keep(sub_trees, ~ .x$node.label[1] == 'g__561')[[1]]
new_tree <- RRphylo::fix.poly(fam_tree, type = 'resolve')

# new_tree <- fam_tree

stats <- x |>
    select(-Evidence) |>
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score'
    ) |>
    filter(NCBI_ID %in% fam_tree$node.label) |>
    # filter(!grepl('^[st]__', NCBI_ID)) |>
    rename(node = NCBI_ID)
    # tibble::column_to_rownames(var = 'NCBI_ID') |>
    # as.matrix()


palette('Tableau 10')
attrs <- colnames(stats)[2:10]
cols <- setNames(palette()[1:length(unique(attrs))],sort(unique(attrs)))
pies <- nodepie(stats, cols = 2:10)
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols))
p <- ggtree(new_tree, branch.length = 'none') +
    geom_tiplab() +
    # geom_tippoint(aes(color = stat)) +
    scale_color_manual(values = cols) +
    theme(legend.position = "right")
    # xlim(NA, 8)
p2 <- p + ggtree::geom_inset(pies, width = 0.1, height = 0.1)


# ggsave(filename = 'test.png', p2)
# fam_tree <- keep(sub_trees, ~ .x$node.label[1] == 'f__543')[[1]]
# fam_tree <- keep(sub_trees, ~ .x$node.label[1] == 'f__543')[[1]]

