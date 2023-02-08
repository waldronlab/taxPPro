library(dplyr)
library(bugphyzz)
library(tidygraph)
library(igraph)
library(taxPPro)
library(data.tree)

## Define some functions ####

## Add children
removePrefix <- function(x) {
    regex <- '^[kpcofgst]__'
    sub(regex, '', x)
}
addChildren <- function(node) {
    current_children <- removePrefix(names(node$children))
    current_taxon <- removePrefix(node$name)
    df <- tryCatch(
        error = function(e) NULL, {
          taxizedb::children(current_taxon)[[1]]
        }
    )
    if (is.null(df)) {
        return(NULL)
    }
    df <- df[df$rank %in% valid_ranks, ]
    children <- df$name
    new_children <- children[!children %in% current_children]
    if (isFALSE(!length(new_children))) {
        for (i in seq_along(new_children)) {
            message(new_children)
            node$AddChild(new_children[i])
        }
    }
}


## Some code ####
valid_ranks <- validRanks()
ncbi_taxonomy <- get_ncbi_taxonomy() |>
    filter(Rank %in% valid_ranks) |>
    select(-Parent_NCBI_ID, -Rank, -Taxon_name)

aer <- physiologies('aerophilicity', remove_false = TRUE)[[1]] |>
    mutate(NCBI_ID = tolower(as.character(NCBI_ID))) |>
    filter(
        Rank %in% valid_ranks,
        !is.na(NCBI_ID),
        NCBI_ID != 'unknown',
        !is.na(Attribute_value)
    ) |>
    preSteps(tax.id.type = 'NCBI_ID')

bac_data <- left_join(aer, ncbi_taxonomy, by = 'NCBI_ID') |>
    filter(kingdom == 'Bacteria')

bac_data$pathString <- paste(
    'k__', bac_data$kingdom,  ## This will be the root
    '|||p__', bac_data$phylum,
    '|||c__', bac_data$class,
    '|||o__', bac_data$order,
    '|||f__', bac_data$family,
    '|||g__', bac_data$genus,
    '|||s__', bac_data$species,
    '|||t__', bac_data$strain,
    sep = ''
)

regex1 <- '(\\|\\|\\|[kpcofgst]__NA)*$'
bac_data$pathString <- sub(regex1, "", bac_data$pathString)
regex2 <- '(\\|\\|\\|[kpcofgst]__[^\\|]*$)'
bac_data$Taxon_name <- stringr::str_extract(bac_data$pathString, regex2) |>
    stringr::str_remove('\\|\\|\\|')

bac_data_tree <- bac_data |>
    select(
        kingdom, phylum, class, order, family, genus, species, strain,
        Taxon_name, pathString
    ) |>
    distinct()

es <- bac_data_tree[grep('Escherichia', bac_data_tree$pathString),]
es <- es[1,]

start_time <- Sys.time()
bac_tree <- as.Node(bac_data_tree, pathDelimiter = '|||' )
end_time <- Sys.time()
(time_big_tree <- difftime(end_time, start_time))
## Stop here


## function for adding full path (taxonomy) to the names
## This is so that every node has a unique name
# addPath <- function(node) {
    # node$name <- paste0(node$name, '---', node$pathString)
# }
# bac_tree$Do(addPath)
# bac_tree$totalCount
# bac_tree$leafCount

## This is for a test with igraph
# new_tbl <- as_tbl_graph(bac_tree, directed = TRUE)


## Let's try with a smaller dataset

es <- bac_data_tree[grep('Escherichia', bac_data_tree$pathString),]
es <- es[1,]

start_time <- Sys.time()
es_tree <- as.Node(es, pathDelimiter = '|||' )
end_time <- Sys.time()
time_es_tree <- difftime(end_time, start_time)


bac_tree$Do(addChildren)













## Add bugphyzz data
split_by_taxname <- split(bac_data, factor(bac_data$Taxon_name))
addBugphyzzAttributes <- function(node, list_of_df) {
    select_cols <- c('NCBI_ID', 'Attribute', 'Evidence', 'Score', 'Frequency')
    node_name <- node$name
    if (node_name %in% names(list_of_df)) {
        message(node_name, ' has attributes in bugphyzz')
        node$bugphyzz_data <- list_of_df[[node_name]][,select_cols]
    }
}

# # tree$Do(addBugphyzzAttributes, list_of_df = split_by_taxname)
#
# ## Escherichia (genus)
# length(names(new_tree$p__Proteobacteria$c__Gammaproteobacteria$o__Enterobacterales$f__Enterobacteriaceae$g__Escherichia$children))
# ## Lactobacillus sekei (species)
# bac_tree$p__Firmicutes$c__Bacilli$o__Lactobacillales$f__Lactobacillaceae$g__Latilactobacillus$`s__Latilactobacillus sakei`$`t__Latilactobacillus sakei subsp. sakei 23K`
#
# View(print(tree$p__Proteobacteria$c__Gammaproteobacteria$o__Enterobacterales$f__Enterobacteriaceae$g__Escherichia, 'NCBI_ID', 'Attribute', 'Attribute_value', 'Attribute_source'))
#
# level_names <- data |>
#     filter(grepl('\\[Actinobacillus\\] rossii', Taxon_name)) |>
#     pull(pathString) |>
#     stringr::str_split(pattern = '\\|\\|\\|') |>
#     purrr::flatten_chr() |>
#     {\(y) y[-length(y)]}() |>
#     {\(y) y[-length(y)]}() |>
#     {\(y) y[-1]}()
# tree$Climb(
#     name = level_names
# )$children[1:3]
#
# taxizedb::children('Escherichia', db = 'ncbi', verbose = FALSE)
#
# # getChildren2 <- function(node) {
# #     taxon <- sub('^[kpcofgst]__', '', node$name)
# #     node$new_children <- tryCatch(
# #         error = function(e) NULL, {
# #             taxizedb::children(taxon, db = 'ncbi')
# #         }
# #     )
# # }
# #
# # tree$Do(getChildren2)
# #
# # taxizedb::children('Bacteria', db = 'ncbi')
#
# # Quick test with ete3r ---------------------------------------------------
#
# ncbi_taxonomy$helper_pathString <- paste(
#     'k__', ncbi_taxonomy$kingdom,  ## This will be the root
#     '|||p__', ncbi_taxonomy$phylum,
#     '|||c__', ncbi_taxonomy$class,
#     '|||o__', ncbi_taxonomy$order,
#     '|||f__', ncbi_taxonomy$family,
#     '|||g__', ncbi_taxonomy$genus,
#     '|||s__', ncbi_taxonomy$species,
#     '|||t__', ncbi_taxonomy$strain,
#     sep = ''
# )
#
# tax <- unique(bac_data$NCBI_ID)
# output_tree <- getTree(taxids = tax)
# data1 <- output_tree@data
# new_names <- unique(data1$name[data1$name != ""])
#
# colnames(ncbi_taxonomy)
#
# x <- ncbi_taxonomy |>
#     filter(NCBI_ID %in% new_names)
# x$pathString <- x$helper_pathString
#
# new_tree <- as.Node(x, pathDelimiter = '|||')
nodes <- sort(unique(c(aer$NCBI_ID, aer$Parent_NCBI_ID)))
node_df <- data.frame(nodes = nodes)

edges <- distinct(aer[,c('Parent_NCBI_ID', 'NCBI_ID')])
colnames(edges) <- c('from', 'to')
edges <- data.frame(

    from = aer$NCBI_ID
)

my_graph <- graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)
tg <- as_tbl_graph(my_graph, directed = TRUE)

