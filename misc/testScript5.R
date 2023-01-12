
library(bugphyzz)
library(taxPPro)
library(dplyr)
library(data.tree)
library(purrr)

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

data <- left_join(aer, ncbi_taxonomy, by = 'NCBI_ID') |>
    filter(kingdom == 'Bacteria')

data <- ncbi_taxonomy
data$pathString <- paste(
    'k__', data$kingdom,  ## This will be the root
    '|||p__', data$phylum,
    '|||c__', data$class,
    '|||o__', data$order,
    '|||f__', data$family,
    '|||g__', data$genus,
    '|||s__', data$species,
    '|||t__', data$strain,
    sep = ''
)
regex1 <- '(\\|\\|\\|[kpcofgst]__NA)*$'
data$pathString <- sub(regex1, "", data$pathString)
regex2 <- '(\\|\\|\\|[kpcofgst]__[^\\|]*$)'
data$Taxon_name <- stringr::str_extract(data$pathString, regex2) |>
    stringr::str_remove('\\|\\|\\|')
data_tree <- data |>
    select(
        kingdom, phylum, class, order, family, genus, species, strain,
        Taxon_name, pathString
    ) |>
    distinct()

bac_data <- data_tree[data_tree$kingdom == 'Bacteria', ]
arc_data <- data_tree[data_tree$kingdom == 'Archaea',]

start_time <- Sys.time()
arc_tree <- as.Node(arc_data, pathDelimiter = '|||' )
end_time <- Sys.time()
difftime(end_time, start_time)

start_time <- Sys.time()
bac_tree <- as.Node(bac_data, pathDelimiter = '|||' )
end_time <- Sys.time()
difftime(end_time, start_time)




tree <- as.Node(data, pathDelimiter = '|||')

## Add bugphyzz data
split_by_taxname <- split(data, factor(data$Taxon_name))
addBugphyzzAttributes <- function(node, list_of_df) {
    select_cols <- c('NCBI_ID', 'Attribute', 'Evidence', 'Score', 'Frequency')
    node_name <- node$name
    if (node_name %in% names(list_of_df)) {
        message(node_name, ' has attributes in bugphyzz')
        node$bugphyzz_data <- list_of_df[[node_name]][,select_cols]
    }
}

tree$Do(addBugphyzzAttributes, list_of_df = split_by_taxname)

## Escherichia (genus)
tree$p__Proteobacteria$c__Gammaproteobacteria$o__Enterobacterales$f__Enterobacteriaceae$g__Escherichia
## Lactobacillus sekei (species)
tree$p__Firmicutes$c__Bacilli$o__Lactobacillales$f__Lactobacillaceae$g__Latilactobacillus$`s__Latilactobacillus sakei`$`t__Latilactobacillus sakei subsp. sakei 23K`

View(print(tree$p__Proteobacteria$c__Gammaproteobacteria$o__Enterobacterales$f__Enterobacteriaceae$g__Escherichia, 'NCBI_ID', 'Attribute', 'Attribute_value', 'Attribute_source'))

level_names <- data |>
    filter(grepl('\\[Actinobacillus\\] rossii', Taxon_name)) |>
    pull(pathString) |>
    stringr::str_split(pattern = '\\|\\|\\|') |>
    purrr::flatten_chr() |>
    {\(y) y[-length(y)]}() |>
    {\(y) y[-length(y)]}() |>
    {\(y) y[-1]}()
tree$Climb(
    name = level_names
)$children[1:3]



taxizedb::children('Escherichia', db = 'ncbi', verbose = FALSE)

# getChildren2 <- function(node) {
#     taxon <- sub('^[kpcofgst]__', '', node$name)
#     node$new_children <- tryCatch(
#         error = function(e) NULL, {
#             taxizedb::children(taxon, db = 'ncbi')
#         }
#     )
# }
#
# tree$Do(getChildren2)
#
# taxizedb::children('Bacteria', db = 'ncbi')

