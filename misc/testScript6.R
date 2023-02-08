library(bugphyzz)
library(taxPPro)
library(ete3r)
library(data.tree)

valid_ranks <- validRanks()
ncbi_taxonomy <- get_ncbi_taxonomy() |>
    filter(Rank %in% valid_ranks) |>
    select(-Parent_NCBI_ID, -Rank, -Taxon_name) |>
    mutate(
        pathString = paste0(
            'k__', kingdom,
            '|||p__', phylum,
            '|||c__', class,
            '|||o__', order,
            '|||f__', family,
            '|||g__', genus,
            '|||s__', species,
            '|||t__', strain,
            sep = ''
        )
    )

regex1 <- '(\\|\\|\\|[kpcofgst]__NA)*$'
ncbi_taxonomy$pathString <- sub(regex1, "", ncbi_taxonomy$pathString)
regex2 <- '(\\|\\|\\|[kpcofgst]__[^\\|]*$)'
ncbi_taxonomy$Taxon_name <- stringr::str_extract(ncbi_taxonomy$pathString, regex2) |>
    stringr::str_remove('\\|\\|\\|')

aer <- physiologies('aerophilicity', remove_false = TRUE)[[1]] |>
    mutate(NCBI_ID = tolower(as.character(NCBI_ID))) |>
    filter(
        Rank %in% valid_ranks,
        !is.na(NCBI_ID),
        NCBI_ID != 'unknown',
        !is.na(Attribute_value)
    ) |>
    preSteps(tax.id.type = 'NCBI_ID')

## In this step I get all taxids with full taxonomy in the ncbi taxonomy
bac_taxids <- ncbi_taxonomy |>
    filter(
        NCBI_ID %in% unique(aer$NCBI_ID),
        kingdom == 'Bacteria'
    ) |>
    pull(NCBI_ID) |>
    unique()

phylo_tree <- getTree(bac_taxids)
tax_data <- phylo_tree@data

new_bac_taxids <- unique(tax_data$name[tax_data$name != ''])

bac_data <- ncbi_taxonomy |>
    filter(NCBI_ID %in% tax_names)

bac_tree <- as.Node(bac_data, pathDelimiter = '|||')

names(bac_tree$p__Proteobacteria$c__Gammaproteobacteria$o__Enterobacterales$f__Enterobacteriaceae$g__Escherichia$children)

test_tree <- getTree(c('561', '562', '2762229'))
print(test_tree)


