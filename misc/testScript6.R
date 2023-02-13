library(taxPPro)
library(bugphyzz)
library(data.tree)
library(purrr)

## Define function
fillNAs <- function(x) {
    counter <- 1
    for (i in seq_along(x)) {
        if (is.na(x[i])) {
            x[i] <- paste0('NA', counter)
            counter <- counter + 1
        }
    }
    return(x)
}

## Get NCBI data
ncbi_taxonomy <- get_ncbi_taxonomy()
ncbi_taxids <- get_ncbi_taxids(keyword = 'b', with_taxids = TRUE)
lgl_vct <- ncbi_taxids$Rank == 'species' & ncbi_taxids$kingdom == 'Bacteria'
bacteria <- ncbi_taxids[lgl_vct,]
no_select_cols <- c(
    'strain', 'Taxon_name', 'species', 'Parent_NCBI_ID', 'Rank'
)
bacteria <- bacteria[,!colnames(bacteria) %in% no_select_cols]

## Convert to tree
# pathString <- paste(
#     'k__', bacteria$kingdom, # root
#     '|||p__', bacteria$phylum,
#     '|||c__', bacteria$class,
#     '|||o__', bacteria$order,
#     '|||f__', bacteria$family,
#     '|||g__', bacteria$genus,
#     '|||s__', bacteria$NCBI_ID,
#     sep = ''
# )

bacteria$NCBI_ID <- paste0('s__', bacteria$NCBI_ID)
pathString <- paste(
    bacteria$kingdom, # root
    '|||', bacteria$phylum,
    '|||', bacteria$class,
    '|||', bacteria$order,
    '|||', bacteria$family,
    '|||', bacteria$genus,
    '|||', bacteria$NCBI_ID,
    sep = ''
)

# pathString <- sub(regex1, '', pathString)
bacteria$pathString <- pathString
bacteria <- modify(bacteria, fillNAs)

a <- table(bacteria$phylum) |> as.data.frame() |> View()

##
# phyl <- bacteria[bacteria$phylum == 'Bacillota',]
start_time <- Sys.time()
bac_tree <- as.Node(bacteria, pathDelimiter = '|||')
end_time <- Sys.time()
(elapsed_time <- difftime(end_time, start_time))

# bac_tree$Do(addChildren)




