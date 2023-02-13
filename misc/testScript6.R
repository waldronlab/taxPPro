library(taxPPro)
library(bugphyzz)
library(data.tree)
library(purrr)

## Define functions ####
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

taxname2taxid <- function(tax_tbl) {

    if ('kingdom' %in% colnames(tax_tbl)) {
        pos <- which(colnames(tax_tbl) == 'kingdom')
        colnames(tax_tbl)[pos] <- 'superkingdom'
    }

    ranks <- c(
        'superkingdom', 'phylum', 'class', 'order', 'family', 'genus'
    )
    for (i in ranks) {
        df <- tax_tbl[tax_tbl[['Rank']] == i, ]
        vct <- df$NCBI_ID
        names(vct) <- df[[i]]
        tax_tbl[[i]] <- purrr::map_chr(tax_tbl[[i]], ~ {vct[.x]})
    }
    return(tax_tbl)
}

## Get NCBI data ####
ncbi_taxids <- get_ncbi_taxids(keyword = 'b', with_taxids = TRUE)
new_ncbi_taxids <- taxname2taxid(tax_tbl = ncbi_taxids)
cond1 <- new_ncbi_taxids$Rank == 'species'
cond2 <- new_ncbi_taxids$superkingdom == '2'
bacteria <- new_ncbi_taxids[cond1 & cond2,]
no_select_cols <- c(
    'species', 'strain', 'Taxon_name', 'Parent_NCBI_ID', 'Rank'
)
bacteria <- bacteria[,!colnames(bacteria) %in% no_select_cols]

pathString <- paste(
    'k__', bacteria$superkingdom, # root
    '|||p__', bacteria$phylum,
    '|||c__', bacteria$class,
    '|||o__', bacteria$order,
    '|||f__', bacteria$family,
    '|||g__', bacteria$genus,
    '|||s__', bacteria$NCBI_ID,
    sep = ''
)

bacteria$pathString <- pathString
bacteria <- modify(bacteria, fillNAs)

## Convert full bacteria data to a data.tree object
start_time <- Sys.time()
bac_tree <- as.Node(bacteria, pathDelimiter = '|||')
end_time <- Sys.time()
(elapsed_time <- difftime(end_time, start_time))

print(bac_tree, limit = 20)

bac_tree$p__1224$c__1236$o__91347$f__543$g__561

e <- bacteria[bacteria$genus == '561',]
e_tree <- as.Node(e, pathDelimiter = '|||')


e_tree$Do(addStrain)

e_tree$p__1224$c__1236$o__91347$f__543$g__561$s__562




