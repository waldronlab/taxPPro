
library(bugphyzz)
library(purrr)
library(taxPPro)
library(dplyr)
library(purrr)
library(ComplexHeatmap)

phys <- physiologies()
phys <- map(phys, as_tibble)
remove_dat <- c(
    'habitat', 'antimicrobial resistance', 'isolation site', 'width',
    'length', 'disease association',
    'animal pathogen', 'biofilm forming', 'growth temperature'
)

phys <- phys[!names(phys) %in% remove_dat]

phys_prop <- map(phys, ~ propagate(.x, 'genus', 'both'))

proc_exclude_all_fname <-
    system.file('extdata/proc_exclude_all.txt', package = 'taxPPro')
proc_exclude_all <- read.table(proc_exclude_all_fname, header = FALSE)[[1]] |>
    as.character()



phys_with_taxids <- phys |>
    map(~ filter(.x, !is.na(NCBI_ID) | NCBI_ID != 'unknown')) |>
    map(~ filter(.x, Rank %in% c('genus', 'species', 'strain'))) |>
    map(~ preSteps(.x, 'NCBI_ID')) |>
    discard(~ !nrow(.x))

phys_prop_with_taxids <- phys_prop |>
    map(~ filter(.x, !is.na(NCBI_ID) | NCBI_ID != 'unknown')) |>
    map(~ filter(.x, Rank %in% c('genus', 'species', 'strain'))) |>
    discard(~ !nrow(.x))



x <- phys_with_taxids$aerophilicity

x <- split(x, factor(x$Rank))

sum(x[[2]]$NCBI_ID %in% proc_exclude_all) / length(proc_exclude_all)


map_dbl(x, ~ sum(.x$NCBI_ID %in% proc_exclude_all) / length(proc_exclude_all) * 100)


calcCompleteness <- function(df) {
    split_df <- split(df, factor(df$Rank))
    data <- map_dbl(
        split_df, ~ {
            sum(.x$NCBI_ID %in% proc_exclude_all) / length(proc_exclude_all) * 100
        }
    )
    output <- matrix(data, nrow = 1)
    colnames(output) <- names(split_df)
    as.data.frame(output)
}

my_ranks <- c('genus', 'species', 'strain')
complete1 <- map(phys_with_taxids, calcCompleteness) |>
    bind_rows(.id = 'Attribute') |>
    tibble::column_to_rownames(var = 'Attribute')
complete1[is.na(complete1)] <- 0
complete1 <- complete1[,my_ranks]


complete2 <- map(phys_prop_with_taxids, calcCompleteness) |>
    bind_rows(.id = 'Attribute') |>
    tibble::column_to_rownames(var = 'Attribute')
complete2[is.na(complete2)] <- 0
complete2 <- complete2[, my_ranks]


row_names <- rownames(complete1)

hp1 <- Heatmap(
    matrix = as.matrix(complete1),
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_names_side = 'bottom', row_names_side = 'left',
    cluster_columns = FALSE, cluster_rows = FALSE,
    name = '%',
    col = circlize::colorRamp2(c(0, 1, 100), c("white", 'yellow', "red"))
)
png(filename = 'complete1.png')
hp1
dev.off()

complete2 <- complete2[row_names,]

hp2 <- Heatmap(
    matrix = as.matrix(complete2),
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_names_side = 'bottom', row_names_side = 'left',
    cluster_columns = FALSE, cluster_rows = FALSE,
    name = '%',
    col = circlize::colorRamp2(c(0, 1, 100), c("white", 'yellow', "red"))
)
png(filename = 'complete2.png')
hp2
dev.off()


## Divide main ncbi set by rank.
## Consider remove strain column if it is empty.


