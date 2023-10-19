library(taxPPro)
library(data.tree)
library(ape)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data

data('tree_list')
data_tree <- as.Node(tree_list)
ncbi_tree_nodes <- data_tree$Get(
    attribute = 'name', filterFun = function(node) {
        grepl('[st]__', node$name)
    }
)
ncbi_tree_nodes <- unname(ncbi_tree_nodes)

tip_data2 <- tip_data |>
    filter(NCBI_ID %in% ncbi_tree_nodes)

tim <- system.time({
    res <- cophenetic(tree)
})

df <- res |>
    as.data.frame() |>
    rownames_to_column(var = 'tip1') |>
    as_tibble() |>
    pivot_longer(
        cols = 2:last_col(), names_to = 'tip2', values_to = 'distance'
    ) |>
    filter(tip2 %in% tip_data2$tip_label) |> ## filter only tips that can be mapped to the NCBI tree
    group_by(tip1) |>
    slice_max(order_by = distance, n = 1) ## allow ties

df2 <- left_join(df, tip_data,by = c('tip2' = 'tip_label'))

write_tsv(
    x = df2, file = 'inst/extdata/longest_distance_between_tips.tsv'
)

#     res2 <- dist.nodes(tree)
# })


