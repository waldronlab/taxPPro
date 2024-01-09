library(taxPPro)
library(bugphyzz)
library(dplyr)
library(phytools)
library(purrr)

ltp <- ltp()
full_tree <- ltp$tree

full_tree$node.label <- ifelse(
    full_tree$node.label == "NA",
    paste0("n", Ntip(full_tree) + 1:Nnode(full_tree)),
    full_tree$node.label
)


tip_data <- ltp$tip_data
phys <- physiologies('coding genes')[[1]]

modifyNumeric <- function(x) {
    x |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::mutate(Rank = taxizedb::taxid2rank(NCBI_ID, db = 'ncbi')) |>
        dplyr::mutate(Taxon_name = taxizedb::taxid2name(NCBI_ID, db = 'ncbi')) |>
        dplyr::mutate(NCBI_ID = addRankPrefix(NCBI_ID, Rank)) |>
        dplyr::filter(!is.na(Rank) & !is.na(Taxon_name) & !is.na(NCBI_ID)) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
            .data$Confidence_in_curation, n = 1, with_ties = FALSE
        ) |>
        dplyr::mutate(
            Attribute_value = mean(.data$Attribute_value_min, .data$Attribute_value_max)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-Attribute_value_min, -Attribute_value_max) |>
        dplyr::distinct()
}

dat <- modifyNumeric(phys)

fdat <- left_join(
    x = tip_data,
    y = select(dat, NCBI_ID, Attribute_value),
    by = "NCBI_ID"
)

trait_table <- fdat |>
    select(tip_label, Attribute_value) |>
    rename(taxa = tip_label, trait1 = Attribute_value) |>
    filter(!is.na(trait1)) |>
    tibble::column_to_rownames(var = 'taxa')
pruned_tree <- keep.tip(full_tree, trait_table$taxa)
data_ordered <- trait_table[pruned_tree$tip.label,, drop = FALSE]
states <- data_ordered$trait1
names(states) <- rownames(data_ordered)

asr <- fastAnc(tree = pruned_tree, x = states, vars = TRUE, CI = TRUE)

res <- data.frame(
    node = pruned_tree$node.label,
    asr = unname(asr$ace),
    var = unname(asr$var),
    CI_lower = unname(asr$CI95[,1, drop = TRUE]),
    CI_upper = unname(asr$CI95[,2, drop = TRUE])

)






