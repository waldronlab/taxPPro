library(bugphyzz)
library(taxPPro)
library(dplyr)
library(castor)
library(tibble)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)


# Functions ---------------------------------------------------------------
removeASR <- function(x) {
    x |>
        dplyr::filter(Evidence != "asr") |>
        mutate(NCBI_ID = as.character(NCBI_ID))
}
addAttributes <- function(attrDat, tipDat) {
    dplyr::left_join(
        tipDat, attrDat,
        by = c("taxid" = "NCBI_ID"), relationship = "many-to-many"
    ) |>
        filter(!is.na(Attribute_value))
}

removeDups <- function(x) {
    x |>
        dplyr::rowwise() |>
        dplyr::mutate(
            sorted = paste(sort(c(node1, node2)), collapse = "_")
        ) |>
        dplyr::ungroup() |>
        dplyr::filter(!duplicated(sorted)) |>
        dplyr::select(-sorted)
}
getClosestTips <- function(tr) {
    purrr::map(seq_along(tr$tip.label), function(i) {
        if ((i %% 5000) == 0)
            message(i)
        res <- castor::find_nearest_tips(tr, target_tips = tr$tip.label[i])
        tip_distances <- res[["nearest_distance_per_tip"]]
        names(tip_distances) <- tr$tip.label
        tip_distances <- tip_distances[-i]
        minValue <- min(tip_distances)
        tip_distances <- tip_distances[which(tip_distances == minValue)]
        data.frame(
            node2 = names(tip_distances),
            distance = unname(tip_distances)
        )
    }) |>
        purrr::set_names(tr$tip.label) |>
        dplyr::bind_rows(.id = "node1") |>
        removeDups()
}
# -------------------------------------------------------------------------

ltp <- ltp()
tr <- ltp$tree
tip_data <- ltp$tip_data |>
    group_by(taxid) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n == 1) |>
    select(-n)
bp <- importBugphyzz()

gt <- bp$`optimal ph`

gt2 <- gt |>
    removeASR() |>
    select(NCBI_ID, Attribute_value)

gt_tip_data <- left_join(tip_data, gt2, by = c("taxid" = "NCBI_ID")) |>
    filter(!is.na(Attribute_value))




# gtDat <- gt_tip_data$Attribute_value
# names(gtDat) <- gt_tip_data$tip_label


# discrete_trait_depth()
#
#
# if (any(gtDat < 0)) {
#     negs <- which(gtDat < 0)
#     message("Found ", length(negs), " negatives. Dropping them.")
#     gtDat <- gtDat[-negs]
# }
# gtDat <- log(gtDat + 1)
#
#
#
# tim <- system.time({
#     res <- phylosig(tr, x = gtDat, method = "K", test = TRUE)
# })



# aerDat <- bp$aerophilicity |>
#     removeASR() |>
#     addAttributes(tipDat = tip_data)


annotatedTips <- unique(gt_tip_data$tip_label)
closeTips <- getClosestTips(tr)


x <- select(gt_tip_data, tip_label, Attribute_value)

y <- left_join(closeTips, x, by = c("node1"  = "tip_label")) |>
    rename(node1Val = Attribute_value)
z <- left_join(y, x, by = c("node2"  = "tip_label")) |>
    rename(node2Val = Attribute_value)
z <- drop_na(z)
z$diff <- abs(z$node1Val - z$node2Val)


z |>
    ggplot(aes(log(distance + 1), diff)) +
    geom_point()

res <- lm(formula = distance ~ diff, data = z)
summary(res)

cor(x = z$distance, z$diff, method = "spearman")



annotationsL <- split(aerDat, aerDat$Attribute_value) |>
    map(~ pull(.x, tip_label))

for (i in seq_along(annotationsL)) {
    annName <- names(annotationsL)[i]

    colName1 <- paste0(annName, "_node1")
    colName2 <- paste0(annName, "_node2")

    closeTips[[colName1]] <- closeTips$node1 %in% annotationsL[[i]]
    closeTips[[colName2]] <- closeTips$node2 %in% annotationsL[[i]]

    closeTips[[annName]] <- paste0(closeTips[[colName1]], "|", closeTips[[colName2]])
}


dat <- closeTips |>
    select(-matches("_node\\d$")) |>
    pivot_longer(
        names_to = "Attribute", values_to = "logical",
        cols = aerobic:last_col()
    )


dat2 <- dat |>
    mutate(
        logical = case_when(
            grepl("FALSE", logical) ~ "FALSE",
            TRUE ~ logical
        )
    ) |>
    filter(Attribute != "aerotolerant")

dat2 |>
    count(Attribute, logical) |>
    View()

dat2 |>
    filter(logical != FALSE) |>
    ggplot(aes(x = Attribute, y = distance)) +
    geom_boxplot(aes(fill = logical))
