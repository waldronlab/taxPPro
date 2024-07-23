library(bugphyzz)
library(taxPPro)
library(dplyr)
library(castor)
library(tibble)
library(readr)

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

# -------------------------------------------------------------------------

ltp <- ltp()
tr <- ltp$tree
tip_data <- ltp$tip_data

bp <- importBugphyzz()

aerDat <- bp$aerophilicity |>
    filter(Attribute_value == "aerobic") |>
    removeASR() |>
    addAttributes(tipDat = tip_data)
annotatedTips <- unique(aerDat$tip_label)
pwd <- get_all_pairwise_distances(tree = tr, only_clades = annotatedTips)

colnames(pwd) <- annotatedTips
rownames(pwd) <- annotatedTips

indices <- which(upper.tri(pwd, diag = FALSE), arr.ind = TRUE)

distances <- pwd[indices]
node1 <- rownames(pwd)[indices[, 1]] # rows
node2 <- colnames(pwd)[indices[, 2]] # columns

distDF <- data.frame(
    node1 = node1,
    node2 = node2,
    distance = distances
)
