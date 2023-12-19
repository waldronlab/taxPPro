library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)
library(castor)
library(cvTools)
library(purrr)
library(tidyr)
library(ggplot2)




n <- getMRCA( phy = tree, tip = c('g__620', 'g__620'))
all_labels <- c(tree$tip.label, tree$node.label)
all_labels[n]


## Tree data ####
ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data
gn_tips <- ltp$gn_tips
label_data <- bind_rows(
    distinct(select(as_tibble(tip_data), label = tip_label, taxid)),
    distinct(select(as_tibble(node_data), label = node_label, taxid))
)

## bugphyzz data ####
gt <- physiologies('growth temperature')[[1]]
# cg <- physiologies('coding genes')[[1]]
# w <- physiologies('width')[[1]]

modifyNumeric <- function(x) {
    x |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = FALSE
        ) |>
        dplyr::mutate(
            Attribute_value = mean(.data$Attribute_value_min, .data$Attribute_value_max)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.data$Attribute_value_min, -.data$Attribute_value_max) |>
        dplyr::distinct()
}

dat <- modifyNumeric(gt)

## Small checks and loggin
any(duplicated(dat$NCBI_ID))
mean(dat$NCBI_ID %in% tip_data$taxid) * 100
mean(tip_data$taxid %in% dat$NCBI_ID) * 100

dat <- dat[dat$NCBI_ID %in% tip_data$taxid,]

values <- dat$Attribute_value
names(values) <- dat$NCBI_ID

## Create train and test sets ####
foldN <- cvFolds(n = length(values), K = 10)
testSets <- vector('list', 10)
trainSets <- vector('list', 10)
for (i in 1:10) {
    foldName <- paste0('Fold', i)
    testSets[[i]] <- values[foldN$subsets[foldN$which == i]]
    names(testSets)[i] <- foldName
    trainSets[[i]] <- values[foldN$subsets[foldN$which != i]]
    names(trainSets)[i] <- foldName
}

## hidden-state-prediction ####
hsp <- vector('list', 10)
for (i in seq_along(trainSets)) {
    names(hsp)[i] <- names(trainSets)[i]
    dat_subset <- dat[dat$NCBI_ID %in% names(trainSets[[i]]),]
    annotated_tips <- left_join(
        tip_data, dat_subset, by = c('taxid' = 'NCBI_ID')
    )
    tip_values <- annotated_tips$Attribute_value
    names(tip_values) <- annotated_tips$tip_label
    tip_values <- tip_values[tree$tip.label]
    # res <- hsp_squared_change_parsimony(
    #     tree = tree, tip_states = tip_values, weighted = TRUE,
    #     check_input = TRUE
    # )
    res <- hsp_independent_contrasts(
        tree = tree, tip_states = tip_values, weighted = TRUE,
        check_input = TRUE
    )

    statesDF <- data.frame(
        label = c(tree$tip.label, tree$node.label),
        value = res$states
    ) |>
        filter(label != 'NA') |> # NAs were introduced as character strings
        filter(!grepl('^n\\d+$', label)) |>
        left_join(label_data, by = 'label') |>
        filter(!label %in% ltp$gn_tips) # remove genus tips to avoid duplicated, they are already in the internal nodes

    states <- statesDF$value
    names(states) <- statesDF$taxid

    commonNames <- intersect(names(states), names(testSets[[i]]))
    hsp[[i]] <- states[commonNames]
}

## Calculate metrics
metrics <- map2(hsp, testSets, ~ {
    predicted_values <- .x
    actual_values <- .y

    mse <- mean((predicted_values - actual_values)^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(predicted_values - actual_values))

    ss_total <- sum((actual_values - mean(actual_values))^2)
    ss_residual <- sum((actual_values - predicted_values)^2)
    r_squared <- 1 - (ss_residual / ss_total)

    mape <- mean(abs((actual_values - predicted_values) / actual_values)) * 100

    evaluation_results <- data.frame(
        Metric = c("MSE", "RMSE", "MAE", "R_squared", "MAPE"),
        Value = c(mse, rmse, mae, r_squared, mape)
    )

    return(evaluation_results)
}) |>
    bind_rows(.id = 'Fold') |>
    pivot_wider(
        names_from = 'Metric', values_from = 'Value'
    )

plotList <- map2(hsp, testSets, ~ {
    df <- data.frame(pred = .x, actual = .y)
    df |>
        ggplot(aes(pred, actual)) +
        geom_point()
})

ggpubr::ggarrange(plotlist = plotList, ncol = 2, nrow = 5)

tbls <- map2(hsp, testSets, ~ {
    data.frame(hsp = .x, test = .y)
})
