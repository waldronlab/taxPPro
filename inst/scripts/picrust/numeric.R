
library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)
library(castor)
library(cvTools)
library(purrr)
library(tidyr)
library(ggplot2)

# n <- getMRCA( phy = tree, tip = c('g__620', 'g__620'))
# all_labels <- c(tree$tip.label, tree$node.label)
# all_labels[n]

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
tip_data_annotated <- left_join(tip_data, dat, by = c('taxid' = 'NCBI_ID')) |>
    select(tip_label, Attribute_value)

allTips <- tip_data_annotated$Attribute_value
names(allTips) <- tip_data_annotated$tip_label

values <- allTips[!is.na(allTips)]

## Create train and test sets ####
set.seed(1234)
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

    x <- testSets[[i]]
    x[!is.na(x)] <- NA
    y <- trainSets[[i]]
    z = c(x, y)
    a = allTips[!names(allTips) %in% names(z)]
    input_vector <- c(z, a)
    input_vector <- input_vector[tree$tip.label]

    res <- hsp_squared_change_parsimony(
        tree = tree, tip_states = input_vector, weighted = TRUE,
        check_input = TRUE
    )

    statesDF <- data.frame(
        label = tree$tip.label,
        value = res$states[1:Ntip(tree)]
    )

    states <- statesDF$value
    names(states) <- statesDF$label

    # commonNames <- intersect(names(states), names(testSets[[i]]))
    # hsp[[i]] <- states[commonNames]
    hsp[[i]] <- states[names(testSets[[i]])]
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
        # Metric = c("MSE", "RMSE", "MAE", "R_squared", "MAPE"),
        Metric = c("MSE", "RMSE", "R_squared"),
        Value = c(mse, rmse, r_squared)
        # Value = c(mse, rmse, mae, r_squared, mape)
    )

    return(evaluation_results)
}) |>
    bind_rows() |>
    group_by(Metric) |>
    summarise(
        mean = mean(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE)
    )

plotList <- map2(hsp, testSets, ~ {
    df <- data.frame(pred = .x, actual = .y)
    df |>
        ggplot(aes(pred, actual)) +
        geom_point()
})

p <- ggpubr::ggarrange(plotlist = plotList, ncol = 5, nrow = 2)


write.table(x = metrics, 'inst/scripts/castor_gt.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

ggsave('inst/scripts/castor_gt.png', width = 12, height = 4)








