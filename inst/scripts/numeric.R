library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data

# gt <- physiologies('growth temperature')[[1]]
# cg <- physiologies('coding genes')[[1]]
w <- physiologies('width')[[1]]

myFun <- function(dat) {
    dat <- dat |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max))
    select_cols <- c(
        'NCBI_ID', 'Attribute', 'Attribute_group',
        'Attribute_value_min', 'Attribute_value_max',
        'Attribute_source', 'Confidence_in_curation',
        'Evidence', 'Frequency', 'Unit'
    )

    change_cols <- c(
        'Attribute_source', 'Confidence_in_curation',
        'Evidence', 'Frequency'
    )

    dat <- dat |>
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

    return(dat)


    datMin <- dplyr::select(dat, tidyselect::any_of(select_cols))
    datMin$Attribute_value_max <- NULL
    datMax <- dplyr::select(dat, tidyselect::any_of(select_cols))
    datMax$Attribute_value_min <- NULL

   datMin <- datMin |>
       dplyr::group_by(.data$NCBI_ID) |>
       dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = TRUE
       ) |>
       dplyr::mutate(
           Attribute_value_min = min(.data$Attribute_value_min, na.rm = TRUE),
       ) |>
       dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
       dplyr::arrange(.data$NCBI_ID) |>
       dplyr::ungroup()

   pos_min <- which(colnames(datMin) %in% change_cols)
   colnames(datMin)[pos_min] <- paste0(colnames(datMin)[pos_min], '_min')

   datMax <- datMax |>
       dplyr::group_by(.data$NCBI_ID) |>
       dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = TRUE
       ) |>
       dplyr::mutate(
           Attribute_value_max = max(.data$Attribute_value_max, na.rm = TRUE),
       ) |>
       dplyr::slice_max(.data$Attribute_source, n = 1, with_ties = FALSE) |>
       dplyr::arrange(.data$NCBI_ID) |>
       dplyr::ungroup()

   pos_max <- which(colnames(datMax) %in% change_cols)
   colnames(datMin)[pos_max] <- paste0(colnames(datMax)[pos_max], '_max')

   list(min = datMin, max = datMax)
}

## Percentage of differing range values was too low (less than 2%)
## I just got the mean value

dat <- myFun(w)
any(duplicated(dat$NCBI_ID))

mean(dat$NCBI_ID %in% tip_data$taxid) * 100
mean(tip_data$taxid %in% dat$NCBI_ID) * 100

annotated_tips <- left_join(tip_data, dat, by = c('taxid' = 'NCBI_ID'))

tip_values <- annotated_tips$Attribute_value
names(tip_values) <- annotated_tips$tip_label

fit <- fastAnc(tree = tree, x = tip_values, vars=TRUE, CI=TRUE)















