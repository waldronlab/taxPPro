library(taxPPro)
library(dplyr)
library(bugphyzz)

gt <- physiologies('growth temperature')[[1]]
cg <- physiologies('coding genes')[[1]]


myFun <- function(dat) {
    dat <- dat |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown')
    return(dat)


    datMin <- dat
    datMin$Attribute_value_max <- NULL
    datMax <- dat
    datMax$Attribute_value_min <- NULL



    dat |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
            .data$Confidence_in_curation, n = 1, with_ties = TRUE
        ) |>
        dplyr::mutate(
            .data$Attribute_value_min = min(.data$Attribute_value_min),
            .data$Attribute_value_min_source = m,
            .data$Attribute_value_max = max(.data$Attribute_value_max),
            .data$Attribute_source = head(sort(.data$Attribute_source), 1)
        )
}

x <- split(gt, gt$Attribute_source)
y <- map(x, ~ unique(.x$NCBI_ID))


y$BacMap %in% y$`Microbial Fatty Acid Compositions`
dups <- gt$NCBI_ID[which(duplicated(gt$NCBI_ID))]
View(gt[gt$NCBI_ID %in% dups,])


dupscg <- cg$NCBI_ID[which(duplicated(cg$NCBI_ID))]
dupscgdf <- cg[cg$NCBI_ID %in% dupscg,]

dupscgdf |>
    filter(NCBI_ID == '96345') |>
    View()


cg2 <- physiologies('coding genes')
