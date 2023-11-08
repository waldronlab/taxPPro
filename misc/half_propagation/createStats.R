library(dplyr)
library(taxPPro)
library(purrr)

dat <- read.table(
    file = 'misc/half_propagation/bugphyzz_export_2023-11-07.tsv',
    header = TRUE, sep = '\t'
) |>
    filter(!is.na(Attribute_group) & !is.na(Attribute)) |>
    mutate(
        data_name = case_when(
            Attribute_type %in% c('multisate-union') ~ paste0(Attribute_group, '|', sub('--.*$', '', Attribute)),
            TRUE ~ Attribute_group
        )

    )

l <- split(dat, dat$data_name)

ltp <- ltp()
tip_data <- ltp$tip_data
sourceCodes <- c('exp', 'igc', 'nas', 'tas')

output <- map(l, ~ {
    source_data <- .x |>
        filter(Evidence %in% sourceCodes)
    before <- mean(tip_data$NCBI_ID %in% unique(source_data$NCBI_ID)) * 100
    after <- mean(tip_data$NCBI_ID %in% unique(dat$NCBI_ID)) * 100
    data.frame(
        before = before,
        after = after
    )
}) |>
    bind_rows(.id = 'data_name')
