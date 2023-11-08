library(dplyr)
library(taxPPro)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggpubr)

ltp <- ltp()
tip_data <- ltp$tip_data
sourceCodes <- c('exp', 'igc', 'nas', 'tas')

# Physiologies ------------------------------------------------------------
p <- read.table(
    file = 'misc/half_propagation/no_habitat.tsv',
    header = TRUE, sep = '\t'
) |>
    filter(!is.na(Attribute_group) & !is.na(Attribute)) |>
    mutate(
        data_name = case_when(
            Attribute_type %in% c('multistate-union') ~ paste0(Attribute_group, '|', sub('--.*$', '', Attribute)),
            TRUE ~ Attribute_group
        )

    )

p_l <- split(p, p$data_name)


p_data <- map(p_l, ~ {
    source_data <- .x |>
        filter(Evidence %in% sourceCodes)
    before <- mean(tip_data$NCBI_ID %in% unique(source_data$NCBI_ID)) * 100
    after <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
    data.frame(
        before = before,
        after = after
    )
}) |>
    bind_rows(.id = 'data_name')

(
    p_p <- p_data |>
        arrange(-after) |>
        mutate(data_name = forcats::fct_inorder(data_name)) |>
        pivot_longer(
            names_to = 'time', values_to = 'per', cols = c(before, after)
        ) |>
        mutate(time = factor(time, levels = c('before', 'after'))) |>
        ggplot(mapping = aes(data_name, per)) +
        geom_col(mapping = aes(fill = time), position = 'dodge') +
        geom_hline(yintercept = 10, linetype = 'dotdash') +
        scale_y_continuous(
            limits = c(0, 100), labels = function(x) paste0(x, "%")
        ) +
        scale_fill_manual(
            name = 'Pooling', labels = c('Before', 'After'),
            values = c('gray60', 'gray80')
        ) +
        labs(
            x = 'Physiology/Attribute',
            y = 'LTP tips annotated'
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        )
)

# habitat -----------------------------------------------------------------

h <- read.table(
    file = 'misc/half_propagation/bugphyzz_export_2023-11-08.tsv',
    header = TRUE, sep = '\t'
) |>
    filter(!is.na(Attribute_group) & !is.na(Attribute)) |>
    mutate(
        data_name = case_when(
            Attribute_type %in% c('multistate-union') ~ paste0(Attribute_group, '|', sub('--.*$', '', Attribute)),
            TRUE ~ Attribute_group
        )

    )

h_l <- split(h, h$data_name)


h_data <- map(h_l, ~ {
    source_data <- .x |>
        filter(Evidence %in% sourceCodes)
    before <- mean(tip_data$NCBI_ID %in% unique(source_data$NCBI_ID)) * 100
    after <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
    data.frame(
        before = before,
        after = after
    )
}) |>
    bind_rows(.id = 'data_name')

(
    h_p <- h_data |>
        arrange(-after) |>
        mutate(data_name = sub('^habitat\\|', '', data_name)) |>
        mutate(data_name = forcats::fct_inorder(data_name)) |>
        head(nrow(p_data)) |>
        pivot_longer(
            names_to = 'time', values_to = 'per', cols = c(before, after)
        ) |>
        mutate(time = factor(time, levels = c('before', 'after'))) |>
        ggplot(mapping = aes(data_name, per)) +
        geom_col(mapping = aes(fill = time), position = 'dodge') +
        geom_hline(yintercept = 10, linetype = 'dotdash') +
        scale_y_continuous(
            limits = c(0, 100), labels = function(x) paste0(x, "%")
        ) +
        scale_fill_manual(
            name = 'Pooling', labels = c('Before', 'After'),
            values = c('gray60', 'gray80')
        ) +
        labs(
            x = 'Attribute',
            y = 'LTP tips annotated'
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        )
)

p1 <- ggarrange(
    p_p, h_p, nrow = 2, common.legend = TRUE, align = 'hv', legend = 'right',
    labels = c('A)', 'B)'), hjust = 0

)

ggsave(
    filename = 'misc/half_propagation/percent_plot.png',
    plot = p1, width = 10, height = 9, dpi = 100, bg = 'white'
)
