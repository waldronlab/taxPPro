library(dplyr)
library(taxPPro)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(data.tree)

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
    file = 'misc/half_propagation/habitat.tsv',
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
    plot = p1, width = 10, height = 9, dpi = 150, bg = 'white'
)


tree <- ltp$tree
node_data <- ltp$node_data
tip_data <- ltp$tip_data
data('tree_list')
ncbi_tree <- as.Node(tree_list)
ncbi_nodes <- ncbi_tree$Get(
    attribute = 'name',
    filterFun = function(node) node$name != 'ArcBac',
    simplify = TRUE
) |>
    unname()

ranks_var <- c('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom')
ranks_ncbi <- count(data.frame(r = sub('__\\d+$', '', ncbi_nodes)), r, name = 'NCBI') |>
    mutate(
        Rank = case_when(
            r == 'c' ~ 'class',
            r == 'f' ~ 'family',
            r == 'g' ~ 'genus',
            r == 'k' ~ 'kingdom',
            r == 'o' ~ 'order',
            r == 'p' ~ 'phylum',
            r == 's' ~ 'species',
            r == 't' ~ 'strain'
        )
    ) |>
    mutate(Rank = factor(Rank, levels = ranks_var)) |>
    arrange(Rank) |>
    select(-r)


ranks <- bind_rows(select(node_data, Rank), select(tip_data, Rank))
rownames(ranks) <- NULL
ranks_ltp <- ranks |>
    mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
    filter(Rank %in% ranks_var) |>
    mutate(Rank = factor(Rank, levels = ranks_var)) |>
    count(Rank, name = 'LTP') |>
    arrange(Rank)

rank_data <- full_join(ranks_ltp, ranks_ncbi, by = 'Rank')
rank_data_long <- rank_data |>
    pivot_longer(
        names_to = 'tree', values_to = 'n', cols = c(LTP, NCBI)
    )

rank_data_long |>
    ggplot(aes(Rank, n)) +
    geom_col(aes(fill = tree), position = 'dodge') +
    facet_wrap(.~Rank, scales = 'free', ncol = 2) +
    # scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
    labs(
        y = 'Number of taxids'
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )
