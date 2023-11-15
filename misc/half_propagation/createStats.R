library(bugphyzz)
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
node_data <- ltp$node_data
sourceCodes <- c('exp', 'igc', 'nas', 'tas')

b <- importBugphyzz() |>
    mutate(
        NCBI_ID = case_when(
            Rank == 'kingdom' ~ paste0('k__', NCBI_ID),
            Rank == 'phylum' ~ paste0('p__', NCBI_ID),
            Rank == 'class' ~ paste0('c__', NCBI_ID),
            Rank == 'order' ~ paste0('o__', NCBI_ID),
            Rank == 'family' ~ paste0('f__', NCBI_ID),
            Rank == 'genus' ~ paste0('g__', NCBI_ID),
            Rank == 'species' ~ paste0('s__', NCBI_ID),
            Rank == 'strain' ~ paste0('t__', NCBI_ID)
        )
    )

source_l <- b |>
    filter(Evidence %in% sourceCodes) |>
    {\(y) split(y, factor(y$Attribute_group))}()
taxPool_l <- b |>
    filter(Evidence %in% c('tax', 'inh')) |>
    {\(y) split(y, factor(y$Attribute_group))}()

# all(names(source_l) == names(taxPool_l))

dat <- map2(
    .x = source_l, .y = taxPool_l, .f =  ~ {
        before <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
        after <- mean(tip_data$NCBI_ID %in% unique(.y$NCBI_ID)) * 100
        data.frame(before = before, after = after)
    }
) |>
    bind_rows(.id = 'data_name')

p1 <- dat |>
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
            name = 'taxPool/inh', labels = c('Before', 'After'),
            values = c('gray60', 'gray80')
        ) +
        labs(
            x = 'Physiology/Attribute_group',
            y = 'LTP tips annotated'
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1)
        )

ggsave(
    filename = 'misc/half_propagation/percent_plot.png',
    plot = p1, width = 10, height = 5, dpi = 150, bg = 'white'
)


## Compare completeness between LTP and NCBI tree



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
    labs(
        y = 'Number of taxids'
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )
