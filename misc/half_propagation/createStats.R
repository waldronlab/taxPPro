library(bugphyzz)
library(dplyr)
library(taxPPro)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(data.tree)

ltp <- ltp3()
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


# all physiologies --------------------------------------------------------
source_l <- b |>
    filter(Evidence %in% sourceCodes) |>
    {\(y) split(y, factor(y$Attribute_group))}()
taxPool_l <- b |>
    filter(Evidence %in% c(sourceCodes, c('tax', 'inh'))) |>
    {\(y) split(y, factor(y$Attribute_group))}()

# all(names(source_l) == names(taxPool_l))

source_l <- source_l[names(taxPool_l)]
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
        geom_hline(yintercept = 5, linetype = 'dotdash') +
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


# Habitat -----------------------------------------------------------------
h1 <- split(source_l$habitat, factor(source_l$habitat$Attribute))
h2 <- split(taxPool_l$habitat, factor(taxPool_l$habitat$Attribute))
# all(names(h) == names(h2))
dat_h <- map2(
    .x = h1, .y = h2, .f =  ~ {
        before <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
        after <- mean(tip_data$NCBI_ID %in% unique(.y$NCBI_ID)) * 100
        data.frame(before = before, after = after)
    }
) |>
    bind_rows(.id = 'data_name')

p2 <- dat_h |>
    arrange(-after) |>
    mutate(data_name = forcats::fct_inorder(data_name)) |>
    pivot_longer(
        names_to = 'time', values_to = 'per', cols = c(before, after)
    ) |>
    head(30) |>
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
    filename = 'misc/half_propagation/percent_plot_habitat.png',
    plot = p2, width = 9, height = 4, dpi = 150, bg = 'white'
)

# antimicrobial resistance ------------------------------------------------

ar1 <- split(source_l$`antimicrobial resistance`, factor(source_l$`antimicrobial resistance`$Attribute))
ar2 <- split(taxPool_l$`antimicrobial resistance`, factor(taxPool_l$`antimicrobial resistance`$Attribute))

ar3 <- ar1[which(names(ar1) %in% names(ar2))]
ar4 <- ar2[names(ar3)]

dat_ar <- map2(
    .x = ar3, .y = ar4, .f =  ~ {
        before <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
        after <- mean(tip_data$NCBI_ID %in% unique(.y$NCBI_ID)) * 100
        data.frame(before = before, after = after)
    }
) |>
    bind_rows(.id = 'data_name')

p3 <- dat_ar |>
    arrange(-after) |>
    mutate(data_name = forcats::fct_inorder(data_name)) |>
    pivot_longer(
        names_to = 'time', values_to = 'per', cols = c(before, after)
    ) |>
    head(30) |>
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
    filename = 'misc/half_propagation/percent_plot_antimicrobial_resistance.png',
    plot = p3, width = 9, height = 4, dpi = 150, bg = 'white'
)

ar5 <- ar2[!names(ar2) %in% names(ar1)]

dat_ar2 <- map(ar5, ~ {
    after <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
    data.frame(after = after)
}) |>
    bind_rows(.id = 'data_name')

p4 <- dat_ar2 |>
    arrange(-after) |>
    mutate(data_name = forcats::fct_inorder(data_name)) |>
    pivot_longer(
        names_to = 'time', values_to = 'per', cols = c(after)
    ) |>
    head(30) |>
    mutate(time = factor(time, levels = c('before', 'after'))) |>
    ggplot(mapping = aes(data_name, per)) +
    geom_col(mapping = aes(fill = time), position = 'dodge') +
    geom_hline(yintercept = 10, linetype = 'dotdash') +
    scale_y_continuous(
        limits = c(0, 100), labels = function(x) paste0(x, "%")
    ) +
    scale_fill_manual(
        name = 'taxPool/inh', labels = c('After'),
        values = c('gray80')
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
    filename = 'misc/half_propagation/percent_plot_antimicrobial_resistance_2.png',
    plot = p4, width = 9, height = 4, dpi = 150, bg = 'white'
)

# disease association -----------------------------------------------------

da1 <- split(source_l$`disease association`, factor(source_l$`disease association`$Attribute))
da2 <- split(taxPool_l$`disease association`, factor(taxPool_l$`disease association`$Attribute))
names_in_common <- intersect(names(da1), names(da2))
da1 <- da1[names_in_common]
da2 <- da2[names_in_common]

dat_da <- map2(
    .x = da1, .y = da2, .f =  ~ {
        before <- mean(tip_data$NCBI_ID %in% unique(.x$NCBI_ID)) * 100
        after <- mean(tip_data$NCBI_ID %in% unique(.y$NCBI_ID)) * 100
        data.frame(before = before, after = after)
    }
) |>
    bind_rows(.id = 'data_name')

p5 <- dat_da |>
    arrange(-after) |>
    mutate(data_name = forcats::fct_inorder(data_name)) |>
    pivot_longer(
        names_to = 'time', values_to = 'per', cols = c(before, after)
    ) |>
    head(30) |>
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
    filename = 'misc/half_propagation/percent_plot_disease_association.png',
    plot = p5, width = 9, height = 4, dpi = 150, bg = 'white'
)

# Completeness of ltp and ncbi trees --------------------------------------

## Compare completeness between LTP and NCBI tree
data('tree_list')
ncbi_tree <- as.Node(tree_list)
ncbi_nodes <- ncbi_tree$Get(
    attribute = 'Rank',
    filterFun = function(node) node$name != 'ArcBac',
    simplify = TRUE
)

ranks_var <- c('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom')

ncbi_df <- data.frame(rank = ncbi_nodes, taxid = names(ncbi_nodes)) |>
    mutate(taxid = sub('\\w__', '', taxid))
rownames(ncbi_df) <- NULL
ranks_ncbi <- ncbi_df |>
    count(rank, name = 'NCBI') |>
    mutate(rank = factor(rank, levels = ranks_var)) |>
    arrange(rank)


ltp_df <- bind_rows(
    select(node_data, rank = Rank, taxid),
    select(tip_data, rank = Rank, taxid),
) |>
    tibble::remove_rownames() |>
    mutate(rank = ifelse(rank == 'superkingdom', 'kingdom', rank)) |>
    filter(rank %in% ranks_var)

ranks_ltp <- ltp_df |>
    count(rank, name = 'LTP') |>
    mutate(rank = factor(rank, levels = ranks_var)) |>
    arrange(rank)

output <- vector('integer', length(ranks_var))
for (i in seq_along(ranks_var)) {
    v1 <- ncbi_df[which(ncbi_df$rank == ranks_var[i]),]$taxid
    v2 <- ltp_df[which(ltp_df$rank == ranks_var[i]),]$taxid
    names(output)[i] <- ranks_var[i]
    output[[i]] <- length(intersect(v1, v2))
}

inter <- data.frame(
    rank = names(output),
    BOTH = output
) |>
    tibble::remove_rownames()

rank_data <- reduce(list(ranks_ltp, ranks_ncbi, inter), full_join) |>
    mutate(rank = factor(rank, levels = ranks_var))
rank_data_long <- rank_data |>
    pivot_longer(
        names_to = 'tree', values_to = 'n', cols = c(LTP, NCBI,  BOTH)
    )

p5 <- rank_data_long |>
    ggplot(aes(rank, n)) +
    geom_col(aes(fill = tree), position = 'dodge') +
    facet_wrap(~rank, ncol = 4, scales = 'free') +
    scale_fill_manual(
        name = '',
        values = c('gray75', 'gray50', 'gray25'),
        labels = c('both', 'ltp', 'ncbi')
    ) +
    scale_y_continuous(
        breaks = scales::pretty_breaks()
    ) +
    labs(
        y = 'Number of taxids'
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave(
    filename = 'misc/half_propagation/ranks_numbers.png',
    plot = p5, width = 7, height = 3, dpi = 150, bg = 'white'
)

