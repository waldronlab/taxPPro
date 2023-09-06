
library(ape)
library(phytools)
library(tidyr)
library(purrr)
library(dplyr)
library(ggtree)
library(plotly)

getMRCA <- function(t, df) {
    res <- phytools::findMRCA(t, tips = df[['tip_label']])
    if (is.null(res))
        res <- NA
    res
}
col_names <- c(
    'genome_id', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
    'species'
)
taxRanks <- c(
    'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
)

tfile <- 'http://cmprod1.cibio.unitn.it/biobakery3/mpa_v31_CHOCOPhlAn_201901_species_tree.nwk'
t <- read.tree(tfile)
df <- data.frame(tip_label = t$tip.label)
df <- df |>
    separate(col = 'tip_label', into = col_names, sep = '\\|', remove = FALSE) |>
    modify(~ sub('^\\w__', '', .x))

mrcas <- map(taxRanks, ~ {
     splitted_df <- split(df, factor(df[[.x]]))
     output <- map_chr(splitted_df, \(y) as.character(getMRCA(t = t, df = y)))
     return(output)
})

v <- unlist(mrcas)
v <- v[which(!is.na(v))]
v <- v[!duplicated(v)]

v2 <- names(v)
names(v2) <- v

nodes <- as.character(length(t$tip.label) + 1:t$Nnode)


v3 <- v2[nodes]
names(v3) <- nodes
v3[is.na(v3)] <- ''


t$node.label <- v3

# p <- ggtree(t, branch.length = 'none', size = 0.025)
p <- ggtree(t, layout = 'circular', branch.length = 'none', size = 0.025)

p <- p +
    geom_nodelab(size = 1)
ggsave(
    filename = 'glab.png', p, width = 10, height = 10, dpi = 600
)

# toWebGL(p)
