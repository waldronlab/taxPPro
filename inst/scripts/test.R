library(ape)
library(phytools)
library(dplyr)
library(tidyr)

data('primate.tree')
data('primate.data')
tree <- primate.tree
data <- primate.data
original <- data |>
    tibble::rownames_to_column(var = 'Taxa') |>
    select(Taxa, Activity_pattern) |>
    mutate(Presence = 1) |>
    pivot_wider(
        names_from = 'Activity_pattern', values_from = 'Presence',
        values_fill = 0
    ) |>
    arrange(Taxa) |>
    tibble::column_to_rownames(var = 'Taxa') |>
    as.matrix()

myFun <- function(mat, per = 0.1) {
    n_row <- nrow(mat)
    n <- round(n_row * per)
    rows <- sample(x = 1:nrow(mat), size = n, replace = FALSE)
    mat[rows,] <- rep(1/ncol(mat), ncol(mat))
    mat
}

input_mat <- myFun(original, 0.95)
input_mat <- input_mat[tree$tip.label,]
fit <- fitMk(tree = primate.tree, x = input_mat,
              model = "ARD", pi = "fitzjohn",
              lik.func = "pruning", logscale = TRUE)
ace <- ancr(fit, tips=TRUE)
plot(ace, args.plotTree = list(direction="upwards"))
