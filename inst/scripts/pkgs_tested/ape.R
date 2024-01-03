library(ape)
## For some reason, this ML option is really fast (maybe I should use this one instead)
## I would need to map taxonomic levels for internal nodes
ace_result <- ape::ace(
    x = x, phy = phy, type = "discrete",
    method = "ML", model = "ER", CI = TRUE,
    marginal = TRUE, kappa = 0
)

pruned_tree <- keep.tip(phy = mpa_tree, tip = rownames(m1))
ace_result <- ape::ace(
    x = states, phy = pruned_tree, type = "discrete",
    method = "ML", model = "ER", CI = TRUE,
    marginal = TRUE, kappa = 0
)


states <- vector('character', nrow(m1))
for (i in seq_along(states)) {
    states[i] <- paste0(sort(names(which(m1[i,] > 0))), collapse = ":")
}
names(states) <- rownames(m1)
states[which(grepl(':', states))] <- NA
